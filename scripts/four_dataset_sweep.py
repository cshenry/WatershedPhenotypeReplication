#!/usr/bin/env python
"""H100 JOB 3 -- four-dataset Biolog phenotype sweep over the .MMGF models.

Simulates the 324-media Biolog panel for every genome across the four datasets the
deliverable is organised around, and scores each against Clayton's 2022 calls and the
real K-12 experimental calls. Produces the per-media call data the dashboard's
`count (% changed vs Original 2022p)` tables are built from.

    python scripts/four_dataset_sweep.py --workers 50
    python scripts/four_dataset_sweep.py --limit 2 --workers 1     # smoke test FIRST

The four datasets (Chris's labels -- used verbatim as dataset keys, lowercased):

    original_2022p   Clayton's OWN 2022 spreadsheet calls. The baseline. NOT simulated
                     -- read straight from notebooks/data/BiologMacTestOutput_On-Off.csv,
                     column "<gid>.fbamodel". Every other dataset's "% changed" is
                     measured against this.
    resim_2022p      our resimulation of the 2022 models  (265353/<gid>_2022p.MMGF)
    resim_2025k      the 2025 KBase re-annotated models    (265353/<gid>_2025k.MMGF)
    resim_2026p      the 2026 PATRIC rebuild               (265353/<gid>_2026p.MMGF)

Why this is NOT three_arm_sweep.py (the retired script): the three .MMGF sets are
ALREADY gapfilled (H100 job 2), with the leaky-directionality fix baked in before the
gapfill. So this sweep does NOT gapfill and does NOT rebuild -- it loads each .MMGF
model and simulates it directly. Re-gapfilling would change the very models being
compared.

The simulation is EXACTLY the path that maximally replicated Clayton, unchanged:
  * media via get_biolog_media(block_co2=True) -- CO2 uptake closed, so the latent
    carbon-fixation route cannot inflate growth calls,
  * simulate_biolog_panel(add_missing_exchanges=True) -- the 2022 run's third
    positional arg to simulate_phenotypes,
  * growth called at biomass flux > GROWTH_THRESHOLD (0.01).

Saves the ACTUAL per-media flux/call for every (genome, dataset, media) to --calls,
never just the scorecard -- so the per-substrate tables, the growth matrices for Aaron,
and any re-scoring at a different threshold are all answerable without re-simulating.

Design mirrors three_arm_sweep.py's proven infrastructure: one task per GENOME (all
four datasets together, so the Original baseline and the three sims are in hand at once
to compute "% changed"); per-worker cache of the shared read-only inputs; init and
per-fetch retry against transient KBase blips; checkpointed + resumable by
(genome, dataset).
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

warnings.filterwarnings("ignore")

_PRJ = Path(__file__).resolve().parent.parent / "notebooks" / "PRJ-watershed_phenotype_replication"
sys.path.insert(0, str(_PRJ))

TARGET_WS = 265353

# (dataset_key, model_suffix). model_suffix None == the Clayton baseline, not simulated.
DATASETS = [
    ("original_2022p", None),
    ("resim_2022p", "_2022p"),
    ("resim_2025k", "_2025k"),
    ("resim_2026p", "_2026p"),
]
DATASET_KEYS = [k for k, _ in DATASETS]

SUMMARY_FIELDS = [
    "genome_id", "dataset", "model_name", "n_reactions",
    "n_scored",              # media with a call for this genome x dataset
    "n_growth",              # growth-positive media -- the "count" in the tables
    "n_changed_vs_orig",     # media whose call differs from Original 2022p
    "pct_changed_vs_orig",   # 100 * n_changed / n_comparable
    "free_lunch",            # biomass from minerals alone (model-quality flag; sim only)
    "clayton_accuracy", "clayton_CP", "clayton_CN", "clayton_FP", "clayton_FN",
    "exp_accuracy", "exp_CP", "exp_CN", "exp_FP", "exp_FN",
    "status", "error",
]

CALL_FIELDS = [
    "genome_id", "dataset", "media_name", "panel",
    "flux",             # simulated biomass flux; "" for original_2022p (not simulated)
    "growth",           # 0/1 call for this dataset
    "original_2022p",   # 0/1 Clayton baseline call, carried on every row for easy diffing
    "experimental",     # 0/1 real K-12 Biolog call ("" where unknown)
]

_CTX: dict = {}


def _kbase_retry(fn, tries=4, base_delay=2.0):
    """Call fn(), retrying transient KBase connection failures with backoff.

    A full sweep makes thousands of workspace calls; an occasional ConnectionError /
    ChunkedEncodingError is expected and should self-heal, not fail the job.
    """
    for attempt in range(tries):
        try:
            return fn()
        except Exception as exc:  # noqa: BLE001
            transient = any(t in type(exc).__name__ for t in
                            ("Connection", "ChunkedEncoding", "Timeout", "Protocol"))
            if not transient or attempt == tries - 1:
                raise
            time.sleep(base_delay * (2 ** attempt))


def _init(media_limit=0):
    """Load the shared, read-only inputs once per worker process.

    Runs inside the ProcessPoolExecutor INITIALIZER: if it raises, the whole pool
    breaks and every job fails without writing a row, so it is retried as a unit --
    not just the per-job fetches. A single SSL handshake blip at startup would
    otherwise nuke an entire resume.

    media_limit > 0 keeps only the first N media -- a SMOKE-TEST knob so the plumbing
    (model load, four-dataset orchestration, %-changed math, output format) can be
    exercised without the multi-minute full-324 simulation. Never set in production;
    the scored numbers are meaningless on a subset.
    """
    from util import (get_msfba, get_biolog_media, load_reference_calls,
                      get_experimental_calls)

    def _load():
        fba = get_msfba()
        media = get_biolog_media(fba=fba)                # 324 media, CO2 uptake blocked
        if media_limit:
            media = dict(list(media.items())[:media_limit])
        return {
            "fba": fba,
            "media": media,
            "clayton": load_reference_calls(),           # 324 x 519, cols "<gid>.fbamodel"
            "experiment": get_experimental_calls(fba=fba),
        }

    for attempt in range(5):
        try:
            _CTX.update(_load())
            return
        except Exception:  # noqa: BLE001 -- retry any startup failure as a unit
            if attempt == 4:
                raise
            time.sleep(3.0 * (2 ** attempt))


def _score(rows, calls, threshold):
    """CP/CN/FP/FN of simulated calls vs a reference {media_name: 0/1-or-bool}."""
    CP = CN = FP = FN = 0
    for row in rows:
        name = row["media_name"]
        if name not in calls:
            continue
        sim = row["flux"] > threshold
        obs = bool(calls[name])
        if sim and obs:
            CP += 1
        elif (not sim) and (not obs):
            CN += 1
        elif sim and not obs:
            FP += 1
        else:
            FN += 1
    n = CP + CN + FP + FN
    return {"CP": CP, "CN": CN, "FP": FP, "FN": FN,
            "accuracy": (CP + CN) / n if n else None}


def _clayton_column(genome_id):
    """Clayton's baseline calls for one genome as {media_name: bool}, or None.

    None when the genome has no column in Clayton's matrix -- the sims still run and
    are emitted; only the "% changed vs original" and the original_2022p rows are
    unavailable for that genome.
    """
    clayton = _CTX["clayton"]
    col = f"{genome_id}.fbamodel"
    if col not in clayton.columns:
        return None
    series = clayton[col]
    return {str(m): bool(series[m]) for m in series.index}


def _run_genome(job):
    """All four datasets for ONE genome. Returns (summary_rows, call_rows).

    Grouping the four datasets into one task means Clayton's baseline column and the
    three simulated arms are all in hand together, so "% changed vs Original 2022p" is
    computed here rather than reconstructed downstream.
    """
    genome_id, datasets = job
    orig = _clayton_column(genome_id)   # {media: bool} or None
    exp = _CTX["experiment"]

    summary_rows, call_rows = [], []
    for dataset in datasets:
        suffix = dict(DATASETS)[dataset]
        if suffix is None:
            s_row, c_rows = _run_original(genome_id, orig, exp)
        else:
            s_row, c_rows = _run_sim(genome_id, dataset, suffix, orig, exp)
        summary_rows.append(s_row)
        call_rows.extend(c_rows)
    return summary_rows, call_rows


def _run_original(genome_id, orig, exp):
    """The Original 2022p dataset -- Clayton's own calls, no simulation."""
    row = {f: "" for f in SUMMARY_FIELDS}
    row["genome_id"], row["dataset"] = genome_id, "original_2022p"
    row["model_name"] = f"{genome_id}.fbamodel"
    calls = []
    if orig is None:
        row["status"], row["error"] = "skipped", "genome absent from Clayton matrix"
        return row, calls

    n_growth = 0
    for media_name, obs in orig.items():
        n_growth += int(obs)
        calls.append({
            "genome_id": genome_id, "dataset": "original_2022p",
            "media_name": media_name, "panel": str(media_name).split("-", 1)[0],
            "flux": "", "growth": int(obs),
            "original_2022p": int(obs),
            "experimental": ("" if media_name not in exp else int(bool(exp[media_name]))),
        })
    row["n_scored"] = len(orig)
    row["n_growth"] = n_growth
    row["n_changed_vs_orig"] = 0            # it IS the reference
    row["pct_changed_vs_orig"] = 0.0
    row["status"] = "ok"
    return row, calls


def _run_sim(genome_id, dataset, suffix, orig, exp):
    """One simulated dataset: load the .MMGF model and simulate the Biolog panel.

    Never raises -- a dead model comes back as a failed row so the sweep continues.
    """
    from util import simulate_biolog_panel, free_lunch_growth, GROWTH_THRESHOLD

    row = {f: "" for f in SUMMARY_FIELDS}
    name = f"{genome_id}{suffix}.MMGF"
    row["genome_id"], row["dataset"], row["model_name"] = genome_id, dataset, name
    calls = []
    try:
        fba = _CTX["fba"]
        mdlutl = _kbase_retry(lambda: fba.get_model(f"{TARGET_WS}/{name}"))
        cobra_model = mdlutl.model if hasattr(mdlutl, "model") else mdlutl
        row["n_reactions"] = len(cobra_model.reactions)

        # Model-quality flag: biomass from minerals alone (no carbon). Recorded, not
        # excluded -- with CO2 blocked the growth calls are still sound; free lunch is a
        # separate defect the analysis can filter on.
        row["free_lunch"] = int(free_lunch_growth(mdlutl) > 1e-6)

        sim_rows = simulate_biolog_panel(mdlutl, _CTX["media"], fba=fba)

        n_growth = n_changed = n_comparable = 0
        for r in sim_rows:
            m = r["media_name"]
            grow = int(r["flux"] > GROWTH_THRESHOLD)
            n_growth += grow
            o = orig.get(m) if orig is not None else None
            if o is not None:
                n_comparable += 1
                n_changed += int(bool(grow) != bool(o))
            calls.append({
                "genome_id": genome_id, "dataset": dataset,
                "media_name": m, "panel": r["panel"],
                "flux": round(r["flux"], 6), "growth": grow,
                "original_2022p": ("" if o is None else int(o)),
                "experimental": ("" if m not in exp else int(bool(exp[m]))),
            })

        row["n_scored"] = len(sim_rows)
        row["n_growth"] = n_growth
        if orig is not None:
            row["n_changed_vs_orig"] = n_changed
            row["pct_changed_vs_orig"] = round(100.0 * n_changed / n_comparable, 3) if n_comparable else ""

        if orig is not None:
            s = _score(sim_rows, orig, GROWTH_THRESHOLD)
            row["clayton_accuracy"] = round(s["accuracy"], 6) if s["accuracy"] is not None else ""
            for k in ("CP", "CN", "FP", "FN"):
                row[f"clayton_{k}"] = s[k]

        e = _score(sim_rows, exp, GROWTH_THRESHOLD)
        row["exp_accuracy"] = round(e["accuracy"], 6) if e["accuracy"] is not None else ""
        for k in ("CP", "CN", "FP", "FN"):
            row[f"exp_{k}"] = e[k]

        row["status"] = "ok"
    except Exception as exc:  # noqa: BLE001 -- a dead model must not kill the sweep
        row["status"] = "failed"
        row["error"] = f"{type(exc).__name__}: {exc}"[:200]
        traceback.print_exc(file=sys.stderr)
    return row, calls


def _done(path):
    """(genome, dataset) pairs already in the output -- the sweep is resumable."""
    if not path.exists():
        return set()
    seen = set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            gi, di = header.index("genome_id"), header.index("dataset")
        except ValueError:
            return set()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(gi, di) and parts[header.index("status")] == "ok":
                seen.add((parts[gi], parts[di]))
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=str(_PRJ / "NBOutput/four_dataset_sweep.tsv"))
    ap.add_argument("--calls", default=None,
                    help="ALSO write the per-media results (flux/call per genome x "
                         "dataset x media). The actual prediction data -- always set it. "
                         "Defaults to <output stem>_calls.tsv.")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2),
                    help="Parallel genome-tasks. 50 on H100; keep low (1-2) on a laptop.")
    ap.add_argument("--datasets", default=",".join(DATASET_KEYS),
                    help="Comma-separated subset of: " + ",".join(DATASET_KEYS))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--media-limit", type=int, default=0,
                    help="SMOKE TEST ONLY: simulate just the first N media. Scored "
                         "numbers are meaningless on a subset; never use in production.")
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    calls_path = Path(args.calls) if args.calls else out.with_name(out.stem + "_calls.tsv")
    calls_path.parent.mkdir(parents=True, exist_ok=True)

    datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]
    bad = [d for d in datasets if d not in DATASET_KEYS]
    if bad:
        ap.error(f"unknown dataset(s): {bad}; choose from {DATASET_KEYS}")

    _init(args.media_limit)  # in the parent, to enumerate genomes from Clayton's matrix
    genomes = sorted(c[: -len(".fbamodel")] for c in _CTX["clayton"].columns
                     if str(c).endswith(".fbamodel"))
    if args.limit:
        genomes = genomes[: args.limit]

    already = _done(out)
    jobs = []
    for g in genomes:
        todo = [d for d in datasets if (g, d) not in already]
        if todo:
            jobs.append((g, todo))
    n_pairs = sum(len(d) for _g, d in jobs)
    print(f"{len(genomes)} genomes x {len(datasets)} datasets = {len(genomes)*len(datasets)} "
          f"pairs; {len(already)} already done; {n_pairs} pairs across {len(jobs)} "
          f"genome-tasks on {args.workers} workers", flush=True)
    if not jobs:
        print("nothing to do", flush=True)
        return

    new_out = not out.exists()
    new_calls = not calls_path.exists()
    with open(out, "a", buffering=1) as fh, open(calls_path, "a", buffering=1) as cfh:
        if new_out:
            fh.write("\t".join(SUMMARY_FIELDS) + "\n")
        if new_calls:
            cfh.write("\t".join(CALL_FIELDS) + "\n")
        done = 0
        with ProcessPoolExecutor(max_workers=args.workers, initializer=_init,
                                 initargs=(args.media_limit,)) as pool:
            futures = {pool.submit(_run_genome, j): j for j in jobs}
            for fut in as_completed(futures):
                summary_rows, call_rows = fut.result()
                for row in summary_rows:
                    fh.write("\t".join(str(row.get(f, "")) for f in SUMMARY_FIELDS) + "\n")
                for c in call_rows:
                    cfh.write("\t".join(str(c.get(f, "")) for f in CALL_FIELDS) + "\n")
                done += 1
                bad_rows = [r for r in summary_rows if r["status"] not in ("ok", "skipped")]
                if done % 10 == 0 or bad_rows:
                    gid = summary_rows[0]["genome_id"] if summary_rows else "?"
                    st = ",".join(f"{r['dataset']}={r['status']}" for r in summary_rows)
                    print(f"[{done}/{len(jobs)} genomes] {gid}: {st}", flush=True)
    print(f"done -- {done} genome-tasks -> {out}", flush=True)
    print(f"         per-media results -> {calls_path}", flush=True)


if __name__ == "__main__":
    main()
