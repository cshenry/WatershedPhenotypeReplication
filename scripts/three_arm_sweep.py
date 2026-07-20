#!/usr/bin/env python
"""Full three-arm sweep: 501 genomes x {old, new, rebuilt}, scored against experiment.

Run on H100 with ~50 workers:

    python scripts/three_arm_sweep.py \
        --output NBOutput/three_arm_sweep.tsv \
        --calls  NBOutput/three_arm_calls.tsv \
        --workers 50

Parallelised BY GENOME (all arms per task) so each genome is fetched once, not once
per arm. With the per-worker media/template cache that means ~50 template fetches for
the whole sweep rather than ~1500 -- which is what makes 50 workers viable against
KBase (6 workers previously tripped its rate limit while re-fetching per job).

Saves the ACTUAL per-media predictions (flux per genome x arm x media) to --calls, not
just the CP/CN/FP/FN scorecard. An earlier version kept only the aggregates and threw
the predictions away.

Models are loaded through util.gapfill_model, which constrains the leaky-reversible
reactions (LEAKY_REVERSIBLE_RXNS) before gapfilling and blocks CO2 uptake at
simulation time, so the CO2-fixation artifact cannot fire in the scored predictions.

    OLD      "<gid>.fbamodel"                 2022 build  x  PATRIC annotations
    NEW      "genomeset__<gid>.contigs..."    2025 build  x  KBase RAST re-annotation
    REBUILT  MSBuilder on the PATRIC genome   current     x  PATRIC annotations

    old  vs rebuilt -> same annotations, different build -> BUILD-CODE effect
    rebuilt vs new  -> same build, different annotations -> ANNOTATION effect

Every arm gets the IDENTICAL downstream treatment (gapfill on Carbon-Pyruvic-Acid,
simulate ecoli_biolog with CO2 uptake BLOCKED, call growth at flux > 0.01), or the
comparison means nothing.

Scored two ways:
  vs CLAYTON     replication fidelity against the 2022 run
  vs EXPERIMENT  accuracy against the real Biolog data -- the metric that matters

Design
------
- ONE MODEL PER TASK.  Models share nothing, so each (genome, arm) is a clean unit of
  work for a separate process.  Each worker fetches its own model and owns its own LP
  solver, sidestepping solver thread-safety and cobra/optlang pickling under macOS spawn.
- CHECKPOINTED.  Every completed row is appended to the TSV immediately and the run is
  resumable: re-running skips (genome, arm) pairs already present.  A 1500-job sweep on
  a laptop WILL get interrupted.
- cumulative_solution=[] is passed explicitly inside gapfill_model() -- modelseedpy's
  default is a shared mutable list that leaks reactions between models in one process.
  Workers are separate processes, but a worker handles many models, so this still bites.

Usage::

    python scripts/three_arm_sweep.py --output NBOutput/three_arm.tsv [--workers N]
                                      [--arms old,new,rebuilt] [--limit N]
"""
from __future__ import annotations

import argparse
import os
import sys
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

warnings.filterwarnings("ignore")

_PRJ = Path(__file__).resolve().parent.parent / "notebooks" / "PRJ-watershed_phenotype_replication"
sys.path.insert(0, str(_PRJ))

from util import FreeLunchModelError  # noqa: E402

FIELDS = [
    "genome_id", "arm", "n_reactions", "n_gapfilled",
    "free_lunch", "mineral_biomass", "autotrophic", "co2_biomass", "directionality_fixed",
    "clayton_accuracy", "clayton_CP", "clayton_CN", "clayton_FP", "clayton_FN",
    "exp_accuracy", "exp_CP", "exp_CN", "exp_FP", "exp_FN",
    "status", "error",
]

CALL_FIELDS = ["genome_id", "arm", "media_name", "panel", "flux", "growth",
               "experimental", "sim_2022"]

_CTX: dict = {}


def _kbase_retry(fn, tries=4, base_delay=2.0):
    """Call fn(), retrying transient KBase connection failures with backoff.

    A full sweep makes thousands of workspace calls; an occasional ConnectionError /
    ChunkedEncodingError is expected and should self-heal, not fail the job. Uses a
    plain time.sleep (Date.now/random are fine here -- this is not a workflow script).
    """
    import time
    for attempt in range(tries):
        try:
            return fn()
        except Exception as exc:  # noqa: BLE001
            transient = any(t in type(exc).__name__ for t in
                            ("Connection", "ChunkedEncoding", "Timeout", "Protocol"))
            if not transient or attempt == tries - 1:
                raise
            time.sleep(base_delay * (2 ** attempt))


def _init():
    """Load the shared, read-only inputs once per worker process.

    Every fetch here hits KBase, and this runs inside the ProcessPoolExecutor
    INITIALIZER: if it raises, the whole pool breaks and every job fails without
    writing a row. A single transient KBase blip at worker startup therefore nukes
    an entire resume (observed: an SSL handshake failure took out a 206-job run and
    produced zero rows). So the init is retried as a unit -- not just the per-job
    fetches.
    """
    from util import (get_msfba, partition_models, load_watershed_sheet,
                      get_biolog_media, load_reference_calls, get_experimental_calls)
    from util import (GAPFILL_MEDIA, GAPFILL_TEMPLATE, GAPFILL_TEMPLATE_WS, KBASEMEDIA_WS)

    def _load():
        fba = get_msfba()
        ref, _mn, _wdf = load_watershed_sheet()
        parts, _rep = partition_models(ref, fba=fba)
        return {
            "fba": fba,
            "parts": parts,
            "media": get_biolog_media(fba=fba),          # CO2 uptake blocked
            "clayton": load_reference_calls(),
            "experiment": get_experimental_calls(fba=fba),
            # fetched ONCE per worker and threaded through every job -- identical for
            # all models; re-fetching per job is what tripped KBase's rate limit
            "gapfill_media": fba.get_media(GAPFILL_MEDIA, KBASEMEDIA_WS),
            "gapfill_template": fba.get_template(GAPFILL_TEMPLATE, GAPFILL_TEMPLATE_WS),
        }

    # Broaden the transient net for init: SSL handshake failures surface as bare
    # SSLError / OSError, which _kbase_retry's name-substring test would miss.
    import time
    for attempt in range(5):
        try:
            _CTX.update(_load())
            return
        except Exception as exc:  # noqa: BLE001
            if attempt == 4:
                raise
            time.sleep(3.0 * (2 ** attempt))


def _score(rows, calls, threshold):
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


def _run_genome(job):
    """ALL arms for ONE genome. Returns (rows, calls) for every arm.

    Parallelising by GENOME rather than by (genome, arm) means the genome is fetched
    ONCE per task instead of once per arm -- a ~3x cut in KBase calls, which is what
    makes a 50-worker run viable. Combined with the per-worker media/template cache,
    50 workers issue ~50 template fetches for the whole sweep, not ~1500.
    """
    genome_id, arms = job
    rows, calls = [], []
    for arm in arms:
        r, c = _run_one((genome_id, arm))
        rows.append(r); calls.extend(c)
    return rows, calls


def _run_one(job):
    """One (genome, arm) pair. Never raises -- failures come back as (row, calls).

    ``calls`` is the ACTUAL simulation output: one entry per media with the biomass
    FLUX. Saving the flux (not just a growth bool, and not just the CP/CN/FP/FN
    scorecard) is the whole point -- it makes every downstream question answerable
    without re-simulating: the per-media tables, the growth matrices, and any
    re-scoring at a different threshold or against a different reference.

    An earlier version of this script kept only the aggregate counts and discarded
    the per-media results. That threw away the actual predictions and forced a full
    re-run. Do not do that again.
    """
    genome_id, arm = job
    row = {f: "" for f in FIELDS}
    row["genome_id"], row["arm"] = genome_id, arm
    calls = []
    try:
        from util import (gapfill_model, simulate_biolog_panel, rebuild_from_patric,
                          _latest_ref, GROWTH_THRESHOLD)

        fba, parts = _CTX["fba"], _CTX["parts"]
        template_obj = _CTX["gapfill_template"]
        if arm == "rebuilt":
            base = rebuild_from_patric(genome_id, fba=fba, template_obj=template_obj)
        else:
            if genome_id not in parts[arm]:
                row["status"], row["error"] = "skipped", f"no {arm} model"
                return row, calls
            base = _kbase_retry(lambda: fba.get_model(_latest_ref(parts[arm][genome_id])))

        cobra_model = base.model if hasattr(base, "model") else base
        row["n_reactions"] = len(cobra_model.reactions)

        mdlutl, info = gapfill_model(base, fba=fba,
                                     media_obj=_CTX["gapfill_media"],
                                     template_obj=template_obj)
        row["n_gapfilled"] = info["n_reactions_added"]
        row["autotrophic"] = int(info["autotrophic"])
        row["co2_biomass"] = round(info["co2_biomass"], 6)
        row["free_lunch"] = 0
        row["mineral_biomass"] = round(info["mineral_biomass"], 6)
        row["directionality_fixed"] = info.get("directionality_fixed", "")

        rows = simulate_biolog_panel(mdlutl, _CTX["media"], fba=fba)

        # THE ACTUAL RESULTS -- keep the flux, not just the call.
        exp_calls = _CTX["experiment"]
        clayton = _CTX["clayton"]
        col = f"{genome_id}.fbamodel"
        for r in rows:
            m = r["media_name"]
            calls.append({
                "genome_id": genome_id, "arm": arm, "media_name": m, "panel": r["panel"],
                "flux": round(r["flux"], 6), "growth": int(r["growth"]),
                "experimental": ("" if m not in exp_calls else int(bool(exp_calls[m]))),
                "sim_2022": ("" if col not in clayton.columns or m not in clayton.index
                             else int(bool(clayton[col][m]))),
            })

        column = f"{genome_id}.fbamodel"
        clayton = _CTX["clayton"]
        if column in clayton.columns:
            s = _score(rows, clayton[column], GROWTH_THRESHOLD)
            row["clayton_accuracy"] = round(s["accuracy"], 6) if s["accuracy"] is not None else ""
            for k in ("CP", "CN", "FP", "FN"):
                row[f"clayton_{k}"] = s[k]

        e = _score(rows, _CTX["experiment"], GROWTH_THRESHOLD)
        row["exp_accuracy"] = round(e["accuracy"], 6) if e["accuracy"] is not None else ""
        for k in ("CP", "CN", "FP", "FN"):
            row[f"exp_{k}"] = e[k]

        row["status"] = "ok"
    except FreeLunchModelError as exc:
        # Not a crash -- a finding. The model builds biomass from salts alone, so its
        # predictions are meaningless and it must be EXCLUDED from any mean, not scored.
        row["status"] = "free_lunch"
        row["free_lunch"] = 1
        row["error"] = str(exc)[:200]
    except Exception as exc:  # noqa: BLE001 -- a dead model must not kill the sweep
        row["status"] = "failed"
        row["error"] = f"{type(exc).__name__}: {exc}"[:200]
        traceback.print_exc(file=sys.stderr)
    return row, calls


def _done(path):
    """(genome, arm) pairs already in the output -- the sweep is resumable."""
    if not path.exists():
        return set()
    seen = set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, ai = header.index("genome_id"), header.index("arm")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(gi, ai):
                seen.add((parts[gi], parts[ai]))
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--calls", default=None,
                    help="ALSO write the per-media simulation results (flux per genome x "
                         "arm x media). This is the actual prediction data -- always set it.")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2),
                    help="Parallel genome-tasks. 50 on H100; keep low (2-3) on a laptop.")
    ap.add_argument("--arms", default="old,new,rebuilt")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    _init()  # in the parent, to enumerate genomes
    genomes = sorted(_CTX["parts"]["old"])
    if args.limit:
        genomes = genomes[:args.limit]
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]

    already = _done(out)
    # one task per GENOME, carrying only the arms still outstanding for it
    jobs = []
    for g in genomes:
        todo = [a for a in arms if (g, a) not in already]
        if todo:
            jobs.append((g, todo))
    n_pairs = sum(len(a) for _g, a in jobs)
    print(f"{len(genomes)} genomes x {len(arms)} arms = {len(genomes)*len(arms)} pairs; "
          f"{len(already)} already done; {n_pairs} pairs across {len(jobs)} genome-tasks "
          f"on {args.workers} workers", flush=True)
    if not jobs:
        print("nothing to do", flush=True)
        return

    calls_path = Path(args.calls) if args.calls else out.with_name(out.stem + "_calls.tsv")
    calls_path.parent.mkdir(parents=True, exist_ok=True)

    new_file = not out.exists()
    new_calls = not calls_path.exists()
    with open(out, "a", buffering=1) as fh, open(calls_path, "a", buffering=1) as cfh:
        if new_file:
            fh.write("\t".join(FIELDS) + "\n")
        if new_calls:
            cfh.write("\t".join(CALL_FIELDS) + "\n")
        done = 0
        with ProcessPoolExecutor(max_workers=args.workers, initializer=_init) as pool:
            futures = {pool.submit(_run_genome, j): j for j in jobs}
            for fut in as_completed(futures):
                rows, calls = fut.result()
                for row in rows:
                    fh.write("\t".join(str(row.get(f, "")) for f in FIELDS) + "\n")
                # stream the actual per-media results straight to disk
                for c in calls:
                    cfh.write("\t".join(str(c.get(f, "")) for f in CALL_FIELDS) + "\n")
                done += 1
                bad = [r for r in rows if r["status"] not in ("ok", "skipped")]
                if done % 10 == 0 or bad:
                    gid = rows[0]["genome_id"] if rows else "?"
                    st = ",".join(f"{r['arm']}={r['status']}" for r in rows)
                    print(f"[{done}/{len(jobs)} genomes] {gid}: {st}", flush=True)
    print(f"done -- {done} rows -> {out}", flush=True)
    print(f"         per-media results -> {calls_path}", flush=True)


if __name__ == "__main__":
    main()
