#!/usr/bin/env python
"""Reduce the four-dataset sweep TSVs to a compact JSON the dashboard embeds.

Reads scripts/four_dataset_sweep.py's output and writes docs/data/dashboard.json:
per-dataset fidelity/accuracy stats, the annotation/build decomposition, and the
per-genome and per-substrate tables in the `count (% changed vs Original 2022p)`
format. Re-run after the sweep finishes to refresh the dashboard.

    python scripts/build_dashboard_data.py [sweep.tsv] [out.json] [calls.tsv]

The four datasets (Chris's labels, used verbatim):
    original_2022p  Clayton's own 2022 spreadsheet calls -- the baseline, not simulated
    resim_2022p     our resimulation of the 2022 models
    resim_2025k     the 2025 KBase re-annotated models
    resim_2026p     the 2026 PATRIC rebuild
Every count in columns 2-4 is reported with the % of growth conditions that changed
relative to Original 2022p.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
NB = REPO / "notebooks/PRJ-watershed_phenotype_replication"
DEFAULT_TSV = NB / "NBOutput/four_dataset_sweep.tsv"
DEFAULT_CALLS = NB / "NBOutput/four_dataset_sweep_calls.tsv"
DEFAULT_OUT = REPO / "docs/data/dashboard.json"

# Order matters -- this is the column order everywhere in the dashboard.
DATASETS = ["original_2022p", "resim_2022p", "resim_2025k", "resim_2026p"]
SIM = ["resim_2022p", "resim_2025k", "resim_2026p"]   # the three simulated arms
BASELINE = "original_2022p"


def _stats(series):
    s = pd.to_numeric(series, errors="coerce").dropna()
    n = int(s.count())
    if n == 0:
        return {"n": 0, "mean": None, "se": None}
    sd = float(s.std()) if n > 1 else 0.0
    return {"n": n, "mean": float(s.mean()), "se": sd / n ** 0.5 if n > 1 else 0.0}


def _effect(paired, a, b):
    """Paired mean delta (b - a) on a per-genome accuracy column, in fractional units."""
    d = (paired[b] - paired[a]).dropna()
    n = int(d.count())
    if n < 2:
        return None
    se = float(d.std()) / n ** 0.5
    mean = float(d.mean())
    return {"delta": mean, "se": se, "n": n, "real": bool(abs(mean) > 2 * se)}


def _media_table(calls_path):
    """Per-substrate summary from the actual per-media calls.

    For each of the 324 Biolog conditions and each dataset: how many models predict
    growth, and what % of those calls changed vs Original 2022p (per-genome flip rate).
    Sorted by the rebuild's change rate -- the substrates the reconstruction moved most.
    """
    if not Path(calls_path).exists():
        return []
    c = pd.read_csv(calls_path, sep="\t", dtype={"genome_id": str})
    c["growth"] = pd.to_numeric(c["growth"], errors="coerce")
    c["original_2022p"] = pd.to_numeric(c["original_2022p"], errors="coerce")
    c["experimental"] = pd.to_numeric(c["experimental"], errors="coerce")

    rows = []
    for (media, panel), grp in c.groupby(["media_name", "panel"], sort=False):
        entry = {"media_name": media, "panel": panel}
        for ds in DATASETS:
            d = grp[grp["dataset"] == ds]
            g = d["growth"].dropna()
            grow = int(g.sum()) if len(g) else None
            pct = None
            if ds != BASELINE:
                cmp = d.dropna(subset=["growth", "original_2022p"])
                if len(cmp):
                    changed = int((cmp["growth"].astype(int) != cmp["original_2022p"].astype(int)).sum())
                    pct = round(100.0 * changed / len(cmp), 1)
            entry[ds] = {"grow": grow, "pct": pct}
        exp = grp["experimental"].dropna()
        entry["experimental"] = int(exp.iloc[0]) if len(exp) else None
        rows.append(entry)
    rows.sort(key=lambda r: -(r["resim_2026p"]["pct"] or 0))
    return rows


def _distribution(media_table):
    """Characterise HOW the per-substrate changes are distributed.

    Two things the aggregate %changed hides and a reviewer will ask about:
      * concentration -- is the ~10% mean spread evenly or carried by a few
        substrates? (median vs the count changing >50%.)
      * direction of the 2026p rebuild -- net gain or loss of growth calls, and
        on which panels, vs Original 2022p.
    """
    if not media_table:
        return {}
    import statistics as st

    concentration = {}
    for ds in SIM:
        p = [r[ds]["pct"] for r in media_table if r[ds]["pct"] is not None]
        if p:
            concentration[ds] = {
                "median": round(st.median(p), 1),
                "over50": int(sum(1 for x in p if x > 50)),
                "n": len(p),
            }

    gained, lost, gains = 0, 0, []
    for r in media_table:
        o, n = r["original_2022p"]["grow"], r["resim_2026p"]["grow"]
        if o is None or n is None:
            continue
        if n - o > 50:
            gained += 1
            gains.append({"media": r["media_name"].split("-", 1)[-1], "from": o, "to": n})
        elif o - n > 50:
            lost += 1
    gains.sort(key=lambda g: g["to"] - g["from"], reverse=True)
    return {
        "concentration": concentration,
        "rebuild_shift": {"gained": gained, "lost": lost, "sample_gains": gains[:3]},
    }


def main(tsv=DEFAULT_TSV, out=DEFAULT_OUT, calls=DEFAULT_CALLS):
    df = pd.read_csv(tsv, sep="\t", dtype={"genome_id": str})
    for col in ("clayton_accuracy", "exp_accuracy", "n_growth", "pct_changed_vs_orig",
                "exp_CP", "exp_CN", "exp_FP", "exp_FN", "free_lunch"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    ok = df[df["status"] == "ok"].copy()

    # ---- per-dataset fidelity + accuracy -----------------------------------
    datasets = {}
    for ds in DATASETS:
        a = ok[ok["dataset"] == ds]
        if a.empty:
            continue
        datasets[ds] = {
            "n_ok": int(len(a)),
            "n_growth": _stats(a["n_growth"]),
            "pct_changed": (None if ds == BASELINE else _stats(a["pct_changed_vs_orig"])),
            "clayton": (None if ds == BASELINE else _stats(a["clayton_accuracy"])),
            "exp": _stats(a["exp_accuracy"]),
            "free_lunch": int(a["free_lunch"].fillna(0).sum()),
        }

    # ---- decomposition on accuracy vs experiment, paired per genome --------
    wide = ok.pivot_table(index="genome_id", columns="dataset", values="exp_accuracy")
    present = [d for d in SIM if d in wide.columns]
    paired = wide.dropna(subset=present)
    decomposition = {}
    if {"resim_2022p", "resim_2025k"} <= set(present):
        decomposition["annotation"] = _effect(paired, "resim_2022p", "resim_2025k")
    if {"resim_2022p", "resim_2026p"} <= set(present):
        decomposition["build"] = _effect(paired, "resim_2022p", "resim_2026p")

    # ---- per-genome table: count (% changed) + accuracies ------------------
    # K-12 growers per genome (exp CP+FN, constant across the simulated datasets).
    # original_2022p carries no exp confusion counts, so read it from a sim row.
    exp_grow = {}
    for gid, g in ok.groupby("genome_id"):
        sim = g[g["dataset"].isin(SIM)].dropna(subset=["exp_CP", "exp_FN"])
        if len(sim):
            r = sim.iloc[0]
            exp_grow[gid] = int(r["exp_CP"]) + int(r["exp_FN"])

    genome_table = []
    for gid, g in ok.groupby("genome_id"):
        entry = {"genome_id": gid, "exp_grow": exp_grow.get(gid)}
        for ds in DATASETS:
            r = g[g["dataset"] == ds]
            if r.empty:
                entry[ds] = None
                continue
            r = r.iloc[0]
            entry[ds] = {
                "grow": int(r["n_growth"]) if pd.notna(r["n_growth"]) else None,
                "pct": (None if ds == BASELINE or pd.isna(r["pct_changed_vs_orig"])
                        else round(float(r["pct_changed_vs_orig"]), 1)),
                "clayton": (None if ds == BASELINE or pd.isna(r["clayton_accuracy"])
                            else round(float(r["clayton_accuracy"]), 4)),
                "exp": (None if pd.isna(r["exp_accuracy"]) else round(float(r["exp_accuracy"]), 4)),
            }
        genome_table.append(entry)
    genome_table.sort(key=lambda e: (e.get("resim_2026p") or {}).get("pct") or 0, reverse=True)

    media_table = _media_table(calls)

    # ---- meta --------------------------------------------------------------
    n_genomes = int(df["genome_id"].nunique())
    n_failed = int((df["status"] == "failed").sum())
    expected = n_genomes * len(DATASETS)
    payload = {
        "meta": {
            "n_genomes": n_genomes,
            "datasets": DATASETS,
            "rows": int(len(df)),
            "expected": expected,
            "failed": n_failed,
            "complete": bool(n_failed == 0 and len(ok) >= expected),
        },
        "datasets": datasets,
        "decomposition": decomposition,
        "genome_table": genome_table,
        "media_table": media_table,
        "distribution": _distribution(media_table),
    }

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    m = payload["meta"]
    print(f"wrote {out}  ({m['rows']}/{m['expected']} rows, {n_genomes} genomes, "
          f"{'COMPLETE' if m['complete'] else 'PARTIAL'}, {n_failed} failed)")
    for ds in DATASETS:
        d = datasets.get(ds)
        if not d:
            continue
        cl = d["clayton"]["mean"] if d["clayton"] else None
        pc = d["pct_changed"]["mean"] if d["pct_changed"] else None
        ng, ex = d["n_growth"]["mean"], d["exp"]["mean"]
        print(f"  {ds:16s} n_growth={'--' if ng is None else f'{ng:.1f}'}  "
              f"clayton_acc={'--' if cl is None else f'{cl:.3f}'}  "
              f"%changed={'--' if pc is None else f'{pc:.1f}'}  "
              f"exp_acc={'--' if ex is None else f'{ex:.3f}'}  free_lunch={d['free_lunch']}")


if __name__ == "__main__":
    tsv = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_TSV
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUT
    calls = Path(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_CALLS
    main(tsv, out, calls)
