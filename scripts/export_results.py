#!/usr/bin/env python
"""Export the per-genome three-arm results as a shareable spreadsheet.

Reduces the sweep TSV to one row per genome with, for each arm (old / new / rebuilt):
predicted growers, accuracy vs the experimental K-12 Biolog data, and the model-quality
flag. Aggregate summary only -- NOT the per-(genome, media) growth matrices (the sweep
did not save per-media calls; those need a re-run of the simulations).

    python scripts/export_results.py [sweep.tsv] [out.xlsx|out.csv]
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ARMS = ("old", "new", "rebuilt")
REPO = Path(__file__).resolve().parent.parent
DEFAULT_TSV = REPO / "notebooks/PRJ-watershed_phenotype_replication/NBOutput/three_arm_sweep.tsv"


def main(tsv, out):
    df = pd.read_csv(tsv, sep="\t", dtype={"genome_id": str})
    for c in ("exp_CP", "exp_CN", "exp_FP", "exp_FN", "exp_accuracy", "clayton_accuracy"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    rows = []
    for gid, g in df.groupby("genome_id"):
        row = {"genome_id": gid}
        for arm in ARMS:
            r = g[g["arm"] == arm]
            if r.empty:
                continue
            r = r.iloc[0]
            status = r["status"]
            if status == "ok":
                cp, cn, fp, fn = (int(r[f"exp_{k}"]) for k in ("CP", "CN", "FP", "FN"))
                row[f"{arm}_predicted_growers"] = cp + fp
                row[f"{arm}_accuracy_vs_K12"] = round(float(r["exp_accuracy"]), 4)
                row[f"{arm}_accuracy_vs_2022"] = (round(float(r["clayton_accuracy"]), 4)
                                                  if pd.notna(r["clayton_accuracy"]) else "")
                row[f"{arm}_quality"] = ("autotrophic" if int(float(r["autotrophic"] or 0))
                                         else "clean")
                if "exp_grow" not in row:
                    row["K12_experimental_growers"] = cp + fn
            else:
                row[f"{arm}_predicted_growers"] = ""
                row[f"{arm}_accuracy_vs_K12"] = ""
                row[f"{arm}_accuracy_vs_2022"] = ""
                row[f"{arm}_quality"] = status  # free_lunch / failed
        rows.append(row)

    cols = ["genome_id", "K12_experimental_growers"]
    for arm in ARMS:
        cols += [f"{arm}_predicted_growers", f"{arm}_accuracy_vs_K12",
                 f"{arm}_accuracy_vs_2022", f"{arm}_quality"]
    out_df = pd.DataFrame(rows).reindex(columns=cols).sort_values("genome_id")

    out = Path(out)
    if out.suffix.lower() in (".xlsx", ".xls"):
        try:
            out_df.to_excel(out, index=False, sheet_name="per_genome_results")
        except ModuleNotFoundError:
            out = out.with_suffix(".csv")
            out_df.to_csv(out, index=False)
            print("(openpyxl not installed -> wrote CSV instead)")
    else:
        out_df.to_csv(out, index=False)
    print(f"wrote {out}  ({len(out_df)} genomes x {len(cols)} columns)")


if __name__ == "__main__":
    tsv = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_TSV
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else (REPO / "per_genome_results.csv")
    main(tsv, out)
