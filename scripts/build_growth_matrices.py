#!/usr/bin/env python
"""Growth matrices for the collaborators — one media x genome sheet per dataset.

Pivots four_dataset_sweep_calls.tsv into 0/1 growth matrices (324 Biolog media rows x
519 genome columns), one sheet per dataset, laid out to match Clayton's original
BiologMacTestOutput_On-Off.csv so it can be overlaid directly. These go to Aaron as an
email attachment; they are NOT part of the public repo.

    python scripts/build_growth_matrices.py [out.xlsx]

Sheets:
    Original_2022p     Clayton's own 2022 calls (identical to his file — the baseline)
    Resimulated_2022p  our resimulation of the 2022 models
    Resimulated_2025k  the 2025 KBase re-annotated models
    Resimulated_2026p  the 2026 PATRIC rebuild
    Experimental       measured K-12 Biolog growth (one column)
    README             legend
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
NB = REPO / "notebooks/PRJ-watershed_phenotype_replication"
CALLS = NB / "NBOutput/four_dataset_sweep_calls.tsv"
CLAYTON = REPO / "notebooks/data/BiologMacTestOutput_On-Off.csv"
DEFAULT_OUT = NB / "NBOutput/growth_matrices_for_collaborators.xlsx"

DATASETS = [
    ("original_2022p", "Original_2022p"),
    ("resim_2022p", "Resimulated_2022p"),
    ("resim_2025k", "Resimulated_2025k"),
    ("resim_2026p", "Resimulated_2026p"),
]


def main(out=DEFAULT_OUT):
    c = pd.read_csv(CALLS, sep="\t", dtype={"genome_id": str})
    c["growth"] = pd.to_numeric(c["growth"], errors="coerce")

    # Row/column order: follow Clayton's original file so the sheets overlay it exactly.
    row_order = col_order = None
    if CLAYTON.exists():
        cl = pd.read_csv(CLAYTON, index_col=0)
        row_order = list(cl.index)
        col_order = [str(x)[:-len(".fbamodel")] if str(x).endswith(".fbamodel") else str(x)
                     for x in cl.columns]

    def matrix(ds):
        m = (c[c["dataset"] == ds]
             .pivot_table(index="media_name", columns="genome_id", values="growth", aggfunc="max"))
        if row_order is not None:
            m = m.reindex(index=[r for r in row_order if r in m.index],
                          columns=[g for g in col_order if g in m.columns])
        return m.astype("Int64")

    with pd.ExcelWriter(out, engine="openpyxl") as xl:
        readme = pd.DataFrame({
            "Watershed E. coli phenotype replication — growth matrices": [
                "Each dataset sheet: rows = 324 Biolog conditions, columns = 519 genomes,",
                "values = 1 (predicted growth) / 0 (no growth). Layout matches the original",
                "2022 BiologMacTestOutput_On-Off.csv so sheets can be overlaid directly.",
                "",
                "Original_2022p     Clayton's own 2022 output (the baseline).",
                "Resimulated_2022p  the same 2022 models, resimulated with current code.",
                "Resimulated_2025k  the KBase RAST re-annotated models (2025).",
                "Resimulated_2026p  a full rebuild from the original genomes (2026, ATP-safe).",
                "Experimental       measured K-12 Biolog growth calls (reference).",
                "",
                "All four datasets: gap-fill on Carbon-Pyruvic-Acid, simulate the Biolog panel",
                "with CO2 uptake blocked, call growth at biomass flux > 0.01. 519 genomes,",
                "0 simulation failures.",
            ]
        })
        readme.to_excel(xl, sheet_name="README", index=False)
        shapes = {}
        for ds, sheet in DATASETS:
            m = matrix(ds)
            m.to_excel(xl, sheet_name=sheet)
            shapes[sheet] = m.shape
        # Experimental: constant per media; take the first non-null per media_name.
        exp = (c.dropna(subset=["experimental"])
               .drop_duplicates("media_name").set_index("media_name")["experimental"])
        exp = pd.to_numeric(exp, errors="coerce").astype("Int64").rename("experimental_growth")
        if row_order is not None:
            exp = exp.reindex([r for r in row_order if r in exp.index])
        exp.to_frame().to_excel(xl, sheet_name="Experimental")

    print(f"wrote {out}")
    for sheet, (nr, nc) in shapes.items():
        print(f"  {sheet:20s} {nr} media x {nc} genomes")
    print(f"  Experimental         {int(exp.notna().sum())} media with a measured call")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT)
