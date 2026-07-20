#!/usr/bin/env python
"""Analyse the three-arm sweep: model quality, accuracy, and the annotation/build split.

Reads the TSV written by three_arm_sweep.py (partial runs are fine -- it reports what
is there and says so) and answers three questions:

  1. MODEL QUALITY   how many models in each arm are thermodynamically broken?
                     free_lunch (biomass from minerals alone -- incurable, excluded)
                     autotrophic (biomass from CO2 + minerals -- curable, kept)

  2. ACCURACY        vs the EXPERIMENTAL Biolog data, which is what "right" means.
                     Agreement with the 2022 run only says we reproduced a pipeline
                     that we now know leaked CO2.

  3. DECOMPOSITION   old vs rebuilt -> same annotations, different build -> BUILD effect
                     rebuilt vs new -> same build, different annotations -> ANNOTATION effect
                     Paired per genome, so a genome that fails in one arm cannot skew it.

Free-lunch models are EXCLUDED from every accuracy number. Their predictions are
meaningless (they grow on all 324 media regardless of content), and averaging them in
would silently drag an arm's accuracy toward the base rate.

Usage::

    python scripts/analyze_three_arm.py NBOutput/three_arm_sweep.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path
from statistics import mean, stdev

import pandas as pd

ARMS = ("old", "new", "rebuilt")


def main(path):
    df = pd.read_csv(path, sep="\t", dtype={"genome_id": str})
    total_expected = 519 * 3
    print(f"{len(df)} rows of an expected {total_expected} "
          f"({len(df)/total_expected:.0%} of the sweep)\n")
    if len(df) < total_expected:
        print("*** PARTIAL RUN -- numbers below are from the rows completed so far ***\n")

    # ---- 1. model quality -------------------------------------------------
    print("=" * 74)
    print("1. MODEL QUALITY -- how many models are thermodynamically broken?\n")
    print(f"{'arm':10s} {'n':>5s} {'FREE LUNCH':>12s} {'autotrophic':>13s} {'clean':>8s}")
    print("-" * 74)
    for arm in ARMS:
        a = df[df["arm"] == arm]
        if a.empty:
            continue
        n = len(a)
        fl = int((a["status"] == "free_lunch").sum())
        auto = int(a[a["status"] == "ok"]["autotrophic"].fillna(0).astype(float).sum())
        clean = n - fl - auto
        print(f"{arm:10s} {n:5d} {fl:8d} ({fl/n:3.0%}) {auto:9d} ({auto/n:3.0%}) {clean:5d} ({clean/n:3.0%})")
    print("\n  FREE LUNCH  = builds biomass from MINERALS ALONE, no carbon of any kind.")
    print("                Thermodynamically impossible. No media-level fix exists.")
    print("                Predictions are meaningless -> EXCLUDED from all accuracy below.")
    print("  autotrophic = builds biomass from CO2 + minerals. Curable by blocking CO2")
    print("                uptake (which the media already specify) -> kept, predictions sound.")

    ok = df[df["status"] == "ok"].copy()
    for col in ("exp_accuracy", "clayton_accuracy"):
        ok[col] = pd.to_numeric(ok[col], errors="coerce")

    # ---- 2. accuracy ------------------------------------------------------
    print("\n" + "=" * 74)
    print("2. ACCURACY vs the EXPERIMENTAL Biolog data (free-lunch models excluded)\n")
    print(f"{'arm':10s} {'n':>5s} {'accuracy':>10s} {'sd':>7s} | {'vs 2022 run':>12s} {'sd':>7s}")
    print("-" * 74)
    for arm in ARMS:
        a = ok[(ok["arm"] == arm) & ok["exp_accuracy"].notna()]
        if a.empty:
            continue
        e, c = a["exp_accuracy"], a["clayton_accuracy"].dropna()
        print(f"{arm:10s} {len(a):5d} {e.mean():10.4f} {(e.std() if len(e)>1 else 0):7.4f} | "
              f"{(c.mean() if len(c) else float('nan')):12.4f} {(c.std() if len(c)>1 else 0):7.4f}")

    # ---- 3. decomposition -------------------------------------------------
    print("\n" + "=" * 74)
    print("3. DECOMPOSITION -- what actually moved the models?\n")
    wide = ok.pivot_table(index="genome_id", columns="arm", values="exp_accuracy")
    have = [a for a in ARMS if a in wide.columns]
    paired = wide.dropna(subset=have)
    print(f"paired on {len(paired)} genomes with a usable model in all {len(have)} arms\n")
    if len(paired) < 2:
        print("  not enough paired genomes yet")
        return

    def delta(a, b, label, meaning):
        d = (paired[b] - paired[a]).dropna()
        if d.empty:
            return
        m, s = d.mean(), (d.std() if len(d) > 1 else 0.0)
        se = s / (len(d) ** 0.5) if len(d) > 1 else 0.0
        sig = "" if se == 0 else ("  (>2 SE -- real)" if abs(m) > 2 * se else "  (within noise)")
        print(f"  {label:22s} {m:+.4f}  +/- {se:.4f} SE   n={len(d)}{sig}")
        print(f"  {'':22s} {meaning}\n")

    if "old" in have and "rebuilt" in have:
        delta("old", "rebuilt", "BUILD-CODE effect",
              "old -> rebuilt: same PATRIC annotations, 2022 build -> current build")
    if "rebuilt" in have and "new" in have:
        delta("rebuilt", "new", "ANNOTATION effect",
              "rebuilt -> new: same modern build, PATRIC -> KBase re-annotation")
    if "old" in have and "new" in have:
        delta("old", "new", "TOTAL old -> new",
              "what the collaborators actually observed (both changes at once)")

    print("-" * 74)
    print("Reading it: if the TOTAL old->new change is ~0 while the individual effects")
    print("are not, the two changes cancelled -- accuracy held while the ERRORS MOVED.")
    print("That is exactly what Clayton reported and could not explain.")


if __name__ == "__main__":
    p = Path(sys.argv[1] if len(sys.argv) > 1 else
             "notebooks/PRJ-watershed_phenotype_replication/NBOutput/three_arm_sweep.tsv")
    if not p.exists():
        sys.exit(f"no sweep output at {p}")
    main(p)
