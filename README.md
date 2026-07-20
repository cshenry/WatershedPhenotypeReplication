# Watershed *E. coli* — Phenotype Model Replication

Reproduces the Hope College team's Summer-2022 *E. coli* Biolog phenotype predictions
with current KBase/ModelSEED code, then runs the identical pipeline over newer model
sets so that any divergence can be attributed to a specific cause.

**📊 [Live results dashboard →](https://cshenry.github.io/WatershedPhenotypeReplication/)**

**🔧 [Setup & run guide →](docs/SETUP.md)** — start here if you want to run this yourself.

---

## The four datasets

Everything is measured over the same **519 genomes × 324 Biolog conditions**. One
baseline and three resimulations, each through the *identical* downstream pipeline —
gap-fill on `Carbon-Pyruvic-Acid`, simulate the Biolog panel with CO₂ uptake blocked,
call growth at biomass flux > 0.01.

| Dataset | What it is | Annotation | Reconstruction |
|---|---|---|---|
| **Original 2022p** | the collaborators' own 2022 output — the baseline | PATRIC (2022) | 2022, as published |
| **Resimulated 2022p** | the same 2022 models, resimulated with current code | PATRIC (2022) | 2022 models |
| **Resimulated 2025k** | the KBase RAST re-annotated models | KBase RAST (2025) | 2025 build |
| **Resimulated 2026p** | a full rebuild from the original genomes | PATRIC (2022) | 2026 (MSBuilder, ATP-safe) |

`2022p → 2025k` isolates the **re-annotation effect**; `2022p → 2026p` isolates the
**rebuild effect**.

## What this study found

- **Current code reproduces ~90% of the 2022 calls.** Resimulating the original 2022
  models matches the collaborators' own growth/no-growth calls on **89.8%** of
  conditions.
- **Re-annotation and rebuilding each shift only ~10% of growth conditions**
  (2025k: 11.0%, 2026p: 10.0%) relative to the 2022 baseline.
- **Neither change costs accuracy.** Scored against the *experimental* K-12 Biolog data,
  all three resimulated sets land at ~0.64–0.65. Both re-annotation (+0.48 pp) and the
  rebuild (+0.50 pp) nudge accuracy up slightly — small, but statistically real across
  519 genomes.
- **The changes are concentrated, not uniform.** The median substrate changes only
  0.5–1.7% of its calls, but a tail of 32–35 substrates (of 324) changes more than half.
- **A methodological finding: the models can fix CO₂ — and shouldn't.** Freshly
  gap-filled *E. coli* models can build biomass from CO₂ + minerals with no organic
  carbon. The cause is **reaction directionality** (`rxn05759`, `rxn03978`), not
  gap-filling — and it slips past both ATP correction (the leak is electrons, not ATP)
  and mass-balance checks (the formulas balance; only the charge is wrong). Both
  reactions are constrained and CO₂ uptake blocked at simulation time, so it cannot
  affect the predictions reported here.

See the [dashboard](https://cshenry.github.io/WatershedPhenotypeReplication/) for the
per-substrate and per-genome breakdowns.

## Repository layout

```
scripts/                     the compute pipeline (run in this order)
  build_2026p.py             build the 2026p arm from the PATRIC genomes
  gapfill_all.py             gap-fill all three model sets to <name>.MMGF
  four_dataset_sweep.py      simulate the Biolog panel across the four datasets
  build_dashboard_data.py    reduce the sweep to docs/data/dashboard.json
  build_growth_matrices.py   export media × genome growth matrices (xlsx)
docs/                        the dashboard (GitHub Pages serves this directory)
notebooks/PRJ-watershed_phenotype_replication/
  util.py                    %run-loadable helpers: media, gap-fill, simulate, score
  0*.ipynb                   the exploratory analysis notebooks (outputs cleared)
tests/                       offline unit tests
```

## Running it

The pipeline needs a KBase account **and access to the source narratives** — see
**[docs/SETUP.md](docs/SETUP.md)** for the full install and run guide, including the
narratives you'll need shared with you and the 2022 data files you supply yourself.

## Provenance

The 2022 method is **not inferred**. It is Clayton Piehl's own notebook and spreadsheet
output, supplied by the Hope College team and used here as ground truth. Those files are
the collaborators' unpublished data and are therefore **not committed to this public
repository** — see [docs/SETUP.md](docs/SETUP.md) for where to place your copies.

---

Argonne National Laboratory, in collaboration with the Hope College Global Water
Research Institute.
