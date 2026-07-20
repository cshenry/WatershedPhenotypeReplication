# Setup & run guide

How to install this pipeline and reproduce the results yourself. Written to be
followed start-to-finish, including by a coding agent working on your behalf.

There are **three things you need that are not in this repository**: access to the
KBase narratives, your own 2022 data files, and a KBase auth token. Steps 1–3 below.

---

## Step 1 — Ask Chris to share the KBase narratives (do this first)

Nothing in this pipeline runs without read access to the KBase workspaces holding the
genomes and models. **They are not public — email Chris Henry (chenry@anl.gov) and ask
him to share these narratives with your KBase username:**

| Workspace | What it holds | Needed for |
|---|---|---|
| **265353** | the four-dataset model corpus — 519 × 3 base models plus 519 × 3 gap-filled `.MMGF` models | `four_dataset_sweep.py`, `gapfill_all.py` |
| **119455** | the source *E. coli* genomes and the original 2022 models | `build_2026p.py`, genome lookups |
| **220871** | the models workspace referenced by the artifact inventory | `01_artifact_inventory.ipynb` |

Ask for **read access** (viewer is enough — you don't need write unless you intend to
save new models). Give him your exact KBase username.

Two more workspaces are already public and need no action: `KBaseMedia` (the Biolog
media) and `93541` (the `ecoli_biolog` phenotype set with the measured K-12 calls).

> **If you only want the numbers, not a re-run:** the
> [dashboard](https://cshenry.github.io/WatershedPhenotypeReplication/) and the growth
> matrices Chris sent by email already contain every result. You only need narrative
> access to *recompute* them.

## Step 2 — Put your 2022 data files in place

The 2022 output is **your team's unpublished data**, so it is deliberately not committed
here. Copy your own files into `notebooks/data/` (create the directory):

```
notebooks/data/
  BiologMacTestOutput_On-Off.csv    # Clayton's 324 × 519 growth-call matrix (the baseline)
  Watershed_Mastersheet_Sum22.csv   # the 519-genome roster
  BioLogMacTest.ipynb               # Clayton's original notebook (reference only)
```

`BiologMacTestOutput_On-Off.csv` must be the on/off matrix: media names down the first
column, one column per model named `<genome_id>.fbamodel`, values `0`/`1`. The pipeline
reads it as the **Original 2022p** baseline that every "% changed" is measured against.

## Step 3 — Install

Python **3.11** is required (3.12+ has not been tested against this stack).

```bash
git clone https://github.com/cshenry/WatershedPhenotypeReplication.git
cd WatershedPhenotypeReplication
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

> ### ⚠️ Known dependency caveat — read this before reporting a bug
>
> `requirements.txt` installs KBUtilLib from its `main` branch. Several fixes this
> pipeline depends on may not have landed on `main` yet — most importantly a
> `NameError: name 'genome' is not defined` raised by
> `kb_gapfill_metabolic_models` on **every** gap-fill call.
>
> If you hit that error, KBUtilLib is too old. Ask Chris which branch or commit to pin,
> and install it directly:
> ```bash
> pip install "git+https://github.com/cshenry/KBUtilLib@<branch-or-commit>#egg=kbutillib"
> ```
> The same applies to `modelseedpy` — the pinned PyPI release may lag the version this
> was developed against. If `MSBuilder` is missing on import, ask Chris for the
> ModelSEEDpy branch.

### KBase token

Every script authenticates to KBase. Get your token from
[narrative.kbase.us](https://narrative.kbase.us) (Account → Developer Tokens) and make
it available either way:

```bash
export KB_AUTH_TOKEN="<your-token>"
# or write it to the file the loader checks:
mkdir -p ~/.kbase && echo "<your-token>" > ~/.kbase/token
```

Verify the install and your access in one shot:

```bash
python -c "
import sys; sys.path.insert(0,'notebooks/PRJ-watershed_phenotype_replication')
from util import get_msfba
fba = get_msfba()
n = len(fba.list_ws_objects(265353, type='KBaseFBA.FBAModel'))
print(f'OK — {n} models visible in narrative 265353')   # expect 3114
"
```

If that prints ~3114 models, you are ready. If it raises a permissions error, Step 1
has not completed.

---

## Running the pipeline

The stages are **checkpointed and resumable** — every one appends to a TSV and skips
work already marked `ok`, so an interrupted run is safe to restart with the same
command. All are parallel; `--workers 50` suits a large machine, `--workers 2` a laptop.

Run them in this order (or skip to stage 3 — the models already exist in 265353):

```bash
# 1. Build the 2026p arm (519 models, ~30 s each)
python scripts/build_2026p.py --workers 50

# 2. Gap-fill all three sets to .MMGF (519 × 3, ~500 s each — this is the long one)
python scripts/gapfill_all.py --workers 50

# 3. Simulate the Biolog panel across the four datasets (~36 s per model)
python scripts/four_dataset_sweep.py --workers 50

# 4. Regenerate the dashboard data + the growth matrices
python scripts/build_dashboard_data.py
python scripts/build_growth_matrices.py
```

Outputs land in `notebooks/PRJ-watershed_phenotype_replication/NBOutput/`
(gitignored): the sweep summary, the per-media calls, and the growth-matrix workbook.

**Smoke-test first.** Every stage takes `--limit N` to run a couple of items before
committing to the full set:

```bash
python scripts/gapfill_all.py --limit 2 --workers 1
python scripts/four_dataset_sweep.py --limit 2 --workers 1 --media-limit 5
```

(`--media-limit` is a smoke-test knob only — it simulates a subset of the 324 media, so
the scored numbers are meaningless. Never use it for a real run.)

### Viewing the dashboard locally

`docs/` is a static site. After `build_dashboard_data.py` regenerates
`docs/data/dashboard.json`:

```bash
cd docs && python3 -m http.server 8801
# open http://localhost:8801/
```

It must be served over HTTP — opening `index.html` from the filesystem fails, because
the page fetches `data/dashboard.json`.

### The notebooks

`notebooks/PRJ-watershed_phenotype_replication/*.ipynb` are the exploratory analyses
(outputs cleared). They follow a `%run util.py` convention: **every cell starts with
`%run util.py`** so it can be re-run independently after a kernel restart. Start
JupyterLab from the repository root:

```bash
pip install jupyterlab && jupyter lab
```

Note that the notebooks reflect an earlier three-arm framing of the analysis; the
four-dataset results in the dashboard come from `scripts/four_dataset_sweep.py`.

---

## How the pipeline works

Each of the three simulated datasets goes through the *identical* treatment, which is
what makes them comparable:

1. **Gap-fill** on `Carbon-Pyruvic-Acid`, ATP-safe, with two leaky-reversible reactions
   (`rxn05759`, `rxn03978`) constrained to their physiological direction *first*. The
   saved 2022 models are auxotrophic — without gap-filling they predict zero growth on
   everything, glucose included.
2. **Simulate** the `ecoli_biolog` phenotype set (324 conditions) with missing
   transporters auto-added and **CO₂ uptake blocked** (restoring a bound the media
   already declare — see the CO₂ finding in the README).
3. **Score** the growth/no-growth calls twice: against the 2022 output (replication
   fidelity) and against the experimental K-12 Biolog data (accuracy).
4. **Verify** each gap-filled model actually grows before it is scored — a model that
   fails to grow on its own gap-fill medium is flagged `failed`, not silently reported.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `NameError: name 'genome' is not defined` | KBUtilLib too old — see the dependency caveat above. |
| `Either set callback URL...` | The native ontology path wasn't enabled. The scripts set it themselves; if you call the library directly, set `native_ontology` on the **delegate** (`k.recon._delegate.native_ontology = True`) — the Impl wrapper proxies `__getattr__` but not `__setattr__`, so assigning on the wrapper silently does nothing. |
| `reference cannot be null or the empty string` | A model was saved without a `template_ref`. Rebuild it with `build_2026p.py`. |
| `No object with name Carbon-Pyruvic-Acid` | Media reference not workspace-qualified — it lives in `KBaseMedia`, not your narrative. |
| Permission / object-not-found on 265353 | Step 1 hasn't completed — ask Chris to share the narrative. |
| Models gap-fill "ok" but predict zero growth everywhere | Sink/demand reactions were dropped on save. Rebuild with the current `build_2026p.py`, which preserves them. |
| Dashboard shows "Could not load results data" | You opened `index.html` from the filesystem. Serve it over HTTP (above). |

## Getting help

Chris Henry — chenry@anl.gov. For narrative access include your KBase username; for a
pipeline failure include the command you ran and the full traceback.
