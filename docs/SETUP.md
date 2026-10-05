# Setup & run guide

How to install this pipeline and reproduce the results yourself. Written to be
followed start-to-finish, including by a coding agent working on your behalf.

Things you need that are not in this repository: access to the KBase narratives
(Step 1), two data repositories fetched at pinned commits (Step 4), and a KBase auth
token (Step 5). Follow the steps in order. **Do not substitute newer versions of
anything** — every dependency is pinned to the exact commit it was tested with, and
unpinned installs are what broke this guide in September 2026.

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

## Step 2 — The 2022 data is already included

The phenotype data this pipeline reproduces is bundled in `notebooks/data/` — you don't
need to supply anything:

```
notebooks/data/
  BiologMacTestOutput_On-Off.csv    # Clayton's 324 × 519 growth-call matrix (the baseline)
  BiologMacTestOutput.xlsx          # the same 2022 output, source workbook
  Watershed_Mastersheet_Sum22.csv   # the 519-genome roster
  BioLogMacTest.ipynb               # Clayton's original 2022 notebook (the method, reference only)
```

`BiologMacTestOutput_On-Off.csv` is the on/off matrix (media names down the first column,
one column per model named `<genome_id>.fbamodel`, values `0`/`1`) that the pipeline reads
as the **Original 2022p** baseline every "% changed" is measured against.

The Watershed strain-metadata files (sampling locations/dates) are *not* included — they
aren't needed to reproduce the phenotype comparison. You still need the KBase narratives
from Step 1 to recompute the model-based datasets (2022p/2025k/2026p); the bundled data
alone reproduces only the Original 2022p baseline.

## Step 3 — Install

Python **3.11** is required (3.12+ has not been tested against this stack).

Use a **fresh** virtual environment. If you installed this repository before October
2026, delete the old `.venv` and start again rather than upgrading it in place.

```bash
git clone https://github.com/cshenry/WatershedPhenotypeReplication.git
cd WatershedPhenotypeReplication
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` pins every code dependency to an exact commit:

| Package | Source | Commit |
|---|---|---|
| KBUtilLib | `github.com/cshenry/KBUtilLib` | `c89703e0f36d843cad1a384b9388890fca6337a9` (main, 2026-08-01) |
| ModelSEEDpy | `github.com/cshenry/ModelSEEDpy` | `b758b44217bca756a5e07a70e986fbdfe73dd1b6` (main, 2026-07-21) |
| cobrakbase | `github.com/cshenry/cobrakbase` | `68444e46fe3b68482da80798642461af2605e349` (branch `filipe-cobra-model`) |
| cobra | PyPI | `0.30.0` |

Confirm pip installed those commits and not something else:

```bash
pip freeze | grep -iE 'kbutillib|modelseedpy|cobrakbase|^cobra='
```

Each git package should show its SHA from the table. Note that **ModelSEEDpy still
reports itself as version 0.4.2** even at this commit, so `modelseedpy.__version__` and
`pip show` cannot tell it apart from the PyPI 0.4.2 release. Only the SHA in `pip freeze`
can. **Do not** install `modelseedpy` from PyPI, and do not install KBUtilLib or
ModelSEEDpy from a branch (`@main`). The PyPI release predates the API KBUtilLib calls,
and a branch head is a combination nobody has tested.

## Step 4 — Fetch the two data repositories (pinned)

KBUtilLib reads two data repositories from disk. They are not Python packages, so pip
cannot install them. Fetch each one at its pinned commit **into the repository root,
next to `.venv/`**. KBUtilLib finds them by searching the directories above its own
installed location. With the virtual environment at `./.venv` as in Step 3, that search
reaches the repository root. (A virtual environment anywhere else, or a conda
environment, will not find them. Use `./.venv`.) From the repository root:

```bash
# ModelSEED biochemistry database, ~720 MB.
# Do NOT plain-`git clone` it: its default branch (master) was last updated in 2021 and
# is not the biochemistry these models were built with.
git init -q ModelSEEDDatabase
git -C ModelSEEDDatabase remote add origin https://github.com/ModelSEED/ModelSEEDDatabase.git
git -C ModelSEEDDatabase fetch --depth 1 origin e507319a968316f07427977866f99ab2486d349a
git -C ModelSEEDDatabase checkout -q FETCH_HEAD

# Annotation-ontology data (reaction filters and ontology dictionaries), ~100 MB.
git init -q cb_annotation_ontology_api
git -C cb_annotation_ontology_api remote add origin https://github.com/kbaseapps/cb_annotation_ontology_api.git
git -C cb_annotation_ontology_api fetch --depth 1 origin 97f9525aec95390f43f8be88491f35bb11bd589c
git -C cb_annotation_ontology_api checkout -q FETCH_HEAD

# Needed by build_2026p.py; harmless elsewhere. Add it to your shell profile, or export
# it in every new shell before running the pipeline.
export KBUTILLIB_ONTOLOGY_DATA_DIR="$PWD/cb_annotation_ontology_api/data"
```

Both directories are listed in `.gitignore`. Do not configure them through
`~/.kbutillib/dependencies.yaml`: at these KBUtilLib commits that file is read but
never applied.

## Step 5 — KBase token, then verify

Every script authenticates to KBase. Get your token from
[narrative.kbase.us](https://narrative.kbase.us) (Account → Developer Tokens) and make
it available either way:

```bash
export KB_AUTH_TOKEN="<your-token>"
# or write it to the file the loader checks:
mkdir -p ~/.kbase && echo "<your-token>" > ~/.kbase/token
```

Verify the install and your access in one shot. The first line printed checks
Steps 3–4 (no network needed); the second checks Step 1 and your token:

```bash
python -c "
import sys; sys.path.insert(0,'notebooks/PRJ-watershed_phenotype_replication')
from util import get_msfba
fba = get_msfba()
print('install OK: database, ontology data and cobrakbase all loaded')
n = len(fba.list_ws_objects(265353, type='KBaseFBA.FBAModel'))
print(f'OK — {n} models visible in narrative 265353')   # expect 3114
"
```

If both lines print, you are ready. If the first one fails, look up the error in
**Troubleshooting** below. `Invalid token` means the token is wrong or has expired; a
permissions error means Step 1 is not done yet.

---

## Running the pipeline

The stages are **checkpointed and resumable** — every one appends to a TSV and skips
work already marked `ok`, so an interrupted run is safe to restart with the same
command. All are parallel; `--workers 50` suits a large machine, `--workers 2` a laptop.

**Start at stage 3.** The models already exist in 265353. Stages 1 and 2 *save models
into narrative 265353*. With read access they fail when they try to save. With write
access they overwrite the published models. Run them only if you mean to rebuild the
corpus, and agree that with Chris first.

The full order:

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
committing to the full set. Stage 3 only reads from KBase, so smoke-test that one:

```bash
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
| `TypeError: ModelSEEDBiochem.get() got an unexpected keyword argument 'path'` | ModelSEEDpy came from PyPI. Use a fresh venv and `pip install -r requirements.txt` (Step 3), then check `pip freeze` for the SHA. |
| `FileNotFoundError: ... ModelSEEDDatabase/Biochemistry/` | Step 4 is not done, the checkout is not in the repository root, or the virtual environment is not at `./.venv`. |
| `FileNotFoundError: ... cb_annotation_ontology_api/data/FilteredReactions.csv` | Same causes as the row above (Step 4). |
| `Ontology dictionaries (SSO_dictionary.json et al) not found` | `KBUTILLIB_ONTOLOGY_DATA_DIR` not exported in this shell (Step 4). |
| `ModuleNotFoundError: No module named 'cobrakbase'` | Installed from an old `requirements.txt`. Re-run Step 3 in a fresh venv. |
| `Token validation failed ... Invalid token` | The KBase token is wrong or expired. Make a new one (Step 5). |
| `NameError: name 'genome' is not defined` | KBUtilLib too old. Re-run Step 3 in a fresh venv. |
| `Either set callback URL...` | The native ontology path wasn't enabled. The scripts set it themselves; if you call the library directly, set `native_ontology` on the **delegate** (`k.recon._delegate.native_ontology = True`) — the Impl wrapper proxies `__getattr__` but not `__setattr__`, so assigning on the wrapper silently does nothing. |
| `reference cannot be null or the empty string` | A model was saved without a `template_ref`. Rebuild it with `build_2026p.py`. |
| `No object with name Carbon-Pyruvic-Acid` | Media reference not workspace-qualified — it lives in `KBaseMedia`, not your narrative. |
| Permission / object-not-found on 265353 | Step 1 hasn't completed — ask Chris to share the narrative. |
| Models gap-fill "ok" but predict zero growth everywhere | Sink/demand reactions were dropped on save. Rebuild with the current `build_2026p.py`, which preserves them. |
| Dashboard shows "Could not load results data" | You opened `index.html` from the filesystem. Serve it over HTTP (above). |

## Getting help

Chris Henry — chenry@anl.gov. For narrative access include your KBase username; for a
pipeline failure include the command you ran and the full traceback.
