"""WatershedPhenotypeReplication — work-notebook util

Run at the top of every notebook cell with::

    %run util.py

After that line, the cell can use the path constants and ``session`` directly.
Do NOT add per-cell imports.

Path layout (all resolved relative to this file so the notebook works from
any cwd)::

    <repo>/
      notebooks/
        PRJ-watershed_phenotype_replication/
          util.py          <- this file
          NBCache/         <- per-PRJ NotebookSession cache (gitignored)
          NBOutput/        <- per-PRJ outputs (gitignored)
        models/
        genomes/
        data/
"""
from __future__ import annotations

import sys as _sys
import warnings as _warnings
from pathlib import Path as _Path


# === sys.path bootstrap =====================================================
# Prepend machine-specific Python paths before any heavy imports.
# Reads ~/.kbu-sys-paths (one path per line; # comments OK).
# Silent no-op if the file is missing or unreadable.

def _bootstrap_sys_paths() -> None:
    user_file = _Path.home() / ".kbu-sys-paths"
    if not user_file.exists():
        return
    try:
        for raw in user_file.read_text().splitlines():
            s = raw.split("#", 1)[0].strip()
            if not s:
                continue
            expanded = str(_Path(s).expanduser())
            if expanded and expanded not in _sys.path:
                _sys.path.insert(0, expanded)
    except Exception:
        pass


_bootstrap_sys_paths()
# ============================================================================

from pathlib import Path

from kbutillib.notebook import NotebookSession

# ---------------------------------------------------------------------------
# Path constants
# Anchored relative to __file__ so they work regardless of cwd.
#
# PRJ layout:
#   <repo>/notebooks/PRJ-watershed_phenotype_replication/util.py  <- __file__
#   <repo>/notebooks/                           <- NOTEBOOKS_DIR
#   <repo>/                                     <- PROJECT_ROOT
# ---------------------------------------------------------------------------

#: This PRJ directory (contains util.py, NBCache/, NBOutput/).
_PRJ_DIR: Path = Path(__file__).resolve().parent

#: The notebooks/ directory shared across all PRJs.
NOTEBOOKS_DIR: Path = _PRJ_DIR.parent

#: Repository root — one level above notebooks/.
PROJECT_ROOT: Path = NOTEBOOKS_DIR.parent

#: Shared models directory (notebooks/models/).
MODELS_DIR: Path = NOTEBOOKS_DIR / "models"

#: Shared genomes directory (notebooks/genomes/).
GENOMES_DIR: Path = NOTEBOOKS_DIR / "genomes"

#: Shared data directory (notebooks/data/).
DATA_DIR: Path = NOTEBOOKS_DIR / "data"

#: Per-PRJ output directory (gitignored).
NBOUTPUT_DIR: Path = _PRJ_DIR / "NBOutput"

# ---------------------------------------------------------------------------
# NotebookSession
# Cache lands in the PRJ-local NBCache/ directory (gitignored).
# ---------------------------------------------------------------------------

session: NotebookSession = NotebookSession.for_notebook(
    __file__,
    project_name="WatershedPhenotypeReplication",
    cache_dir="NBCache",
)

# === project-specific helpers below ===

import json as _json
import re as _re
import subprocess as _subprocess
import uuid as _uuid

import pandas as pd

# ===========================================================================
# Watershed phenotype-replication study
#
# Goal: reproduce the collaborators' 2022 phenotype-simulation results with
# current tooling, then run the SAME pipeline over the newer model sets so the
# old/new divergence can be attributed to a cause.  The 2022 method is known,
# not inferred -- it is Clayton Piehl's own notebook, notebooks/data/
# BioLogMacTest.ipynb.  See "The 2022 pipeline" section below.
#
# MODEL PROVENANCE -- the old/new split is by OBJECT NAME, not by workspace:
#
#     old  ->  "<genome_id>.fbamodel"                    (2022-09-05; built from
#                                                         PATRIC-imported, RAST-
#                                                         annotated genomes)
#     new  ->  "genomeset__<genome_id>.contigs.fbamodel" (2025-10-03; built from
#                                                         the KBase bulk-RAST
#                                                         re-annotated genome set)
#
# An earlier version of this module claimed the two sets "live in DIFFERENT KBase
# workspaces, so the old/new split is by workspace provenance ... not by parsing
# the object name for a 'genomeset' prefix."  That is FALSE, and it was the exact
# inversion of the truth: NARRATIVE_WS (119455) contains BOTH generations -- 519
# <gid>.fbamodel objects from 2022 AND 524 genomeset__* objects from 2025.
# Splitting by workspace silently labelled the 2025 rebuilds as "old", which
# contaminated the reference arm of the entire study.  Split by name.
# ===========================================================================

#: The Hope College team's genome narrative — source workspace for the "old"
#: (direct genome-model) FBA models.
NARRATIVE_WS: int = 119455

#: Dedicated models workspace holding the "new" (genomeset-rebuilt) FBA models
#: for the Hope genomes.  Numeric workspace ID -- kept as int so
#: list_ws_objects/get_model route it through the ``ids``/``wsid`` branch (a
#: str would be treated as a ws *name*).
MODELS_WS: int = 220871

#: KBase workspace holding the Biolog-style media used as phenotype conditions.
KBASEMEDIA_WS: str = "KBaseMedia"

#: Media used for the baseline carbon growth sanity check.
GLUCOSE_MEDIA: str = "Carbon-D-Glucose"

#: The collaborators' prior phenotype-simulation reference matrix.
WATERSHED_CSV: Path = DATA_DIR / "Watershed_Mastersheet_Sum22.csv"

#: Provenance-correction note shared verbatim between the artifact-inventory
#: notebook's markdown cell and the artifact_inventory.tsv header line -- see
#: the "Artifact inventory" section below for the full explanation.
PROVENANCE_NOTE: str = (
    "PROVENANCE: old/new models are split by OBJECT NAME, not by workspace. "
    "'<genome_id>.fbamodel' (saved 2022-09-05) are the OLD models, built from "
    "PATRIC-imported RAST-annotated genomes. 'genomeset__<genome_id>.contigs."
    "fbamodel' (saved 2025-10-03) are the NEW models, built from the KBase "
    "bulk-RAST re-annotated genome set. NARRATIVE_WS (119455) contains BOTH "
    "generations -- 519 old and 524 new -- so a workspace-based split mislabels "
    "the 2025 rebuilds as 'old'. Two earlier conclusions were wrong and are "
    "superseded here: (a) that only 3/519 genomes have an old model (a "
    "wrong-workspace artifact; the true count is 519/519), and (b) that the "
    "old/new split is by workspace at all. Both old and new models are present "
    "for all 519 reference genomes."
)


def ensure_kbase_token() -> str | None:
    """Ensure KBASE_AUTH_TOKEN is set in the environment, returning the token.

    KBUtilLib's SharedEnvUtils knows how to source the token from
    ``~/.kbase/token`` (and bridges KB_AUTH_TOKEN/KBASE_AUTH_TOKEN), but it
    loads the token into its own token hash, not the environment. FBA/media
    calls elsewhere in this notebook (and the dispatched sweep script) expect
    ``KBASE_AUTH_TOKEN`` to be set, so this reads it via ``get_token('kbase')``
    and exports it. No-op if already set.
    """
    import os

    existing = os.environ.get("KBASE_AUTH_TOKEN")
    if existing:
        return existing
    try:
        from kbutillib.core.shared_env_utils import SharedEnvUtils

        tok = SharedEnvUtils().get_token("kbase")
    except Exception:
        tok = None
    if tok:
        os.environ["KBASE_AUTH_TOKEN"] = tok
    return tok


def get_msfba():
    """Authenticated MSFBAUtils facade (model + media + FBA ops)."""
    tok = ensure_kbase_token()
    from kbutillib.domains.modeling.ms_fba_utils import MSFBAUtils

    return MSFBAUtils(token=tok)


def load_watershed_sheet():
    """Load the Watershed mastersheet into (genome_ids, media_names, value_df).

    Rows are media names (first column → the index); columns are model object
    names of the form ``<genome_id>.fbamodel``.  The ``.fbamodel`` suffix is
    stripped so columns are bare genome IDs.

    Returns:
        genome_ids: list[str] in column order.
        media_names: list[str] in row order.
        value_df: pandas.DataFrame indexed by media, columns = genome IDs,
            values coerced to float (NaN where unparseable).
    """
    df = pd.read_csv(WATERSHED_CSV, index_col=0)
    df.index = [str(i).strip() for i in df.index]
    df = df[[c for c in df.columns if not str(c).startswith("Unnamed")]]
    df.columns = [
        str(c)[: -len(".fbamodel")] if str(c).endswith(".fbamodel") else str(c)
        for c in df.columns
    ]
    df = df.apply(pd.to_numeric, errors="coerce")
    return list(df.columns), list(df.index), df


#: PATRIC-style genome id token (e.g. "562.55367") embedded in a model name.
_PATRIC_ID = _re.compile(r"\d+\.\d+")


def _extract_genome_id(name, ref_set):
    """Pull the genome id from a (possibly nonstandard) model name.

    Returns the first ``\\d+.\\d+`` token in the name that is a reference
    genome; failing that, the first such token; ``None`` if there is none.
    This handles names like ``562.55510Old2022.fbamodel`` (-> 562.55510),
    ``562.85951set_562.85951.fbamodel`` (-> 562.85951), and
    ``genomeset_5_WS_562.55367.fa_assembly.fbamodel`` (-> 562.55367).
    """
    toks = _PATRIC_ID.findall(name)
    ref_hits = [t for t in toks if t in ref_set]
    if ref_hits:
        return ref_hits[0]
    return toks[0] if toks else None


def _latest_ref(refs):
    """Pick the highest-version ``wsid/objid/version`` ref from a list.

    Used wherever a single model reference is needed per genome (e.g. before
    dispatching a sweep) but multiple saved versions/variants exist.
    """
    return max(refs, key=lambda r: int(r.rsplit("/", 1)[-1]))


#: Marks a model rebuilt from the KBase genome-SET pipeline (i.e. from the
#: re-annotated genomes) rather than built in 2022 from the PATRIC import.
#:
#: Matched CASE-INSENSITIVELY and ANYWHERE in the name, because the naming is not
#: consistent.  A ``genomeset__`` prefix test is too narrow and silently mislabels
#: 13 models as "old" -- putting 2025 rebuilds into the 2022 reference arm, which is
#: exactly the contamination this split exists to prevent.  The real names include:
#:
#:     genomeset__<gid>.contigs.fbamodel              (519 + 519)
#:     genomeset__<gid>.fa_assembly.fbamodel          (307)
#:     genomeset_WS5_10-25_<gid>.fa_assembly.fbamodel (5)   <- single underscore
#:     genomeset_5_WS_<gid>.fa_assembly.fbamodel      (5)   <- single underscore
#:     GenomeSetATCC_<gid>.fa_assembly.fbamodel       (3)   <- different case
#:
#: Old models are the plain ``<gid>.fbamodel`` objects (saved 2022-09-05).
GENOMESET_PATTERN = _re.compile(r"genomeset", _re.IGNORECASE)


def is_genomeset_model(name: str) -> bool:
    """True if a model object name marks it as a genome-set ("new") rebuild."""
    return bool(GENOMESET_PATTERN.search(str(name)))


#: Models saved on/after this date are NOT the 2022 build, whatever they are called.
#:
#: The NAME ALONE IS NOT SUFFICIENT, and trusting it puts 2025 models into the 2022
#: reference arm.  Three real counter-examples:
#:
#:   562.61170   "562.61170.fbamodel"        119455/1656/1  saved 2022-06-17  <- true old
#:               "562.61170.fbamodel"        220871/5579/3  saved 2025-11-13  <- SAME NAME, 2025
#:               _latest_ref picks the highest VERSION -> v3 -> the 2025 model. Wrong arm.
#:
#:   562.55510   "562.55510Old2022.fbamodel" 220871/5597/1  saved 2026-06-15
#:               named "Old2022", saved 2026.  The name lies.
#:
#:   562.85951   "562.85951.fbamodel"        119455/2080/2  saved 2025-10-03
#:               plain <gid>.fbamodel in the narrative -- but built in 2025, not 2022.
#:
#: The save date is the only reliable discriminator.  Use it.
OLD_MODEL_CUTOFF: str = "2023-01-01"


def _save_date(info) -> str:
    """ISO save date from a list_ws_objects info tuple (index 3)."""
    return str(info[3])[:10] if len(info) > 3 else ""


def _model_generation(name, info, cutoff=OLD_MODEL_CUTOFF):
    """Classify a model object as "old" (2022) or "new" (genomeset / post-cutoff).

    A model is OLD only if BOTH hold:
      * its name carries no genome-set marker, AND
      * it was saved before ``cutoff``.

    Either signal alone is insufficient -- see OLD_MODEL_CUTOFF for the three real
    objects that break a name-only rule.
    """
    if is_genomeset_model(name):
        return "new"
    saved = _save_date(info)
    if saved and saved >= cutoff:
        return "new"
    return "old"


def partition_models(ref_genomes, narrative_ws=NARRATIVE_WS, models_ws=MODELS_WS, fba=None):
    """Partition FBA models into old/new by OBJECT-NAME provenance.

    Provenance correction (2026-07-12) -- the old/new split is by object *name*,
    NOT by workspace:

        old  ->  "<genome_id>.fbamodel"                  (built 2022-09-05 from
                                                          PATRIC-imported, RAST-
                                                          annotated genomes)
        new  ->  "genomeset__<genome_id>.contigs.fbamodel" (built 2025-10-03 from
                                                          the KBase bulk-RAST
                                                          re-annotated genome set)

    An earlier version of this function split on workspace (narrative=old,
    models=new).  That is WRONG: the narrative workspace (119455) contains BOTH
    generations -- 519 ``<gid>.fbamodel`` objects saved 2022-09-05 *and* 524
    ``genomeset__*`` objects saved 2025-10-03.  Splitting by workspace silently
    labelled the 2025 genomeset models as "old", contaminating the reference arm
    of the replication study.  Both tracked workspaces are scanned here and every
    model is classified by its name.

    Where a genome has multiple saved refs for the same class (repeat saves /
    versions), all refs are collected; use ``_latest_ref(refs)`` to pick one.

    Args:
        ref_genomes: iterable of reference genome ids to keep.
        narrative_ws: Hope College genome narrative (holds both generations).
        models_ws: dedicated models workspace (holds "new" models only).
        fba: optional MSFBAUtils facade (defaults to ``get_msfba()``).

    Returns:
        (partitions, report) where:
          partitions: {"old": {genome_id: [ref, ...]},
                       "new": {genome_id: [ref, ...]}}
          report: {"excluded": {"old": [...], "new": [...]},
                   "unmappable": {"old": [...], "new": [...]}}
    """
    fba = fba or get_msfba()
    ref_set = set(str(g) for g in ref_genomes)
    partitions = {"old": {}, "new": {}}
    report = {
        "excluded": {"old": [], "new": []},
        "unmappable": {"old": [], "new": []},
    }
    saved_dates: dict = {}

    for wsid in (narrative_ws, models_ws):
        objs = fba.list_ws_objects(wsid, type="KBaseFBA.FBAModel")
        for name, info in objs.items():
            name = str(name)
            label = _model_generation(name, info)
            gid = _extract_genome_id(name, ref_set)
            if gid is None:
                report["unmappable"][label].append(name)
                continue
            if gid not in ref_set:
                report["excluded"][label].append(name)
                continue
            objid, version, obj_wsid = info[0], info[4], info[6]
            ref = f"{obj_wsid}/{objid}/{version}"
            partitions[label].setdefault(gid, []).append(ref)
            saved_dates[(label, gid, ref)] = _save_date(info)

    # Collapse each genome to ONE ref per arm, deterministically.
    #
    # _latest_ref (highest VERSION) is wrong for the old arm: 562.61170 has a 2022 model
    # at v1 and a 2025 model at v3 under the SAME NAME, and version-ordering picks the
    # 2025 one. For "old" take the EARLIEST-saved object (the actual 2022 build); for
    # "new" take the latest.
    for label in ("old", "new"):
        for gid, refs in partitions[label].items():
            if len(refs) <= 1:
                continue
            dated = [(saved_dates.get((label, gid, r), ""), r) for r in refs]
            dated.sort()
            chosen = dated[0][1] if label == "old" else dated[-1][1]
            report.setdefault("ambiguous", {}).setdefault(label, []).append(
                {"genome_id": gid, "chosen": chosen,
                 "candidates": [{"ref": r, "saved": d} for d, r in dated]}
            )
            partitions[label][gid] = [chosen]

    return partitions, report



# ---------------------------------------------------------------------------
# Artifact inventory
#
# Provenance-correction note: an earlier pass concluded "only 3/519 genomes
# have an old model" by listing KBaseFBA.FBAModel objects in the WRONG
# workspace (MODELS_WS / 220871, which only holds the *new* genomeset-rebuilt
# models). The "old" models were built directly against the Hope College
# genome narrative and live in NARRATIVE_WS (119455). census_artifacts()
# walks all involved workspaces and splits old/new by workspace provenance
# (see partition_models), so this class of wrong-workspace undercount cannot
# recur silently -- the counts below now come from the correct workspace per
# artifact class.
# ---------------------------------------------------------------------------

#: KBase type (module.Type, version suffix stripped) -> artifact class, for
#: the classes whose type alone determines the class (i.e. everything except
#: KBaseFBA.FBAModel, which additionally depends on source workspace --
#: see _classify_artifact).
_TYPE_ARTIFACT_CLASS: dict[str, str] = {
    "KBaseGenomes.Genome": "genome",
    "KBasePhenotypes.PhenotypeSet": "observed_phenotype_set",
    "KBasePhenotypes.PhenotypeSimulationSet": "prior_simulation_set",
    "KBaseBiochem.Media": "media",
    "KBaseGenomeAnnotations.Assembly": "assembly",
}

#: Artifact classes that are normalized to (and expected once per) a single
#: reference genome.
GENOME_SCOPED_CLASSES: tuple[str, ...] = (
    "genome",
    "old_model",
    "new_model",
    "observed_phenotype_set",
    "prior_simulation_set",
    "assembly",
)

#: Artifact classes that are shared across all genomes (not genome-scoped).
#: Media definitions (Biolog-style conditions) live once in KBaseMedia and
#: apply to every reference genome, so they are censused as a pool rather
#: than normalized per genome.
GLOBAL_CLASSES: tuple[str, ...] = ("media",)

#: All artifact classes tracked by census_artifacts, in a stable order.
ARTIFACT_CLASSES: tuple[str, ...] = GENOME_SCOPED_CLASSES + GLOBAL_CLASSES

#: Column order for census_artifacts' table_df.
TABLE_COLUMNS: list[str] = [
    "genome_id",
    "artifact_class",
    "workspace",
    "object_name",
    "ref",
    "version",
    "save_date",
]

#: Metadata keys checked (in order) for an explicit genome linkage before
#: falling back to reference-token extraction from the object name.
_GENOME_REF_META_KEYS: tuple[str, ...] = (
    "genome_ref",
    "Genome ref",
    "genome_id",
    "Genome ID",
)


def _classify_artifact(base_type, wsid, narrative_ws=NARRATIVE_WS, models_ws=MODELS_WS,
                       name="", save_date=""):
    """Map a (version-stripped) KBase type to an artifact class.

    KBaseFBA.FBAModel splits into "old_model" / "new_model" by OBJECT NAME --
    the ``genomeset__`` prefix marks the 2025 re-annotated rebuild -- not by
    workspace, because the narrative workspace holds both generations (see
    partition_models).  An FBAModel in neither tracked workspace is not
    classified (returns ``None``).  All other types map 1:1 via
    ``_TYPE_ARTIFACT_CLASS``.
    """
    if base_type == "KBaseFBA.FBAModel":
        if wsid not in (narrative_ws, models_ws):
            return None
        # Same rule as partition_models: name AND save date. Neither alone is enough
        # -- see OLD_MODEL_CUTOFF for the objects that break a name-only test.
        gen = _model_generation(name, (None, None, None, save_date))
        return "new_model" if gen == "new" else "old_model"
    return _TYPE_ARTIFACT_CLASS.get(base_type)


def _normalize_genome_id(name, meta, ref_set):
    """Normalize an artifact to a reference genome id.

    Preference order:
      1. Object metadata (``genome_ref`` / ``Genome ref`` / ``genome_id`` /
         ``Genome ID``), tokenized for reference-genome ids the same way as
         ``_extract_genome_id``.
      2. Reference-token-first extraction from the object name.

    Unlike ``_extract_genome_id`` (which silently takes the first reference
    hit), this reports ambiguity: if more than one *distinct* reference
    genome id is found in the chosen source, the artifact is unmappable
    rather than guessed.

    Args:
        name: object name.
        meta: object metadata dict (may be ``None``/empty).
        ref_set: set of reference genome id strings.

    Returns:
        (genome_id, reason) where genome_id is ``None`` on failure and
        reason is one of "metadata", "name", "ambiguous", "no_match".
    """
    meta = meta or {}
    for key in _GENOME_REF_META_KEYS:
        val = meta.get(key)
        if not val:
            continue
        toks = _PATRIC_ID.findall(str(val))
        ref_hits = sorted(set(t for t in toks if t in ref_set))
        if len(ref_hits) == 1:
            return ref_hits[0], "metadata"
        if len(ref_hits) > 1:
            return None, "ambiguous"

    toks = _PATRIC_ID.findall(str(name))
    ref_hits = sorted(set(t for t in toks if t in ref_set))
    if len(ref_hits) == 1:
        return ref_hits[0], "name"
    if len(ref_hits) > 1:
        return None, "ambiguous"
    return None, "no_match"


def census_artifacts(
    workspaces,
    ref_genomes=None,
    fba=None,
    narrative_ws=NARRATIVE_WS,
    models_ws=MODELS_WS,
):
    """Census KBase artifacts across ``workspaces``, normalized to reference genomes.

    Walks every workspace in ``workspaces`` via ``MSFBAUtils.list_ws_objects``
    (typically ``[NARRATIVE_WS, MODELS_WS, KBASEMEDIA_WS]``), classifies each
    object by KBase type (see ``_classify_artifact``) into one of
    ``ARTIFACT_CLASSES``, and normalizes genome-scoped artifacts to a
    reference genome id via ``_normalize_genome_id``. Objects whose type is
    not tracked are silently skipped; objects of a tracked type that cannot
    be normalized to a reference genome are recorded in
    ``missing_report["unmappable"]`` and excluded from ``table_df``.

    See the module-level provenance-correction note: this function replaces
    an earlier single-workspace FBAModel scan that undercounted "old"
    models by looking in the wrong workspace.

    Args:
        workspaces: iterable of workspace ids/names to scan (passed through
            to ``MSFBAUtils.list_ws_objects``).
        ref_genomes: iterable of reference genome ids. Defaults to the
            genome id columns of the Watershed mastersheet
            (``load_watershed_sheet()``).
        fba: optional MSFBAUtils facade (defaults to ``get_msfba()``).
        narrative_ws: workspace id for "old" (direct genome) FBA models.
        models_ws: workspace id for "new" (genomeset-rebuilt) FBA models.

    Returns:
        (table_df, missing_report):
          table_df: pandas.DataFrame with columns
              ``["genome_id", "artifact_class", "workspace", "object_name",
              "ref", "version", "save_date"]`` -- one row per censused
              artifact. ``genome_id`` is ``None`` for GLOBAL_CLASSES
              (media) rows, since those artifacts are shared, not
              per-genome.
          missing_report: {
              "per_genome": {genome_id: [missing artifact_class, ...]},
                  -- GENOME_SCOPED_CLASSES with no row for that genome;
                  "media" is appended to every genome's list if and only if
                  zero media objects were found anywhere (media is a pool,
                  not per-genome, so a single found media object clears it
                  for all genomes).
              "counts": {artifact_class: total objects found},
              "missing_counts": {artifact_class: number of genomes missing it},
              "unmappable": [{"workspace", "object_name", "type", "reason"}, ...],
              "media_available": bool,
          }
    """
    if ref_genomes is None:
        ref_genomes, _unused_media, _unused_df = load_watershed_sheet()
    ref_set = set(str(g) for g in ref_genomes)
    fba = fba or get_msfba()

    rows = []
    unmappable = []
    counts = {c: 0 for c in ARTIFACT_CLASSES}

    for ws in workspaces:
        objs = fba.list_ws_objects(ws)
        for name, info in objs.items():
            name = str(name)
            full_type = info[2]
            base_type = full_type.split("-", 1)[0]
            obj_wsid = info[6]
            objid, version, save_date = info[0], info[4], info[3]

            artifact_class = _classify_artifact(
                base_type, obj_wsid, narrative_ws, models_ws,
                name=name, save_date=_save_date(info),
            )
            if artifact_class is None:
                continue  # not a tracked artifact type

            meta = info[10] if len(info) > 10 else {}
            ref = f"{obj_wsid}/{objid}/{version}"

            if artifact_class in GLOBAL_CLASSES:
                counts[artifact_class] += 1
                rows.append(
                    {
                        "genome_id": None,
                        "artifact_class": artifact_class,
                        "workspace": obj_wsid,
                        "object_name": name,
                        "ref": ref,
                        "version": version,
                        "save_date": save_date,
                    }
                )
                continue

            gid, reason = _normalize_genome_id(name, meta, ref_set)
            if gid is None:
                unmappable.append(
                    {
                        "workspace": obj_wsid,
                        "object_name": name,
                        "type": base_type,
                        "reason": reason,
                    }
                )
                continue

            counts[artifact_class] += 1
            rows.append(
                {
                    "genome_id": gid,
                    "artifact_class": artifact_class,
                    "workspace": obj_wsid,
                    "object_name": name,
                    "ref": ref,
                    "version": version,
                    "save_date": save_date,
                }
            )

    table_df = pd.DataFrame(rows, columns=TABLE_COLUMNS)

    present = set()
    for gid, cls in zip(table_df["genome_id"], table_df["artifact_class"]):
        if cls in GENOME_SCOPED_CLASSES:
            present.add((gid, cls))

    media_available = counts["media"] > 0

    per_genome = {}
    missing_counts = {c: 0 for c in ARTIFACT_CLASSES}
    for gid in sorted(ref_set):
        missing = [c for c in GENOME_SCOPED_CLASSES if (gid, c) not in present]
        if not media_available:
            missing.append("media")
        for c in missing:
            missing_counts[c] += 1
        per_genome[gid] = missing

    missing_report = {
        "per_genome": per_genome,
        "counts": counts,
        "missing_counts": missing_counts,
        "unmappable": unmappable,
        "media_available": media_available,
    }

    return table_df, missing_report

# ===========================================================================
# Compare — mastersheet scoring (SUPERSEDED; retained + tested, not on the main path)
#
# NOTE: the study's scoring now runs through score_calls() against Clayton's own
# 2022 calls and the experimental Biolog data (see "The 2022 pipeline" section).
# The helpers below score against the mastersheet the earlier free-transport sweep
# produced. They are pure, still unit-tested (tests/test_compare.py), and kept for
# reference, but nothing in the current notebooks or scripts calls them.
#
# Scores long-format simulated growth calls against two independent references:
#   - the collaborators' prior mastersheet (Watershed_Mastersheet_Sum22.csv,
#     via load_watershed_sheet -- a WIDE media x genome_id matrix of growth
#     rates), and
#   - the narrative's observed KBasePhenotypes.PhenotypeSet objects (the
#     actual Biolog-style experimental calls the mastersheet and the models
#     are both trying to replicate).
#
# ``score_against_reference`` is the general-purpose long-format scorer used
# for both; ``integrity_check`` is a narrower sanity check that the
# mastersheet matches the narrative's own KBasePhenotypes.PhenotypeSimulationSet
# ("..._SimulatedGrowth") objects -- i.e. that the collaborators' spreadsheet
# is a faithful export of what they simulated in KBase, not a transcription
# drift.
#
# ``score_against_reference`` and ``integrity_check`` are pure (no KBase
# calls) and fully unit-tested offline (tests/test_compare.py). The
# ``resolve_media`` / ``fetch_observed_phenotypes`` /
# ``fetch_narrative_simulated_growth`` helpers below DO make KBase calls
# (list_ws_objects / get_object / get_object_info) -- they exist to feed the
# notebook and are exercised in tests only via an in-memory FakeFBA stub
# (never a live session), matching the pattern already used by
# partition_models / census_artifacts.
# ===========================================================================


def score_against_reference(
    pred_df,
    ref_df,
    *,
    growth_threshold=0.01,
    observed_threshold=1e-6,
    model_col="genome_id",
    media_col="media_name",
    pred_value_col="value",
    ref_value_col="value",
):
    """Score predicted growth values against a reference, per-model and aggregate.

    Both inputs are long-format tables sharing a (``model_col``, ``media_col``)
    key -- e.g. ``pred_df`` from a filtered slice of the sweep results TSV
    (one ``model_set``/``free_transport`` combo) and ``ref_df`` from either
    the melted mastersheet or ``fetch_observed_phenotypes``.

    Binarization (independent per side, both growth-inclusive at the
    threshold -- i.e. a value exactly equal to its threshold counts as
    growth):
        pred_growth = pred_value >= growth_threshold
        ref_growth  = ref_value  >= observed_threshold
    A NaN value on either side binarizes to ``False`` (no growth) via normal
    NaN-comparison semantics (``NaN >= x`` is always ``False`` in numpy/
    pandas) -- this implements "NaN prediction -> no-growth" without special
    casing. The one explicit exception: a cell that is NaN on **both** sides
    (no prediction row and no reference row/value for that
    (model, media) pair) is excluded entirely from CP/CN/FP/FN/accuracy
    rather than silently counted as a correct negative.

    The continuous correlation (Pearson + Spearman, on the raw un-binarized
    values) is computed separately and only over cells with a real
    (non-NaN) value on *both* sides.

    Args:
        pred_df: long-format DataFrame with at least
            [model_col, media_col, pred_value_col].
        ref_df: long-format DataFrame with at least
            [model_col, media_col, ref_value_col].
        growth_threshold: prediction-side binarization cutoff (inclusive).
        observed_threshold: reference-side binarization cutoff (inclusive).
        model_col: grouping key for the "per-model" breakdown (default
            "genome_id" -- each row is one model/genome's call on one media).
        media_col: the media identity column shared by both inputs.
        pred_value_col: value column name in ``pred_df``.
        ref_value_col: value column name in ``ref_df``.

    Returns:
        dict with:
          "merged": the full outer-joined long DataFrame, one row per
              (model, media) pair seen in either input, with pred_value,
              ref_value, pred_growth, ref_growth, class (CP/CN/FP/FN or
              None if excluded).
          "per_model": DataFrame indexed by ``model_col`` with columns
              [CP, CN, FP, FN, n, accuracy, pearson_r, pearson_p,
              spearman_r, spearman_p].
          "aggregate": same columns as a single dict, computed across all
              models combined.
    """
    p = pred_df[[model_col, media_col, pred_value_col]].rename(
        columns={pred_value_col: "pred_value"}
    )
    r = ref_df[[model_col, media_col, ref_value_col]].rename(
        columns={ref_value_col: "ref_value"}
    )
    merged = p.merge(r, on=[model_col, media_col], how="outer")

    pred_vals = pd.to_numeric(merged["pred_value"], errors="coerce")
    ref_vals = pd.to_numeric(merged["ref_value"], errors="coerce")
    merged["pred_value"] = pred_vals
    merged["ref_value"] = ref_vals

    both_nan = pred_vals.isna() & ref_vals.isna()
    # NaN >= threshold is False -- this is the "NaN -> no-growth" rule for
    # both sides, no special-casing needed beyond the both-NaN exclusion.
    merged["pred_growth"] = pred_vals >= growth_threshold
    merged["ref_growth"] = ref_vals >= observed_threshold

    def _classify(pred_growth, ref_growth, excluded):
        if excluded:
            return None
        if ref_growth and pred_growth:
            return "CP"
        if (not ref_growth) and (not pred_growth):
            return "CN"
        if (not ref_growth) and pred_growth:
            return "FP"
        return "FN"

    merged["class"] = [
        _classify(pg, rg, ex)
        for pg, rg, ex in zip(merged["pred_growth"], merged["ref_growth"], both_nan)
    ]

    def _summarize(df):
        scored = df[df["class"].notna()]
        counts = scored["class"].value_counts()
        CP = int(counts.get("CP", 0))
        CN = int(counts.get("CN", 0))
        FP = int(counts.get("FP", 0))
        FN = int(counts.get("FN", 0))
        n = CP + CN + FP + FN
        accuracy = (CP + CN) / n if n else float("nan")

        corr_df = df.dropna(subset=["pred_value", "ref_value"])
        if (
            len(corr_df) >= 2
            and corr_df["pred_value"].nunique() > 1
            and corr_df["ref_value"].nunique() > 1
        ):
            from scipy import stats

            pearson_r, pearson_p = stats.pearsonr(corr_df["pred_value"], corr_df["ref_value"])
            spearman_r, spearman_p = stats.spearmanr(corr_df["pred_value"], corr_df["ref_value"])
        else:
            pearson_r = pearson_p = spearman_r = spearman_p = float("nan")

        return {
            "CP": CP, "CN": CN, "FP": FP, "FN": FN, "n": n, "accuracy": accuracy,
            "pearson_r": float(pearson_r), "pearson_p": float(pearson_p),
            "spearman_r": float(spearman_r), "spearman_p": float(spearman_p),
        }

    aggregate = _summarize(merged)

    per_model_rows = []
    for key, grp in merged.groupby(model_col, dropna=False):
        row = _summarize(grp)
        row[model_col] = key
        per_model_rows.append(row)
    per_model_cols = [
        model_col, "CP", "CN", "FP", "FN", "n", "accuracy",
        "pearson_r", "pearson_p", "spearman_r", "spearman_p",
    ]
    per_model_df = pd.DataFrame(per_model_rows, columns=per_model_cols)
    if not per_model_df.empty:
        per_model_df = per_model_df.set_index(model_col)

    return {"merged": merged, "per_model": per_model_df, "aggregate": aggregate}


def integrity_check(
    mastersheet_df,
    simset_df,
    *,
    tolerance=1e-6,
    genome_col="genome_id",
    media_col="media_name",
    value_col="value",
    version_col="version",
):
    """Check the mastersheet against the narrative's own simulated-growth objects.

    ``mastersheet_df`` is the WIDE growth-rate matrix returned by
    ``load_watershed_sheet()`` (index = media name, columns = genome id).
    ``simset_df`` is the LONG table returned by
    ``fetch_narrative_simulated_growth()`` (one row per (genome, media,
    saved version) triple, from the narrative's
    KBasePhenotypes.PhenotypeSimulationSet "..._SimulatedGrowth" objects).

    Joins on (genome_id, media_name); when ``simset_df`` has multiple saved
    ``version`` rows for the same (genome_id, media_name), only the highest
    version is kept (the mastersheet should reflect the *latest* simulation,
    not a stale one). Flags cells where
    ``abs(mastersheet_value - simset_value) > tolerance``.

    Args:
        mastersheet_df: wide DataFrame, index=media_name, columns=genome_id.
        simset_df: long DataFrame with [genome_col, media_col, value_col,
            version_col] (version_col optional -- if absent, every row is
            treated as already deduplicated).
        tolerance: absolute-difference flag threshold (default 1e-6).

    Returns:
        (compared_df, flagged_df):
          compared_df: long DataFrame with one row per (genome_id,
              media_name) pair present in BOTH inputs, columns
              [genome_col, media_col, "mastersheet_value", "simset_value",
              "abs_diff", "flagged"].
          flagged_df: the subset of compared_df where "flagged" is True.
    """
    master_long = (
        mastersheet_df.rename_axis(media_col)
        .reset_index()
        .melt(id_vars=media_col, var_name=genome_col, value_name="mastersheet_value")
    )

    sim = simset_df.rename(columns={value_col: "simset_value"}).copy()
    if version_col in sim.columns:
        sim = sim.sort_values(version_col).drop_duplicates(
            subset=[genome_col, media_col], keep="last"
        )

    compared = master_long.merge(
        sim[[genome_col, media_col, "simset_value"]],
        on=[genome_col, media_col],
        how="inner",
    )
    compared["mastersheet_value"] = pd.to_numeric(compared["mastersheet_value"], errors="coerce")
    compared["simset_value"] = pd.to_numeric(compared["simset_value"], errors="coerce")
    compared["abs_diff"] = (compared["mastersheet_value"] - compared["simset_value"]).abs()
    compared["flagged"] = compared["abs_diff"] > tolerance

    flagged_df = compared[compared["flagged"]].reset_index(drop=True)
    return compared, flagged_df


def resolve_media(media_names, fba=None, media_ws=KBASEMEDIA_WS):
    """Resolve mastersheet media row names against the KBaseMedia workspace.

    Args:
        media_names: iterable of media row names (e.g. from
            ``load_watershed_sheet()``'s ``media_names``).
        fba: optional MSFBAUtils facade (defaults to ``get_msfba()``).
        media_ws: workspace holding the media objects (default KBaseMedia).

    Returns:
        (resolved, unresolved) where:
          resolved: {media_name: media object data} for every name found by
              exact match against object names in ``media_ws``.
          unresolved: sorted list of media_names with no match -- flag these
              rather than silently dropping them from the sweep.
    """
    fba = fba or get_msfba()
    objs = fba.list_ws_objects(media_ws, type="KBaseBiochem.Media")

    resolved = {}
    unresolved = []
    for media_name in media_names:
        media_name = str(media_name)
        info = objs.get(media_name)
        if info is None:
            unresolved.append(media_name)
            continue
        objid, version, obj_wsid = info[0], info[4], info[6]
        ref = f"{obj_wsid}/{objid}/{version}"
        resolved[media_name] = fba.get_object(ref)

    return resolved, sorted(unresolved)


# ===========================================================================
# The 2022 pipeline — reconstructed from Clayton Piehl's own notebook
#
# Source of truth: notebooks/data/BioLogMacTest.ipynb (supplied 2026-07-13).
# This is the code that actually produced the mastersheet, so it -- not our
# inference -- defines the method:
#
#     gmm      = get_media("Carbon-Pyruvic-Acid", "KBaseMedia")
#     template = get_template("GramNegModelTemplateV4", "NewKBaseModelTemplates")
#     biolog   = get_phenotypeset("ecoli_biolog", 93541)     # 324 phenotypes
#
#     for model in old_models:                               # ws 119455
#         model.objective = "bio1"
#         gf  = MSGapfill(model, [template], [], {}, {}, [])
#         sol = gf.run_gapfilling(gmm, "bio1")
#         gf.integrate_gapfill_solution(sol)
#         simulate(gf.mdlutl, biolog, add_missing_exchanges=True)
#
# They GAPFILLED.  They did not supplement the medium.  An earlier version of
# this module inferred a media-supplement mechanism and reached ~85% by
# approximating what gapfilling does; it was wrong and has been deleted rather
# than kept alongside, so that nobody can pick it up by mistake.
#
# The Watershed mastersheet is the OLD MODELS' OUTPUT (Clayton, 2026-06-23:
# "Here is the csv output from the old models from Summer of 2022"), not
# experimental data.  Reproducing it is REPLICATION FIDELITY, never "accuracy".
# The scoring target is his own binarised call matrix
# (notebooks/data/BiologMacTestOutput_On-Off.csv), which we verified reproduces
# the mastersheet exactly: 0 disagreements across 168,156 cells, with growth
# calls above flux 0.048 and no-growth below 1e-10 -- a clean gap, so any cutoff
# in [1e-9, 0.01] recovers his calls.
#
# ---------------------------------------------------------------------------
# THREE DEFECTS IN THE UPSTREAM TOOLS THAT THIS MODULE WORKS AROUND
# ---------------------------------------------------------------------------
#
# 1. CO2 MAKES E. COLI AUTOTROPHIC.  The KBase media declare CO2 (cpd00011) with
#    maxFlux=0 -- NO uptake.  MSGrowthPhenotype.build_media() discards that and
#    rebuilds every compound at -100, CO2 included.  Gapfilling then closes a
#    carbon-fixation route in some models, and they grow on CO2 + minerals with
#    NO organic carbon at all (562.55367: biomass 0.727).  E. coli is not
#    autotrophic.  Those models then "grow" on every medium carrying CO2 --
#    gelatin, inulin, adonitol, 2,3-butanediol, 56 substrates in total -- while
#    glucose (whose medium has no CO2) stays correct, so nothing looks broken.
#    Measured across 8 genomes: blocking CO2 uptake removes ALL 235 false
#    positives and lifts agreement 0.909 -> 0.965.
#    => block_co2_uptake() restores the bound the medium already specified.
#    => assert_not_autotrophic() is the guard that would have caught it.
#
# 2. GAPFILL SOLUTIONS LEAK BETWEEN MODELS.  MSGapfill.integrate_gapfill_solution
#    declares `cumulative_solution=[]` -- a MUTABLE DEFAULT, shared across every
#    call in the process.  Model N's solution gets contaminated with reactions
#    from models 1..N-1, and eventually raises KeyError when a leaked reaction is
#    looked up in a model that never had it.  The crash is the lucky case; the
#    silent corruption is the dangerous one.
#    => ALWAYS pass cumulative_solution=[] explicitly.  gapfill_model() does.
#
# 3. integrate_gapfill_solution RETURNS A DICT, NOT A MODEL.  It mutates
#    gf.mdlutl in place.  Clayton's 2022 line `model = msgapfill.integrate_...`
#    is an older API and now silently rebinds `model` to a dict -- which is why
#    his notebook cannot be re-run as-is on current modelseedpy.
#    => use gf.mdlutl.
#
# See agent-io/research/2026-07-13-*.md for reproductions.
# ===========================================================================

#: Media the 2022 run gapfilled against, before simulating phenotypes.
GAPFILL_MEDIA: str = "Carbon-Pyruvic-Acid"

#: Reconstruction template used for gapfilling.
GAPFILL_TEMPLATE: str = "GramNegModelTemplateV4"
GAPFILL_TEMPLATE_WS: str = "NewKBaseModelTemplates"

#: The Biolog phenotype set. Its 324 phenotype ids match the mastersheet's 324
#: media rows exactly (verified: zero difference either way).
BIOLOG_PHENOTYPESET: str = "ecoli_biolog"
BIOLOG_PHENOTYPESET_WS: int = 93541

#: CO2 -- must never be importable. See defect 1 above.
CO2_COMPOUND: str = "cpd00011"

#: Simulated biomass flux above which we call growth.
GROWTH_THRESHOLD: float = 0.01

#: Mastersheet value above which the 2022 run called growth.  Clayton's growth
#: calls sit above 0.048 and his no-growth below 1e-10, so anything in
#: [1e-9, 0.01] reproduces his binarisation exactly.
OBSERVED_THRESHOLD: float = 1e-6


def get_biolog_media(fba=None, block_co2=True):
    """The 324 Biolog phenotype media, keyed by media name.

    Args:
        fba: optional MSFBAUtils facade.
        block_co2: close CO2 uptake, honouring the ``maxFlux=0`` the KBase media
            already declare.  Leave this True unless you are deliberately
            reproducing the autotrophy bug -- with CO2 open, gapfilled models can
            fix carbon and report growth on substrates they cannot eat.

    Returns:
        {media_name: MSMedia}
    """
    fba = fba or get_msfba()
    phenoset = fba.get_phenotypeset(BIOLOG_PHENOTYPESET, BIOLOG_PHENOTYPESET_WS)
    media = {}
    for pheno in phenoset.phenotypes:
        medium = pheno.build_media()
        if block_co2:
            block_co2_uptake(medium)
        media[pheno.id] = medium
    return media


def block_co2_uptake(medium, compound=CO2_COMPOUND):
    """Close CO2 uptake on a medium, leaving excretion open.

    ``lower_bound = 0`` forbids import while still permitting the model to emit
    CO2 -- which is exactly what the KBase media mean by ``maxFlux = 0``, and what
    MSGrowthPhenotype.build_media() throws away.  Without this, a gapfilled model
    that acquired a carbon-fixation route grows on any CO2-bearing medium.
    """
    for cpd in medium.mediacompounds:
        if cpd.id == compound:
            cpd.lower_bound = 0.0
    return medium


def gapfill_model(model, fba=None, gapfill_media=GAPFILL_MEDIA,
                  template=GAPFILL_TEMPLATE, template_ws=GAPFILL_TEMPLATE_WS,
                  biomass="bio1", check_autotrophy=True, strict=False,
                  media_obj=None, template_obj=None):
    """Gapfill one model exactly as the 2022 pipeline did, and verify the result.

    The saved 2022 models are auxotrophic: simulated on a bare Biolog medium every
    one of them returns 0.0 biomass, glucose included.  Gapfilling against
    ``Carbon-Pyruvic-Acid`` is what makes them viable -- it is not optional, and it
    is the step that was missing from every earlier reconstruction attempt.

    Args:
        model: a cobra.Model / MSModelUtil.  A PRIVATE COPY is taken: ``get_model``
            caches, and ``add_missing_exchanges`` mutates in place, so simulating a
            model before gapfilling it corrupts the gapfill (KeyError rxn05119_c0).
        check_autotrophy: measure whether the gapfilled model CAN fix carbon, and
            report it.  This does not raise: with CO2 uptake correctly blocked (see
            get_biolog_media) the capability never fires and the predictions are
            sound, so refusing to proceed would be wrong.  But it is a genuine
            model-quality defect and must not pass silently -- 4 of 8 Watershed
            models can build biomass from CO2 + minerals alone.  Pass
            ``strict=True`` to turn it into an error.

    Returns:
        (MSModelUtil, info) where info = {n_reactions_added, autotrophic,
        co2_biomass}.

    Note:
        The gapfiller does NOT exploit CO2 -- ``Carbon-Pyruvic-Acid`` already
        declares cpd00011 with maxFlux=0, so CO2 is not importable while gapfilling.
        The fixation capability is EMERGENT: reactions added to build biomass from
        pyruvate incidentally close a carbon-fixing loop, which only becomes
        exploitable later because build_media() hands the model CO2 the medium
        forbids.  Blocking CO2 at simulation time is therefore the correct fix for
        the predictions; the latent route remains a model-quality issue.
    """
    from modelseedpy import MSGapfill

    fba = fba or get_msfba()
    cobra_model = model.model if hasattr(model, "model") else model
    cobra_model = cobra_model.copy()
    cobra_model.objective = biomass

    # Constrain the leaky-reversible reactions BEFORE gapfilling, so the gapfill
    # solves against a model that cannot conjure reducing power. See
    # LEAKY_REVERSIBLE_RXNS.
    n_fixed = fix_leaky_directionality(cobra_model)

    before = {r.id for r in cobra_model.reactions}
    # The gapfill media and template are identical for every model, and each is a
    # KBase fetch. Accept pre-fetched objects so a sweep pulls them ONCE per worker
    # rather than once per model -- ~2 x N redundant workspace calls otherwise, which
    # is enough to trip KBase's rate limit on a full 519-genome run.
    if media_obj is None:
        media_obj = fba.get_media(gapfill_media, KBASEMEDIA_WS)
    if template_obj is None:
        template_obj = fba.get_template(template, template_ws)

    gf = MSGapfill(cobra_model, [template_obj], [], {}, {}, [])
    solution = gf.run_gapfilling(media_obj, biomass)
    if not solution:
        raise RuntimeError(f"no gapfill solution on {gapfill_media}")

    # cumulative_solution MUST be explicit -- modelseedpy's default is a shared
    # mutable list that leaks reactions between models (see defect 2 above).
    gf.integrate_gapfill_solution(solution, cumulative_solution=[])
    mdlutl = gf.mdlutl  # integrate mutates in place and returns a dict, not a model

    n_added = len([
        r for r in mdlutl.model.reactions
        if r.id not in before and not r.id.startswith(("EX_", "DM_", "SK_", "bio"))
    ])
    info = {
        "n_reactions_added": n_added,
        "autotrophic": False, "co2_biomass": 0.0,
        "free_lunch": False, "mineral_biomass": 0.0,
        "directionality_fixed": n_fixed,
    }
    if check_autotrophy:
        # TWO checks, and they are NOT the same failure.
        #
        #   free lunch -- biomass from minerals ALONE, no carbon of any kind.
        #                 Thermodynamically impossible; predictions are worthless;
        #                 NO media-level fix exists (there is nothing to block).
        #                 -> FATAL. The model must not be scored.
        #
        #   autotrophy -- biomass from CO2 + minerals.
        #                 Curable: block CO2 uptake (which the media already specify
        #                 via maxFlux=0) and the route never fires.
        #                 -> warn; predictions remain sound.
        #
        # Testing them together (minerals + CO2 in one shot) cannot tell them apart,
        # and reports an incurable model as mitigated. Test separately. Always.
        mineral_biomass = free_lunch_growth(mdlutl, biomass=biomass)
        info["mineral_biomass"] = mineral_biomass
        info["free_lunch"] = mineral_biomass > 1e-6

        co2_biomass = autotrophic_growth(mdlutl, biomass=biomass)
        info["co2_biomass"] = co2_biomass
        info["autotrophic"] = co2_biomass > 1e-6

        if info["free_lunch"]:
            raise FreeLunchModelError(
                f"gapfilled model grows at {mineral_biomass:.4f} on MINERALS ALONE -- no "
                "carbon source of any kind. It manufactures biomass out of salts, will "
                "'grow' on every medium regardless of content, and cannot be fixed by "
                "blocking anything. Its predictions are meaningless. Exclude it."
            )
        if info["autotrophic"]:
            message = (
                f"gapfilled model grows at {co2_biomass:.4f} on CO2 + minerals -- an "
                "autotrophic E. coli. Predictions stay sound while CO2 uptake is blocked "
                "(get_biolog_media does this by default), but the model carries a "
                "carbon-fixation route it should not have."
            )
            if strict:
                raise AutotrophicModelError(message)
            _warnings.warn(message, RuntimeWarning, stacklevel=2)
    return mdlutl, info


#: Charge-imbalanced reaction in the ModelSEED template that makes E. coli models
#: autotrophic. THE root cause of the CO2-fixation artifact in this study.
#:
#:     rxn05759  hydrogen:ferredoxin oxidoreductase
#:               2 H+ + 2 Reducedferredoxin <=> 2 Oxidizedferredoxin + H2
#:               mass balance: {'charge': 2}      <- charge off by 2; MASS balances
#:               bounds: (-1000, 1000)            <- reversible
#:
#: ModelSEED's ferredoxin pair differs by 2 charges (cpd11620 Reducedferredoxin
#: Fe2R4S6 charge=4; cpd11621 Oxidizedferredoxin Fe2R4S6 charge=6), i.e. 2 electrons
#: per ferredoxin. Oxidising 2 Fd_red therefore releases 4 electrons, but the written
#: stoichiometry makes only one H2 (2 electrons) -- 2 electrons vanish. Run in REVERSE
#: (observed flux -101.7, at the bound) it CREATES 2 free electrons. The balanced
#: pyruvate:ferredoxin oxidoreductase next to it (rxn05938) then fixes CO2 using that
#: free reducing power -> an autotrophic E. coli.
#:
#: It slips past BOTH standard guardrails:
#:   * ATP correction  -- the leak is electrons/reducing power, not ATP.
#:   * mass balance    -- Fe2R4S6 balances; only the CHARGE is wrong.
#:
#: Measured on as-delivered KBase models (562.55367/55368/55380), identical across
#: genomes because it is a TEMPLATE defect, not genome-specific:
#:     as delivered          CO2+minerals 0.3974   glucose 4.8532
#:     rxn05759 forward-only CO2+minerals 0.0000   glucose 3.3043
#: Forward-only removes the autotrophy completely and keeps real growth. Only the
#: reverse direction is the leak (H2 evolution is the physiological direction), so
#: constraining rather than deleting is the correct, minimal fix.
#:
#: NOTE the glucose drop (4.85 -> 3.30): the leak was inflating REAL growth too, so
#: this changes predictions generally -- it is not only a CO2 edge case.
FERREDOXIN_HYDROGENASE_RXN: str = "rxn05759"

#: Reactions declared REVERSIBLE that must not run both ways, keyed to the direction
#: to BLOCK. Both are the same class of bug: a reaction that is one-way in life is
#: marked reversible, and the solver runs it backwards to manufacture reducing power,
#: which then drives CO2 fixation.
#:
#:   rxn05759  hydrogen:ferredoxin oxidoreductase   BLOCK REVERSE
#:       2 H+ + 2 Reducedferredoxin <=> 2 Oxidizedferredoxin + H2   [charge off by 2]
#:       ModelSEED's ferredoxin pair differs by 2 charges, so oxidising 2 Fd_red frees
#:       4 electrons while the written product (one H2) carries only 2. Run in reverse
#:       it CREATES 2 free electrons. Forward (H2 evolution) is physiological.
#:
#:   rxn03978  assimilatory nitrite reductase        BLOCK FORWARD
#:       2 H2O + NH3 + 6 Fd_ox <=> 8 H+ + NO2- + 6 Fd_red
#:       Physiologically this runs NO2- -> NH3, CONSUMING reductant, for nitrogen
#:       assimilation. Forward it becomes NH3 -> NO2- + reduced ferredoxin, i.e.
#:       ammonia oxidation yielding free electrons. E. coli is not a nitrifier.
#:       (rxn00568/rxn00569, the NADH/NADPH nitrite reductases, are already
#:       assimilatory-only in the template and need no fix.)
#:
#: MEASURED, as-delivered KBase 'new' models (562.55367/55380/55368), growth on
#: CO2 + minerals vs glucose:
#:     as delivered                      0.3974   /  4.8532
#:     rxn03978 forward blocked          0.3076   /  4.8532   (reduced, not eliminated)
#:     NH3 removed from the medium       0.0000   /     n.a.  (NH3 is the donor)
#:
#: NOTE -- these fixes REDUCE but do not fully eliminate the autotrophy: at least one
#: further NH3-oxidation route exists that we did not chase down. That is acceptable
#: here because the phenotype scoring runs with CO2 uptake BLOCKED (get_biolog_media),
#: and the Biolog media supply no CO2 -- so the artifact cannot fire in the scored
#: predictions. The residual is a model-quality issue to raise with the ModelSEED
#: team, not a threat to these results.
LEAKY_REVERSIBLE_RXNS: dict = {
    "rxn05759": "reverse",   # blocking reverse => lower_bound = 0
    "rxn03978": "forward",   # blocking forward => upper_bound = 0
}


def fix_leaky_directionality(model, rxn_directions=None):
    """Constrain reactions whose reversibility lets the model manufacture electrons.

    See LEAKY_REVERSIBLE_RXNS for the mechanism and the measured effect. Matches a
    reaction in any compartment index (rxn05759_c0, _c1, ...). Returns the number of
    reactions constrained. Safe to call on any model.
    """
    rxn_directions = rxn_directions or LEAKY_REVERSIBLE_RXNS
    cobra_model = model.model if hasattr(model, "model") else model
    fixed = 0
    for rxn in cobra_model.reactions:
        block = rxn_directions.get(rxn.id.split("_")[0])
        if block == "reverse" and rxn.lower_bound < 0:
            rxn.lower_bound = 0.0
            fixed += 1
        elif block == "forward" and rxn.upper_bound > 0:
            rxn.upper_bound = 0.0
            fixed += 1
    return fixed


#: Backwards-compatible alias.
fix_ferredoxin_hydrogenase = fix_leaky_directionality


class AutotrophicModelError(RuntimeError):
    """A heterotroph that grows without an organic carbon source."""


class FreeLunchModelError(RuntimeError):
    """A model that builds biomass with NO carbon source at all -- not even CO2."""


#: The mineral background: water, O2, N, P, S and trace metals. NO carbon of any kind.
_MINERALS: tuple[str, ...] = (
    "cpd00001", "cpd00007", "cpd00009", "cpd00013", "cpd00048", "cpd00030",
    "cpd00034", "cpd00058", "cpd00063", "cpd00067", "cpd00099", "cpd00149",
    "cpd00205", "cpd00254", "cpd00971", "cpd10515", "cpd10516",
)


def _growth_on(model, compounds, biomass="bio1"):
    """Biomass flux with ONLY ``compounds`` importable."""
    cobra_model = model.model if hasattr(model, "model") else model
    cobra_model.objective = biomass
    exchanges = {r.id for r in cobra_model.reactions if r.id.startswith("EX_")}
    with cobra_model:
        for rid in exchanges:
            cobra_model.reactions.get_by_id(rid).lower_bound = 0.0
        for cid in compounds:
            rid = f"EX_{cid}_e0"
            if rid in exchanges:
                cobra_model.reactions.get_by_id(rid).lower_bound = -100.0
        flux = cobra_model.slim_optimize()
    return 0.0 if flux is None or flux != flux else float(flux)


def free_lunch_growth(model, biomass="bio1"):
    """Biomass on MINERALS ALONE -- no carbon source, not even CO2. Must be 0.

    A model that grows here is manufacturing biomass out of salts. It is
    thermodynamically impossible and its predictions are worthless: it will "grow"
    on all 324 media regardless of what they contain, and NO media-level fix can
    help -- there is nothing to block.

    This is a STRICTLY STRONGER failure than autotrophy, and it must be tested
    SEPARATELY.  An earlier version of this guard opened minerals AND CO2 together,
    which cannot distinguish "fixes CO2" (curable by blocking CO2) from "needs no
    carbon whatsoever" (incurable) -- and so reported the latter as mitigated when
    it was not.  562.55368 rebuilt grows at 4.81 on minerals alone and produced 113
    false positives with CO2 already blocked.

    MSBuilder rebuilds fail this materially more often than the 2022/2025 KBase
    builds do.  Check every model; exclude the ones that fail.
    """
    return _growth_on(model, _MINERALS, biomass=biomass)


def autotrophic_growth(model, biomass="bio1"):
    """Biomass on CO2 + minerals -- should be 0 for any E. coli.

    E. coli cannot fix carbon.  A gapfilled model that can is defective, and the
    defect is invisible in ordinary phenotype scoring: real growth still works (no
    false negatives), so the model merely looks like it over-predicts.

    Unlike ``free_lunch_growth`` this one IS curable: block CO2 uptake (which the
    media already specify via maxFlux=0) and the route never fires.  Measure it,
    report it, but it does not invalidate the model's predictions.
    """
    return _growth_on(model, list(_MINERALS) + [CO2_COMPOUND], biomass=biomass)


def simulate_biolog_panel(model, media, fba=None, biomass="bio1",
                          add_missing_exchanges=True):
    """Simulate the Biolog panel and return one row per medium.

    ``add_missing_exchanges=True`` mirrors the 2022 run's third positional argument
    to ``simulate_phenotypes(model, "bio1", True)``: it fabricates uptake routes for
    media compounds the model cannot otherwise import.

    Returns:
        [{media_name, panel, flux, growth}, ...]
    """
    fba = fba or get_msfba()
    result = fba.simulate_growth_phenotypes(
        model, media, add_missing_exchanges=add_missing_exchanges, biomass=biomass
    )
    rows = []
    for detail in result["details"]:
        flux = detail["simulated_flux"]
        flux = 0.0 if flux is None or flux != flux else float(flux)
        rows.append({
            "media_name": detail["media_id"],
            "panel": media_panel(detail["media_id"]),
            "flux": flux,
            "growth": flux > GROWTH_THRESHOLD,
        })
    return rows


def media_panel(media_name):
    """Biolog panel a medium belongs to: ``Carbon-D-Glucose`` -> ``Carbon``."""
    return str(media_name).split("-", 1)[0]


def get_experimental_calls(fba=None):
    """Experimental Biolog growth calls carried by the ``ecoli_biolog`` phenotype set.

    {media_name: bool} for all 324 media -- real measured growth for E. coli K-12.

    This is the referee for MODEL QUALITY: agreement with the 2022 run only tells us
    we reproduced a pipeline that we now know had a CO2 leak, whereas this tells us
    whether a prediction is actually right.

    Two boundaries on its use:

      * It is K-12, not the 32 watershed strains.  Chris's judgement is that the
        per-strain deviation is small, which makes it a sound referee for aggregate
        model quality ("does this model set over-predict growth?").
      * It is NOT valid for the diversity question.  The collaborators' paper turns
        on the 32 strains differing FROM EACH OTHER; scoring every strain against a
        single K-12 reference would erase exactly the variation being measured.
        That still needs their 32-strain Biolog data, which we do not have.

    It is also the check that caught a wrong conclusion here: a claim that a third of
    the 2022 matrix was substrate-independent artifact collapsed once the calls were
    scored against experiment -- the repeated values turned out to be the easy sugars
    (glucose/fructose/maltose/trehalose), which legitimately share a biomass yield
    because they share a pathway, and are 80% experimentally correct.  Score against
    experiment before believing any story about a signature in the numbers.
    """
    fba = fba or get_msfba()
    phenoset = fba.get_phenotypeset(BIOLOG_PHENOTYPESET, BIOLOG_PHENOTYPESET_WS)
    return {
        p.id: (float(p.experimental_value) > OBSERVED_THRESHOLD)
        for p in phenoset.phenotypes
        if p.experimental_value is not None
    }


def load_reference_calls(path=None):
    """Clayton's 2022 growth/no-growth calls (324 media x 519 models).

    This is the scoring target -- his own binarised output, not our binarisation of
    the mastersheet (though we verified the two agree on all 168,156 cells).
    """
    path = path or (DATA_DIR / "BiologMacTestOutput_On-Off.csv")
    return pd.read_csv(path, index_col=0)


def score_calls(sim_rows, reference, growth_threshold=GROWTH_THRESHOLD):
    """Score simulated growth CALLS against a reference column.

    The replication target is the calls, not the biomass magnitudes: the 2022
    absolute fluxes depend on uptake bounds that cannot be recovered from the saved
    artifacts, while the calls are what the study uses.

    Args:
        sim_rows: rows from ``simulate_biolog_panel``.
        reference: a column of ``load_reference_calls()`` (0/1), OR a
            {media_name: value} mapping scored against ``OBSERVED_THRESHOLD``.

    Returns:
        {CP, CN, FP, FN, n, accuracy, disagreements} -- disagreements lists every
        medium where the call differs, so they can be read rather than merely counted.
    """
    def _ref_call(name):
        value = reference[name]
        return bool(value) if value in (0, 1, True, False) else float(value) > OBSERVED_THRESHOLD

    CP = CN = FP = FN = 0
    disagreements = []
    for row in sim_rows:
        name = row["media_name"]
        if name not in reference:
            continue
        sim = row["flux"] > growth_threshold
        obs = _ref_call(name)
        if sim and obs:
            CP += 1
        elif (not sim) and (not obs):
            CN += 1
        else:
            if sim and not obs:
                FP += 1
            else:
                FN += 1
            disagreements.append({
                "media_name": name,
                "panel": row["panel"],
                "call": "FP" if sim else "FN",
                "simulated_flux": row["flux"],
            })
    n = CP + CN + FP + FN
    return {
        "CP": CP, "CN": CN, "FP": FP, "FN": FN, "n": n,
        "accuracy": (CP + CN) / n if n else None,
        "disagreements": disagreements,
    }


# ===========================================================================
# The three model arms
#
# The whole point of the study is to attribute the old/new divergence to a cause,
# which needs a THIRD arm.  Comparing only old vs new confounds two changes that
# happened at once -- the genomes were re-annotated AND the reconstruction code
# moved on.  Rebuilding from the ORIGINAL PATRIC genomes with CURRENT code holds
# the annotations fixed and varies only the build:
#
#     OLD      "<gid>.fbamodel"                  2022 build  x  PATRIC annotations
#     NEW      "genomeset__<gid>.contigs..."     2025 build  x  KBase RAST re-annotation
#     REBUILT  MSBuilder on the PATRIC genome    current     x  PATRIC annotations
#
#     old  vs rebuilt -> same annotations, different build -> BUILD-CODE effect
#     rebuilt vs new  -> same build, different annotations -> ANNOTATION effect
#
# Every arm then gets the IDENTICAL downstream treatment (gapfill_model ->
# simulate_biolog_panel), or the comparison means nothing.
# ===========================================================================

def get_patric_genome_ref(genome_id, fba=None, narrative_ws=NARRATIVE_WS):
    """Workspace ref of the ORIGINAL PATRIC-imported genome for ``genome_id``.

    These are the plain ``<gid>`` Genome objects in the narrative (saved
    2022-06-14) -- imported from PATRIC with the RAST annotations PATRIC already
    had, never re-annotated in KBase.  The ``genomeset__<gid>.contigs`` objects
    are the 2025 KBase bulk-RAST re-annotation and are NOT what this arm wants.
    """
    fba = fba or get_msfba()
    objs = fba.list_ws_objects(narrative_ws, type="KBaseGenomes.Genome")
    info = objs.get(str(genome_id))
    if info is None:
        raise KeyError(f"no PATRIC genome object named {genome_id!r} in ws {narrative_ws}")
    return f"{info[6]}/{info[0]}/{info[4]}"


def rebuild_from_patric(genome_id, fba=None, template=GAPFILL_TEMPLATE,
                        template_ws=GAPFILL_TEMPLATE_WS, narrative_ws=NARRATIVE_WS,
                        template_obj=None):
    """ARM 3 -- rebuild a model from the original PATRIC genome with current code.

    Gapfilling is deliberately NOT done here: every arm is gapfilled downstream by
    ``gapfill_model`` so all three receive identical treatment.  Building with
    ``gapfill_model=True`` would give this arm a different gapfill from the others
    and quietly invalidate the comparison.

    ``annotate_with_rast=False`` keeps the PATRIC annotations as they are -- the
    entire point of the arm.  Re-annotating here would reproduce the very drift we
    are trying to isolate.
    """
    from modelseedpy import MSBuilder

    fba = fba or get_msfba()
    genome = fba.get_msgenome(get_patric_genome_ref(genome_id, fba=fba,
                                                    narrative_ws=narrative_ws))
    if template_obj is None:  # accept a pre-fetched template -- see gapfill_model
        template_obj = fba.get_template(template, template_ws)
    return MSBuilder(genome, template=template_obj).build(
        f"{genome_id}.rebuilt", index="0",
        allow_all_non_grp_reactions=False, annotate_with_rast=False,
    )


def replicate_2022(genome_id, model_ref, media=None, fba=None, reference=None):
    """The whole pipeline for one model: gapfill -> simulate -> score.

    Returns:
        {genome_id, n_reactions_added, rows, score}
    """
    fba = fba or get_msfba()
    media = media if media is not None else get_biolog_media(fba=fba)
    mdlutl, info = gapfill_model(fba.get_model(model_ref), fba=fba)
    rows = simulate_biolog_panel(mdlutl, media, fba=fba)
    out = {"genome_id": genome_id, "rows": rows, **info}
    if reference is not None:
        column = f"{genome_id}.fbamodel"
        if column in reference.columns:
            out["score"] = score_calls(rows, reference[column])
    return out


def _media_ref_to_name(media_ref, fba, cache):
    """Resolve a media_ref ('wsid/objid/version') to its object name, cached."""
    if media_ref is None:
        return None
    if media_ref in cache:
        return cache[media_ref]
    try:
        info = fba.get_object_info(media_ref)
        name = info[1] if info else None
    except Exception:
        name = None
    cache[media_ref] = name
    return name


def fetch_observed_phenotypes(
    ref_genomes, fba=None, narrative_ws=NARRATIVE_WS, observed_threshold=1e-6
):
    """Fetch + normalize the narrative's observed KBasePhenotypes.PhenotypeSet objects.

    Each PhenotypeSet.phenotypes entry (see KBasePhenotypes.spec) carries a
    ``media_ref`` and a ``normalizedGrowth`` float; media identity is taken
    from the phenotype's own ``name`` field when present (collaborators
    typically stash the media name there for readability), falling back to
    resolving ``media_ref`` -> object name via ``get_object_info`` (cached
    per call).

    Args:
        ref_genomes: iterable of reference genome ids to keep.
        fba: optional MSFBAUtils facade (defaults to ``get_msfba()``).
        narrative_ws: workspace holding the observed PhenotypeSet objects.
        observed_threshold: growth-call binarization cutoff (inclusive).

    Returns:
        DataFrame with columns [genome_id, media_name, growth_call] -- one
        row per (genome, media) observed phenotype.
    """
    fba = fba or get_msfba()
    ref_set = set(str(g) for g in ref_genomes)
    objs = fba.list_ws_objects(narrative_ws, type="KBasePhenotypes.PhenotypeSet")

    media_name_cache = {}
    rows = []
    for name, info in objs.items():
        meta = info[10] if len(info) > 10 else {}
        gid, _reason = _normalize_genome_id(str(name), meta, ref_set)
        if gid is None:
            continue
        objid, version, obj_wsid = info[0], info[4], info[6]
        ref = f"{obj_wsid}/{objid}/{version}"
        data = fba.get_object(ref) or {}
        for pheno in data.get("phenotypes", []):
            media_name = pheno.get("name") or _media_ref_to_name(
                pheno.get("media_ref"), fba, media_name_cache
            )
            if media_name is None:
                continue
            growth = pheno.get("normalizedGrowth")
            rows.append({
                "genome_id": gid,
                "media_name": media_name,
                "growth_call": bool(growth is not None and float(growth) >= observed_threshold),
            })

    return pd.DataFrame(rows, columns=["genome_id", "media_name", "growth_call"])


def fetch_narrative_simulated_growth(
    ref_genomes, fba=None, narrative_ws=NARRATIVE_WS, name_suffix="_SimulatedGrowth"
):
    """Fetch + normalize the narrative's "..._SimulatedGrowth" simulation sets.

    Lists KBasePhenotypes.PhenotypeSimulationSet objects in ``narrative_ws``
    whose name ends with ``name_suffix`` (the collaborators' own prior
    simulation export -- this is the reference ``integrity_check`` compares
    the mastersheet against). Each PhenotypeSimulationSet.phenotypeSimulations
    entry carries a ``simulatedGrowth`` float and a ``phenotype_ref`` pointing
    back into the linked ``phenotypeset_ref`` (a KBasePhenotypes.PhenotypeSet);
    media identity is recovered by fetching that PhenotypeSet once (cached per
    call) and matching phenotype ids.

    Args:
        ref_genomes: iterable of reference genome ids to keep.
        fba: optional MSFBAUtils facade (defaults to ``get_msfba()``).
        narrative_ws: workspace holding the simulation-set objects.
        name_suffix: only objects whose name ends with this suffix are
            scanned (set to "" / None to scan all PhenotypeSimulationSet
            objects in the workspace).

    Returns:
        DataFrame with columns [genome_id, media_name, value, version] --
        one row per (genome, media, saved-object-version) triple. Multiple
        ``version`` rows can exist for the same (genome_id, media_name) if
        the object was re-saved; ``integrity_check`` keeps only the latest.
    """
    fba = fba or get_msfba()
    ref_set = set(str(g) for g in ref_genomes)
    objs = fba.list_ws_objects(narrative_ws, type="KBasePhenotypes.PhenotypeSimulationSet")

    phenoset_cache = {}
    rows = []
    for name, info in objs.items():
        name = str(name)
        if name_suffix and not name.endswith(name_suffix):
            continue
        meta = info[10] if len(info) > 10 else {}
        gid = _extract_genome_id(name, ref_set)
        if gid is None:
            gid, _reason = _normalize_genome_id(name, meta, ref_set)
        if gid is None or gid not in ref_set:
            continue

        objid, version, obj_wsid = info[0], info[4], info[6]
        ref = f"{obj_wsid}/{objid}/{version}"
        data = fba.get_object(ref) or {}

        phenoset_ref = data.get("phenotypeset_ref")
        media_by_phenotype_id = phenoset_cache.get(phenoset_ref)
        if media_by_phenotype_id is None:
            media_by_phenotype_id = {}
            if phenoset_ref:
                pdata = fba.get_object(phenoset_ref) or {}
                for pheno in pdata.get("phenotypes", []):
                    pid = pheno.get("id")
                    if pid is not None:
                        media_by_phenotype_id[pid] = pheno.get("name") or _media_ref_to_name(
                            pheno.get("media_ref"), fba, {}
                        )
            phenoset_cache[phenoset_ref] = media_by_phenotype_id

        for sim in data.get("phenotypeSimulations", []):
            phenotype_ref = sim.get("phenotype_ref")
            pid = str(phenotype_ref).rsplit("/", 1)[-1] if phenotype_ref else None
            media_name = media_by_phenotype_id.get(pid) or media_by_phenotype_id.get(phenotype_ref)
            if media_name is None:
                continue
            rows.append({
                "genome_id": gid,
                "media_name": media_name,
                "value": sim.get("simulatedGrowth"),
                "version": version,
            })

    return pd.DataFrame(rows, columns=["genome_id", "media_name", "value", "version"])
