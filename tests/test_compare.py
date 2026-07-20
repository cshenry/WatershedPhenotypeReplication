"""Tests for util.score_against_reference / util.integrity_check and the
KBase-facing normalization helpers (resolve_media, fetch_observed_phenotypes,
fetch_narrative_simulated_growth).

score_against_reference / integrity_check are pure -- no KBase calls, no
fixtures beyond plain DataFrames. The fetch/resolve helpers DO shell out to
an MSFBAUtils facade, so they are exercised here against an in-memory
FakeFBA stub (list_ws_objects / get_object / get_object_info), matching the
convention in test_inventory.py / test_models_partition.py -- no live KBase
session is ever touched.
"""
from __future__ import annotations

import math

import pandas as pd
import pytest

import util

GROWTH_THRESHOLD = 0.01
OBSERVED_THRESHOLD = 1e-6


# ---------------------------------------------------------------------------
# score_against_reference
# ---------------------------------------------------------------------------


def _long_df(rows):
    return pd.DataFrame(rows, columns=["genome_id", "media_name", "value"])


def test_score_against_reference_counts_cp_cn_fp_fn():
    # g1: media A (CP: both grow), media B (CN: both no-grow)
    # g2: media A (FP: pred grows, ref doesn't), media B (FN: pred doesn't, ref does)
    pred = _long_df([
        {"genome_id": "g1", "media_name": "A", "value": 0.5},
        {"genome_id": "g1", "media_name": "B", "value": 0.0},
        {"genome_id": "g2", "media_name": "A", "value": 0.5},
        {"genome_id": "g2", "media_name": "B", "value": 0.0},
    ])
    ref = _long_df([
        {"genome_id": "g1", "media_name": "A", "value": 0.3},
        {"genome_id": "g1", "media_name": "B", "value": 0.0},
        {"genome_id": "g2", "media_name": "A", "value": 0.0},
        {"genome_id": "g2", "media_name": "B", "value": 0.4},
    ])

    result = util.score_against_reference(
        pred, ref, growth_threshold=GROWTH_THRESHOLD, observed_threshold=OBSERVED_THRESHOLD
    )
    agg = result["aggregate"]

    assert agg["CP"] == 1
    assert agg["CN"] == 1
    assert agg["FP"] == 1
    assert agg["FN"] == 1
    assert agg["n"] == 4
    assert agg["accuracy"] == pytest.approx(0.5)

    per_model = result["per_model"]
    assert per_model.loc["g1", "CP"] == 1
    assert per_model.loc["g1", "CN"] == 1
    assert per_model.loc["g1", "accuracy"] == pytest.approx(1.0)
    assert per_model.loc["g2", "FP"] == 1
    assert per_model.loc["g2", "FN"] == 1
    assert per_model.loc["g2", "accuracy"] == pytest.approx(0.0)


def test_score_against_reference_value_exactly_at_threshold_is_growth():
    # pred value exactly equal to growth_threshold must binarize to growth (>=)
    pred = _long_df([{"genome_id": "g1", "media_name": "A", "value": GROWTH_THRESHOLD}])
    ref = _long_df([{"genome_id": "g1", "media_name": "A", "value": 1.0}])  # ref grows

    result = util.score_against_reference(
        pred, ref, growth_threshold=GROWTH_THRESHOLD, observed_threshold=OBSERVED_THRESHOLD
    )
    merged = result["merged"]
    row = merged.iloc[0]
    assert row["pred_growth"] == True  # noqa: E712 -- explicit bool check
    assert row["class"] == "CP"
    assert result["aggregate"]["CP"] == 1

    # symmetric check on the reference side: ref value exactly at
    # observed_threshold must also binarize to growth (>=).
    ref_edge = _long_df([{"genome_id": "g1", "media_name": "A", "value": OBSERVED_THRESHOLD}])
    pred_grows = _long_df([{"genome_id": "g1", "media_name": "A", "value": 1.0}])
    result2 = util.score_against_reference(
        pred_grows, ref_edge, growth_threshold=GROWTH_THRESHOLD, observed_threshold=OBSERVED_THRESHOLD
    )
    assert result2["merged"].iloc[0]["ref_growth"] == True  # noqa: E712
    assert result2["aggregate"]["CP"] == 1


def test_score_against_reference_nan_prediction_scores_as_no_growth():
    # NaN prediction (e.g. a model_load_failed sweep row) with a known
    # reference that DOES grow -> scored as FN (no-growth prediction vs a
    # growing reference), not excluded.
    pred = _long_df([{"genome_id": "g1", "media_name": "A", "value": float("nan")}])
    ref = _long_df([{"genome_id": "g1", "media_name": "A", "value": 1.0}])

    result = util.score_against_reference(
        pred, ref, growth_threshold=GROWTH_THRESHOLD, observed_threshold=OBSERVED_THRESHOLD
    )
    agg = result["aggregate"]
    assert agg["n"] == 1
    assert agg["FN"] == 1
    assert agg["CP"] == 0


def test_score_against_reference_both_nan_excluded():
    # Both prediction and reference NaN/missing for a (genome, media) pair
    # -> excluded entirely, not counted as a correct negative.
    pred = _long_df([{"genome_id": "g1", "media_name": "A", "value": float("nan")}])
    ref = _long_df([{"genome_id": "g1", "media_name": "A", "value": float("nan")}])

    result = util.score_against_reference(
        pred, ref, growth_threshold=GROWTH_THRESHOLD, observed_threshold=OBSERVED_THRESHOLD
    )
    agg = result["aggregate"]
    assert agg["n"] == 0
    assert agg["CP"] == agg["CN"] == agg["FP"] == agg["FN"] == 0
    assert math.isnan(agg["accuracy"])
    assert result["merged"].iloc[0]["class"] is None


def test_score_against_reference_missing_key_on_one_side_treated_as_nan():
    # A (genome, media) pair present in pred but entirely absent from ref
    # (no row at all) behaves the same as an explicit NaN reference value --
    # the pair is included with ref_growth=False (no-growth reference).
    pred = _long_df([
        {"genome_id": "g1", "media_name": "A", "value": 0.5},
        {"genome_id": "g1", "media_name": "B", "value": 0.5},
    ])
    ref = _long_df([{"genome_id": "g1", "media_name": "A", "value": 1.0}])

    result = util.score_against_reference(
        pred, ref, growth_threshold=GROWTH_THRESHOLD, observed_threshold=OBSERVED_THRESHOLD
    )
    agg = result["aggregate"]
    # A: CP (both grow); B: ref missing entirely -> ref_growth False, pred
    # grows -> FP (not excluded, since pred is NOT NaN).
    assert agg["CP"] == 1
    assert agg["FP"] == 1
    assert agg["n"] == 2


def test_score_against_reference_correlation_drops_nan_pairs():
    pred = _long_df([
        {"genome_id": "g1", "media_name": "A", "value": 0.1},
        {"genome_id": "g1", "media_name": "B", "value": 0.2},
        {"genome_id": "g1", "media_name": "C", "value": 0.3},
        {"genome_id": "g1", "media_name": "D", "value": float("nan")},  # dropped from corr
    ])
    ref = _long_df([
        {"genome_id": "g1", "media_name": "A", "value": 1.0},
        {"genome_id": "g1", "media_name": "B", "value": 2.0},
        {"genome_id": "g1", "media_name": "C", "value": 3.0},
        {"genome_id": "g1", "media_name": "D", "value": 4.0},  # pred NaN -> pair dropped
    ])

    result = util.score_against_reference(
        pred, ref, growth_threshold=GROWTH_THRESHOLD, observed_threshold=OBSERVED_THRESHOLD
    )
    agg = result["aggregate"]
    # Perfectly linear over the 3 complete pairs -> pearson/spearman == 1.0
    assert agg["pearson_r"] == pytest.approx(1.0)
    assert agg["spearman_r"] == pytest.approx(1.0)


def test_score_against_reference_custom_column_names():
    pred = pd.DataFrame([
        {"genome_id": "g1", "media_name": "A", "simulated_flux": 0.5},
    ])
    ref = pd.DataFrame([
        {"genome_id": "g1", "media_name": "A", "mastersheet_value": 1.0},
    ])
    result = util.score_against_reference(
        pred, ref, pred_value_col="simulated_flux", ref_value_col="mastersheet_value"
    )
    assert result["aggregate"]["CP"] == 1


# ---------------------------------------------------------------------------
# integrity_check
# ---------------------------------------------------------------------------


def test_integrity_check_flags_cells_above_tolerance():
    mastersheet_df = pd.DataFrame(
        {"g1": [0.5, 0.0], "g2": [0.3, 0.1]},
        index=["media_A", "media_B"],
    )
    simset_df = pd.DataFrame([
        {"genome_id": "g1", "media_name": "media_A", "value": 0.5, "version": 1},
        {"genome_id": "g1", "media_name": "media_B", "value": 0.0000001, "version": 1},  # within tolerance
        {"genome_id": "g2", "media_name": "media_A", "value": 0.9, "version": 1},  # flagged
        {"genome_id": "g2", "media_name": "media_B", "value": 0.1, "version": 1},
    ])

    compared, flagged = util.integrity_check(mastersheet_df, simset_df, tolerance=1e-6)

    assert len(compared) == 4
    assert len(flagged) == 1
    assert flagged.iloc[0]["genome_id"] == "g2"
    assert flagged.iloc[0]["media_name"] == "media_A"
    assert flagged.iloc[0]["abs_diff"] == pytest.approx(0.6)


def test_integrity_check_keeps_only_latest_version():
    mastersheet_df = pd.DataFrame({"g1": [0.5]}, index=["media_A"])
    simset_df = pd.DataFrame([
        # stale version disagrees; latest version matches -> not flagged
        {"genome_id": "g1", "media_name": "media_A", "value": 0.9, "version": 1},
        {"genome_id": "g1", "media_name": "media_A", "value": 0.5, "version": 2},
    ])

    compared, flagged = util.integrity_check(mastersheet_df, simset_df, tolerance=1e-6)

    assert len(compared) == 1
    assert compared.iloc[0]["simset_value"] == pytest.approx(0.5)
    assert len(flagged) == 0


def test_integrity_check_no_row_for_missing_pairs():
    mastersheet_df = pd.DataFrame({"g1": [0.5, 0.2]}, index=["media_A", "media_B"])
    simset_df = pd.DataFrame([
        {"genome_id": "g1", "media_name": "media_A", "value": 0.5, "version": 1},
        # media_B has no simset row at all -> not in compared (inner join)
    ])

    compared, _flagged = util.integrity_check(mastersheet_df, simset_df, tolerance=1e-6)

    assert len(compared) == 1
    assert set(compared["media_name"]) == {"media_A"}


# ---------------------------------------------------------------------------
# resolve_media / fetch_observed_phenotypes / fetch_narrative_simulated_growth
# (KBase-facing helpers, exercised offline against an in-memory FakeFBA)
# ---------------------------------------------------------------------------


class FakeFBA:
    """Stands in for MSFBAUtils: list_ws_objects / get_object / get_object_info
    answered from canned in-memory dicts, mirroring the tuple shape used by
    test_inventory.py / test_models_partition.py.
    """

    def __init__(self, by_workspace, objects_by_ref=None, info_by_ref=None):
        self._by_workspace = by_workspace
        self._objects_by_ref = objects_by_ref or {}
        self._info_by_ref = info_by_ref or {}

    def list_ws_objects(self, wsid_or_ref, type=None, include_metadata=True):
        return dict(self._by_workspace.get(wsid_or_ref, {}))

    def get_object(self, ref):
        return self._objects_by_ref.get(ref)

    def get_object_info(self, ref):
        return self._info_by_ref.get(ref)


def _obj(objid, name, ktype, save_date, version, wsid, ws_name, meta=None):
    return (objid, name, ktype, save_date, version, "chenry", wsid, ws_name, "chsum", 1234, meta or {})


def test_resolve_media_flags_unresolved():
    media_objs = {
        "Carbon-D-Glucose": _obj(1, "Carbon-D-Glucose", "KBaseBiochem.Media-1.0", "2019-01-01", 1, 99000, "KBaseMedia"),
        "Carbon-Acetate": _obj(2, "Carbon-Acetate", "KBaseBiochem.Media-1.0", "2019-01-02", 1, 99000, "KBaseMedia"),
    }
    fba = FakeFBA(
        {"KBaseMedia": media_objs},
        objects_by_ref={
            "99000/1/1": {"mediacompounds": []},
            "99000/2/1": {"mediacompounds": []},
        },
    )

    resolved, unresolved = util.resolve_media(
        ["Carbon-D-Glucose", "Carbon-Acetate", "Carbon-Not-A-Real-Media"], fba=fba
    )

    assert set(resolved) == {"Carbon-D-Glucose", "Carbon-Acetate"}
    assert unresolved == ["Carbon-Not-A-Real-Media"]
    assert resolved["Carbon-D-Glucose"] == {"mediacompounds": []}


def test_fetch_observed_phenotypes_normalizes_to_media_name_growth_call():
    ref_genomes = ["562.55367"]
    narrative_objs = {
        "562.55367.phenotypeset": _obj(
            1, "562.55367.phenotypeset", "KBasePhenotypes.PhenotypeSet-1.0",
            "2020-01-01", 1, util.NARRATIVE_WS, "hope_narrative",
        ),
    }
    phenoset_data = {
        "phenotypes": [
            {"id": "ph1", "name": "Carbon-D-Glucose", "media_ref": "99000/1/1", "normalizedGrowth": 0.8},
            {"id": "ph2", "name": "Carbon-Acetate", "media_ref": "99000/2/1", "normalizedGrowth": 0.0},
        ],
    }
    fba = FakeFBA(
        {util.NARRATIVE_WS: narrative_objs},
        objects_by_ref={f"{util.NARRATIVE_WS}/1/1": phenoset_data},
    )

    df = util.fetch_observed_phenotypes(ref_genomes, fba=fba)

    assert list(df.columns) == ["genome_id", "media_name", "growth_call"]
    assert len(df) == 2
    glucose_row = df[df["media_name"] == "Carbon-D-Glucose"].iloc[0]
    assert glucose_row["genome_id"] == "562.55367"
    assert glucose_row["growth_call"] == True  # noqa: E712
    acetate_row = df[df["media_name"] == "Carbon-Acetate"].iloc[0]
    assert acetate_row["growth_call"] == False  # noqa: E712


def test_fetch_narrative_simulated_growth_only_matches_suffix_and_resolves_media():
    ref_genomes = ["562.55367"]
    narrative_objs = {
        "562.55367_SimulatedGrowth": _obj(
            1, "562.55367_SimulatedGrowth", "KBasePhenotypes.PhenotypeSimulationSet-1.0",
            "2020-01-01", 2, util.NARRATIVE_WS, "hope_narrative",
        ),
        # does not match the "_SimulatedGrowth" suffix -> skipped
        "562.55367_other": _obj(
            2, "562.55367_other", "KBasePhenotypes.PhenotypeSimulationSet-1.0",
            "2020-01-02", 1, util.NARRATIVE_WS, "hope_narrative",
        ),
    }
    phenoset_ref = f"{util.NARRATIVE_WS}/100/1"
    simset_data = {
        "phenotypeset_ref": phenoset_ref,
        "phenotypeSimulations": [
            {"phenotype_ref": f"{phenoset_ref}/phenotypes/id/ph1", "simulatedGrowth": 0.42},
        ],
    }
    phenoset_data = {
        "phenotypes": [
            {"id": "ph1", "name": "Carbon-D-Glucose", "media_ref": "99000/1/1"},
        ],
    }
    fba = FakeFBA(
        {util.NARRATIVE_WS: narrative_objs},
        objects_by_ref={
            f"{util.NARRATIVE_WS}/1/2": simset_data,
            phenoset_ref: phenoset_data,
        },
    )

    df = util.fetch_narrative_simulated_growth(ref_genomes, fba=fba)

    assert list(df.columns) == ["genome_id", "media_name", "value", "version"]
    assert len(df) == 1
    row = df.iloc[0]
    assert row["genome_id"] == "562.55367"
    assert row["media_name"] == "Carbon-D-Glucose"
    assert row["value"] == pytest.approx(0.42)
    assert row["version"] == 2
