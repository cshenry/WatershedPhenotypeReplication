"""Tests for util.census_artifacts (see conftest.py for the import bootstrap).

All fixtures are in-memory -- no live KBase workspace calls. A FakeFBA stands
in for MSFBAUtils and answers list_ws_objects() from a canned dict keyed by
workspace id/name, using the same KBase list_objects tuple shape that
kbutillib.kb_ws_utils.list_ws_objects returns:

    (objid, name, type, save_date, version, saved_by, wsid, ws_name, chsum,
     size, meta)
"""
from __future__ import annotations

import util

NARRATIVE_WS = util.NARRATIVE_WS  # 119455
MODELS_WS = util.MODELS_WS  # 220871
MEDIA_WS = "KBaseMedia"

REF_GENOMES = ["562.55367", "562.55510", "83333.111"]


class FakeFBA:
    """Stands in for MSFBAUtils; answers list_ws_objects from canned data."""

    def __init__(self, by_workspace):
        self._by_workspace = by_workspace

    def list_ws_objects(self, wsid_or_ref, type=None, include_metadata=True):
        return dict(self._by_workspace.get(wsid_or_ref, {}))


def _obj(objid, name, ktype, save_date, version, wsid, ws_name, meta=None):
    return (
        objid,
        name,
        ktype,
        save_date,
        version,
        "chenry",
        wsid,
        ws_name,
        "chsum",
        1234,
        meta or {},
    )


def _make_fba_with_media():
    """Scenario: genome A (562.55367) has every artifact class; genome B
    (562.55510) has only old_model + new_model (via messy names / metadata
    genome_ref); genome C (83333.111) has nothing. Media is present (global).
    Also includes three unmappable objects: one with no genome token at all,
    one with an ambiguous metadata genome_ref, and one with an ambiguous
    name.
    """
    narrative_objs = {
        "562.55367": _obj(1, "562.55367", "KBaseGenomes.Genome-8.0", "2020-01-01T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative"),
        "562.55367.fbamodel": _obj(2, "562.55367.fbamodel", "KBaseFBA.FBAModel-8.0", "2020-01-02T00:00:00+0000", 3, NARRATIVE_WS, "hope_narrative"),
        "562.55510Old2022.fbamodel": _obj(3, "562.55510Old2022.fbamodel", "KBaseFBA.FBAModel-8.0", "2020-01-03T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative"),
        "562.55367.phenotypeset": _obj(4, "562.55367.phenotypeset", "KBasePhenotypes.PhenotypeSet-1.0", "2020-01-04T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative"),
        "562.55367.simset": _obj(5, "562.55367.simset", "KBasePhenotypes.PhenotypeSimulationSet-1.0", "2020-01-05T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative"),
        "562.55367.assembly": _obj(6, "562.55367.assembly", "KBaseGenomeAnnotations.Assembly-1.0", "2020-01-06T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative"),
        # unmappable: no genome token anywhere in name/meta
        "molten.fbamodel": _obj(7, "molten.fbamodel", "KBaseFBA.FBAModel-8.0", "2020-01-07T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative"),
        # unmappable: ambiguous metadata genome_ref (two distinct ref hits)
        "combo.fbamodel": _obj(
            8, "combo.fbamodel", "KBaseFBA.FBAModel-8.0", "2020-01-08T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative",
            meta={"genome_ref": "562.55367;562.55510"},
        ),
        # unmappable: ambiguous name (two distinct ref hits, no metadata)
        "cross_562.55367_562.55510.fbamodel": _obj(9, "cross_562.55367_562.55510.fbamodel", "KBaseFBA.FBAModel-8.0", "2020-01-09T00:00:00+0000", 1, NARRATIVE_WS, "hope_narrative"),
    }
    models_objs = {
        "genomeset_5_WS_562.55367.fa_assembly.fbamodel": _obj(10, "genomeset_5_WS_562.55367.fa_assembly.fbamodel", "KBaseFBA.FBAModel-8.0", "2021-01-01T00:00:00+0000", 1, MODELS_WS, "hope_models"),
        # normalized via metadata genome_ref rather than name
        "model_A": _obj(
            11, "model_A", "KBaseFBA.FBAModel-8.0", "2025-10-03T00:00:00+0000", 2, MODELS_WS, "hope_models",
            meta={"genome_ref": "562.55510"},
        ),
    }
    media_objs = {
        "Carbon-D-Glucose": _obj(1, "Carbon-D-Glucose", "KBaseBiochem.Media-1.0", "2019-01-01T00:00:00+0000", 1, 99000, "KBaseMedia"),
    }
    return FakeFBA({NARRATIVE_WS: narrative_objs, MODELS_WS: models_objs, MEDIA_WS: media_objs})


def test_census_artifacts_table_rows():
    fba = _make_fba_with_media()
    table_df, missing_report = util.census_artifacts(
        [NARRATIVE_WS, MODELS_WS, MEDIA_WS], ref_genomes=REF_GENOMES, fba=fba
    )

    assert list(table_df.columns) == [
        "genome_id",
        "artifact_class",
        "workspace",
        "object_name",
        "ref",
        "version",
        "save_date",
    ]

    # 7 mappable objects: genome, old_model x2, observed_phenotype_set,
    # prior_simulation_set, assembly, new_model x2, media -> 1+2+1+1+1+2+1 = 9
    assert len(table_df) == 9

    row = table_df[table_df["object_name"] == "562.55367.fbamodel"].iloc[0]
    assert row["genome_id"] == "562.55367"
    assert row["artifact_class"] == "old_model"
    assert row["workspace"] == NARRATIVE_WS
    assert row["ref"] == f"{NARRATIVE_WS}/2/3"
    assert row["version"] == 3

    # messy name -> old_model for 562.55510
    messy = table_df[table_df["object_name"] == "562.55510Old2022.fbamodel"].iloc[0]
    assert messy["genome_id"] == "562.55510"
    assert messy["artifact_class"] == "old_model"

    # new_model normalized via name token inside a longer messy name
    new_a = table_df[table_df["object_name"] == "genomeset_5_WS_562.55367.fa_assembly.fbamodel"].iloc[0]
    assert new_a["genome_id"] == "562.55367"
    assert new_a["artifact_class"] == "new_model"

    # new_model normalized via metadata genome_ref (name has no ref token)
    new_b = table_df[table_df["object_name"] == "model_A"].iloc[0]
    assert new_b["genome_id"] == "562.55510"
    assert new_b["artifact_class"] == "new_model"

    # media is global -- genome_id is None
    media_row = table_df[table_df["object_name"] == "Carbon-D-Glucose"].iloc[0]
    assert media_row["artifact_class"] == "media"
    assert media_row["genome_id"] is None

    # the three unmappable objects are excluded from the table
    assert "molten.fbamodel" not in set(table_df["object_name"])
    assert "combo.fbamodel" not in set(table_df["object_name"])
    assert "cross_562.55367_562.55510.fbamodel" not in set(table_df["object_name"])


def test_census_artifacts_unmappable_reasons():
    fba = _make_fba_with_media()
    _table_df, missing_report = util.census_artifacts(
        [NARRATIVE_WS, MODELS_WS, MEDIA_WS], ref_genomes=REF_GENOMES, fba=fba
    )

    unmappable_by_name = {u["object_name"]: u for u in missing_report["unmappable"]}
    assert set(unmappable_by_name) == {
        "molten.fbamodel",
        "combo.fbamodel",
        "cross_562.55367_562.55510.fbamodel",
    }
    assert unmappable_by_name["molten.fbamodel"]["reason"] == "no_match"
    assert unmappable_by_name["combo.fbamodel"]["reason"] == "ambiguous"
    assert unmappable_by_name["cross_562.55367_562.55510.fbamodel"]["reason"] == "ambiguous"


def test_census_artifacts_flags_exactly_the_absent_classes_per_genome():
    fba = _make_fba_with_media()
    _table_df, missing_report = util.census_artifacts(
        [NARRATIVE_WS, MODELS_WS, MEDIA_WS], ref_genomes=REF_GENOMES, fba=fba
    )

    per_genome = missing_report["per_genome"]
    assert set(per_genome) == set(REF_GENOMES)

    # A has every genome-scoped class -- nothing missing (media is present globally)
    assert per_genome["562.55367"] == []

    # B has old_model + new_model only
    assert set(per_genome["562.55510"]) == {
        "genome",
        "observed_phenotype_set",
        "prior_simulation_set",
        "assembly",
    }

    # C has nothing at all
    assert set(per_genome["83333.111"]) == {
        "genome",
        "old_model",
        "new_model",
        "observed_phenotype_set",
        "prior_simulation_set",
        "assembly",
    }

    assert missing_report["media_available"] is True
    assert missing_report["counts"] == {
        "genome": 1,
        "old_model": 2,
        "new_model": 2,
        "observed_phenotype_set": 1,
        "prior_simulation_set": 1,
        "assembly": 1,
        "media": 1,
    }
    assert missing_report["missing_counts"] == {
        "genome": 2,
        "old_model": 1,
        "new_model": 1,
        "observed_phenotype_set": 2,
        "prior_simulation_set": 2,
        "assembly": 2,
        "media": 0,
    }


def test_census_artifacts_no_media_flags_all_genomes_missing_media():
    fba = _make_fba_with_media()
    # Drop the media workspace entirely.
    _table_df, missing_report = util.census_artifacts(
        [NARRATIVE_WS, MODELS_WS], ref_genomes=REF_GENOMES, fba=fba
    )

    assert missing_report["media_available"] is False
    for gid in REF_GENOMES:
        assert "media" in missing_report["per_genome"][gid]
    assert missing_report["missing_counts"]["media"] == len(REF_GENOMES)
    assert missing_report["counts"]["media"] == 0
