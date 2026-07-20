"""Tests for util.partition_models (see conftest.py for the import bootstrap).

No live KBase calls -- MSFBAUtils.list_ws_objects is replaced with an in-memory
FakeFBA fixture.

WHAT THIS FILE IS DEFENDING
---------------------------
The old/new split decides which models form the 2022 REFERENCE arm of the study.
Getting it wrong silently puts 2025 models into that arm, and every downstream number
is then quietly comparing 2025 against 2025.

It has been got wrong twice, in two different ways, and both are pinned here:

  1. Splitting by WORKSPACE.  The narrative (119455) holds BOTH generations -- 519
     <gid>.fbamodel from 2022 AND 524 genomeset__* from 2025 -- so a workspace rule
     labels the 2025 rebuilds "old".

  2. Splitting by NAME alone.  Also insufficient:
       * the genome-set marker appears as genomeset__, genomeset_5_WS_,
         genomeset_WS5_10-25_ and GenomeSetATCC_ -- a "genomeset__" prefix test
         misses 13 real models;
       * 562.61170 has TWO objects with the IDENTICAL name "562.61170.fbamodel",
         one saved 2022 (v1) and one saved 2025 (v3) -- picking the highest VERSION
         selects the 2025 one;
       * 18 genomes have plain <gid>.fbamodel objects in the narrative saved
         2025-10-03 -- they look old and are not;
       * "562.55510Old2022.fbamodel" is saved 2026-06-15. The name lies.

A model is OLD only if it carries no genome-set marker AND was saved before
OLD_MODEL_CUTOFF. Where a genome still has several candidates, the old arm takes the
EARLIEST-saved object and the new arm the latest -- never "highest version".
"""
from __future__ import annotations

import util

NARRATIVE_WS = util.NARRATIVE_WS  # 119455
MODELS_WS = util.MODELS_WS  # 220871

REF_GENOMES = ["562.55367", "562.55510", "562.85951", "562.61170"]

OLD_DATE = "2022-06-15T00:00:00+0000"
NEW_DATE = "2025-10-03T00:00:00+0000"


class FakeFBA:
    """Stands in for MSFBAUtils.list_ws_objects(wsid, type=...)."""

    def __init__(self, by_workspace):
        self._by_workspace = by_workspace

    def list_ws_objects(self, wsid_or_ref, type=None, include_metadata=True):
        return dict(self._by_workspace.get(wsid_or_ref, {}))


def _obj(objid, name, save_date, version, wsid):
    return (objid, name, "KBaseFBA.FBAModel-8.0", save_date, version,
            "chenry", wsid, "ws", "chsum", 1234, {})


def _make_fba():
    narrative = {
        # a genuine 2022 model
        "562.55367.fbamodel": _obj(1, "562.55367.fbamodel", OLD_DATE, 1, NARRATIVE_WS),
        # messy name, still a 2022 build
        "562.55510Old2022.fbamodel": _obj(2, "562.55510Old2022.fbamodel", OLD_DATE, 1, NARRATIVE_WS),
        # no extractable genome token -> unmappable
        "notaref.fbamodel": _obj(3, "notaref.fbamodel", OLD_DATE, 1, NARRATIVE_WS),
        # genome-shaped token, not a reference genome -> excluded
        "999.999.fbamodel": _obj(4, "999.999.fbamodel", OLD_DATE, 1, NARRATIVE_WS),
        # the narrative ALSO holds 2025 genomeset rebuilds -- a WORKSPACE split calls
        # these "old"
        "genomeset__562.55367.contigs.fbamodel": _obj(5, "genomeset__562.55367.contigs.fbamodel", NEW_DATE, 1, NARRATIVE_WS),
        # single-underscore variant -- a "genomeset__" PREFIX test misses this
        "genomeset_WS5_10-25_562.55510.fa_assembly.fbamodel": _obj(6, "genomeset_WS5_10-25_562.55510.fa_assembly.fbamodel", NEW_DATE, 1, NARRATIVE_WS),
        # plain name, but built in 2025 -- looks old, is not
        "562.85951.fbamodel": _obj(7, "562.85951.fbamodel", NEW_DATE, 2, NARRATIVE_WS),
        # the 2022 model for 562.61170, at a LOWER version than its 2025 twin below
        "562.61170.fbamodel": _obj(8, "562.61170.fbamodel", OLD_DATE, 1, NARRATIVE_WS),
    }
    models = {
        "genomeset__562.85951.contigs.fbamodel": _obj(10, "genomeset__562.85951.contigs.fbamodel", NEW_DATE, 1, MODELS_WS),
        # different capitalisation -- also missed by a "genomeset__" prefix test
        "GenomeSetATCC_562.55367.fa_assembly.fbamodel": _obj(11, "GenomeSetATCC_562.55367.fa_assembly.fbamodel", NEW_DATE, 1, MODELS_WS),
        "999.998.fbamodel": _obj(12, "999.998.fbamodel", OLD_DATE, 1, MODELS_WS),
        "randomname.fbamodel": _obj(13, "randomname.fbamodel", OLD_DATE, 1, MODELS_WS),
        # THE TRAP: identical name to the 2022 object, saved 2025, HIGHER version.
        # "latest version" picks this one and calls it old.
        "562.61170.fbamodel": _obj(14, "562.61170.fbamodel", NEW_DATE, 3, MODELS_WS),
    }
    return FakeFBA({NARRATIVE_WS: narrative, MODELS_WS: models})


def _partition():
    return util.partition_models(
        REF_GENOMES, narrative_ws=NARRATIVE_WS, models_ws=MODELS_WS, fba=_make_fba()
    )


def _ref_to_name():
    fba = _make_fba()
    names = {}
    for ws in (NARRATIVE_WS, MODELS_WS):
        for name, info in fba.list_ws_objects(ws).items():
            names[f"{info[6]}/{info[0]}/{info[4]}"] = name
    return names


def test_genomeset_marker_matched_case_insensitively_and_anywhere():
    assert util.is_genomeset_model("genomeset__562.1.contigs.fbamodel")
    assert util.is_genomeset_model("genomeset_5_WS_562.1.fa_assembly.fbamodel")
    assert util.is_genomeset_model("genomeset_WS5_10-25_562.1.fa_assembly.fbamodel")
    assert util.is_genomeset_model("GenomeSetATCC_562.1.fa_assembly.fbamodel")
    assert not util.is_genomeset_model("562.1.fbamodel")
    assert not util.is_genomeset_model("562.1Old2022.fbamodel")


def test_old_arm_never_contains_a_genomeset_model():
    """The reference arm must not be contaminated by 2025 rebuilds."""
    partitions, _ = _partition()
    names = _ref_to_name()
    for refs in partitions["old"].values():
        for ref in refs:
            assert not util.is_genomeset_model(names[ref]), \
                f"genomeset model leaked into the old arm: {names[ref]}"


def test_old_arm_never_contains_a_post_cutoff_model():
    """A plain <gid>.fbamodel saved in 2025 is NOT an old model, whatever it is called."""
    partitions, _ = _partition()
    assert "562.85951" not in partitions["old"]
    assert "562.85951" in partitions["new"]


def test_identical_names_disambiguated_by_save_date_not_version():
    """562.61170: same name, 2022 @ v1 and 2025 @ v3. Highest-version picks WRONG."""
    partitions, _ = _partition()
    old_refs = partitions["old"]["562.61170"]
    assert old_refs == [f"{NARRATIVE_WS}/8/1"], (
        "the old arm must take the EARLIEST-saved object (2022, v1), not the highest "
        f"version (2025, v3) -- got {old_refs}"
    )


def test_one_ref_per_genome_per_arm():
    partitions, _ = _partition()
    for arm in ("old", "new"):
        for gid, refs in partitions[arm].items():
            assert len(refs) == 1, f"{arm}/{gid} has {len(refs)} refs, expected exactly 1"


def test_messy_but_genuine_2022_names_still_land_in_old():
    partitions, _ = _partition()
    assert "562.55367" in partitions["old"]      # 562.55367.fbamodel
    assert "562.55510" in partitions["old"]      # 562.55510Old2022.fbamodel


def test_excludes_and_reports_non_reference_ids():
    partitions, report = _partition()
    excluded = report["excluded"]["old"] + report["excluded"]["new"]
    unmappable = report["unmappable"]["old"] + report["unmappable"]["new"]
    assert "999.999.fbamodel" in excluded
    assert "999.998.fbamodel" in excluded
    assert "notaref.fbamodel" in unmappable
    assert "randomname.fbamodel" in unmappable

    all_old = {r for refs in partitions["old"].values() for r in refs}
    assert f"{NARRATIVE_WS}/3/1" not in all_old   # notaref (unmappable)
    assert f"{NARRATIVE_WS}/4/1" not in all_old   # 999.999 (excluded)
