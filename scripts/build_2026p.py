#!/usr/bin/env python
"""H100 JOB 1 -- build the 2026p model set with the OFFICIAL KBase reconstruction.

The 2026p arm = original PATRIC (2022) annotations + CURRENT reconstruction code.
Compared against 2022p (same annotations, 2022 build) it isolates the BUILD effect.

This replaces util.rebuild_from_patric(), which hand-rolled the build as a raw
MSBuilder(genome, template).build(...) with NO ATP safety. That arm came out 97%
autotrophic with 10 free-lunch models and only 7 clean -- by far the worst of the
three -- because it was never the reconstruction pipeline, just a piece of it. The
2025k models, built by real KBase with atp_safe, had ZERO free-lunch models.

Why build_metabolic_model() and not kb_build_metabolic_models():
    kb_build_metabolic_models is the KBase SDK-APP wrapper. It fetches the genome via
    get_msgenome_from_ontology -> anno_client, which REQUIRES an SDK callback URL; the
    advertised native_python_api fallback is threaded through three call sites and then
    ignored by anno_client, so it does not exist. Outside an SDK container _callback_url
    is None and the call raises. build_metabolic_model is the callback-free core that
    wrapper delegates to (ms_reconstruction_utils.py:568) -- same official MSBuilder,
    same official MSATPCorrection, same template selection. We hand it a genome fetched
    directly and save the result ourselves, which is all the wrapper adds.

gs_template stays 'auto' deliberately: "current reconstruction code" means the template
it picks TODAY (GramNegModelTemplateV6) rather than the 2022 V4. The template shift is
part of the build effect we are trying to measure, not a confound to pin away.

Gapfilling is NOT run here. All three sets are gapfilled uniformly in job 2, so the
comparison is not confounded by differing gapfills.

    python scripts/build_2026p.py --workers 50
    python scripts/build_2026p.py --limit 2 --workers 1    # smoke test
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

NB = Path(__file__).resolve().parent.parent / "notebooks/PRJ-watershed_phenotype_replication"
sys.path.insert(0, str(NB))

TARGET_WS = 265353
SUFFIX = "_2026p"
FIELDS = ["genome_id", "status", "model_name", "genome_ref", "seconds", "error"]

_KBU = None
_CLASSIFIER = None
_CORE_TEMPLATE = None


def _init():
    """Pool initializer. Retries: a KBase blip here kills every worker at once.

    The classifier and core template are fetched ONCE per worker, not per genome --
    at 50 workers x 519 genomes the per-genome path trips KBase's rate limit.
    """
    global _KBU, _CLASSIFIER, _CORE_TEMPLATE
    import warnings
    warnings.filterwarnings("ignore")
    from util import ensure_kbase_token
    ensure_kbase_token()
    from kbutillib import KBUtilLib
    last = None
    for attempt in range(5):
        try:
            _KBU = KBUtilLib()
            _CLASSIFIER = _KBU.recon.get_classifier()
            _CORE_TEMPLATE = _KBU.recon.get_template(_KBU.recon.templates["core"], None)
            return
        except Exception as e:                       # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"worker init failed after 5 tries: {last}")


def _build(genome_id):
    t0 = time.time()
    try:
        from util import get_patric_genome_ref
        fba = _KBU.fba
        gref = get_patric_genome_ref(str(genome_id), fba=fba)
        if not gref:
            return {"genome_id": genome_id, "status": "no_genome", "model_name": "",
                    "genome_ref": "", "seconds": 0, "error": "no PATRIC genome ref"}
        name = f"{genome_id}{SUFFIX}"
        # get_msgenome: the PATRIC annotations AS THEY ARE. Not the ontology path
        # (needs an SDK callback) and not RAST -- re-annotating here would reproduce
        # the very drift this arm exists to isolate.
        genome = fba.get_msgenome(gref)
        _, mdlutl = _KBU.recon.build_metabolic_model(
            genome=genome,
            genome_classifier=_CLASSIFIER,
            model_id=name,
            gs_template="auto",          # -> GramNegModelTemplateV6 for E. coli
            core_template=_CORE_TEMPLATE,
            atp_safe=True,               # the whole point -- what my build skipped
        )
        # KBaseFBA.FBAModel requires a non-null genome_ref. The SDK wrapper sets this
        # from genome.annoont.info (ms_reconstruction_utils.py:618) -- the ontology
        # object we cannot build without a callback. We already hold the same ref.
        #
        # Capture the template ref BEFORE conversion. build_metabolic_model attaches
        # the chosen GS template as mdlutl.model.template (for E. coli, "auto" ->
        # GramNegModelTemplateV6 -> 12998/24/4). The conversion below drops it exactly
        # as it drops genome_ref, so save_model would serialize template_ref="" and
        # gapfill would then raise at get_template('') with "reference cannot be null
        # or the empty string" -- which failed 100% of the 2026p arm on the first run.
        _tmpl = getattr(mdlutl.model, "template", None)
        tmpl_ref = getattr(getattr(_tmpl, "info", None), "reference", None)
        # Capture the SK_/DM_ (sink/demand) reactions BEFORE conversion too.
        # CobraModelConverter.build() silently drops EVERY sink and demand reaction,
        # and several of them are load-bearing: the byproduct sinks SK_cpd03091
        # (5'-deoxyadenosine) and SK_cpd02701 (S-adenosyl-4-methylthio-2-oxobutanoate)
        # drain dead-end products of radical-SAM cofactor-biosynthesis reactions, so
        # without them those biosynthesis reactions cannot carry flux (mass balance
        # forces them to zero), the biomass cofactors (ACP, pyridoxal-P, ubiquinone,
        # ...) become unmakeable, and gapfill is INFEASIBLE -> "No gapfilling solution
        # found" -> a saved model that does not grow. Confirmed causal: re-adding these
        # flips gapfill from infeasible back to feasible. This is why the first 2026p
        # sweep produced zero growth across all 519 models.
        _sink_demand = [r.copy() for r in mdlutl.model.reactions
                        if r.id.startswith(("SK_", "DM_"))]
        # Convert BEFORE setting the refs: MSBuilder returns a plain cobra model, and
        # save_model swaps in CobraModelConverter(...).build() when the model is not
        # already an FBAModel -- which would silently drop a genome_ref set here.
        # Doing the conversion ourselves makes save_model's isinstance check pass, so
        # the refs survive to get_data().
        if not isinstance(mdlutl.model, _KBU.recon.FBAModel):
            mdlutl.model = _KBU.recon.CobraModelConverter(mdlutl.model).build()
            existing = {r.id for r in mdlutl.model.reactions}
            readd = [r for r in _sink_demand if r.id not in existing]
            if readd:
                mdlutl.model.add_reactions(readd)
        mdlutl.model.genome_ref = gref
        if tmpl_ref:                          # same post-conversion re-set as genome_ref
            mdlutl.model.template_ref = tmpl_ref
            mdlutl.model.template_refs = [tmpl_ref]
        # save_model reads self.provenance; the kb_* SDK wrappers set it for you, so
        # bypassing them means setting it here. Wanted anyway: the collaborators get
        # a model that records which genome and which code produced it.
        _KBU.recon.set_provenance(
            method="build_metabolic_model",
            description=f"2026p arm: current ModelSEED reconstruction (atp_safe) "
                        f"on the original PATRIC annotation of {genome_id}",
            input_objects=[gref],
            params={"atp_safe": True, "gs_template": "auto", "run_gapfilling": False},
            service="WatershedPhenotypeReplication",
        )
        _KBU.recon.save_model(mdlutl, str(TARGET_WS), name)
        return {"genome_id": genome_id, "status": "ok", "model_name": name,
                "genome_ref": gref, "seconds": round(time.time() - t0, 1), "error": ""}
    except Exception as e:                           # noqa: BLE001
        return {"genome_id": genome_id, "status": "failed", "model_name": "",
                "genome_ref": "", "seconds": round(time.time() - t0, 1),
                "error": f"{type(e).__name__}: {e}"[:300]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=str(NB / "NBOutput/build_2026p.tsv"))
    ap.add_argument("--workers", type=int, default=50)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    import warnings
    warnings.filterwarnings("ignore")
    from util import load_watershed_sheet

    ref, _, _ = load_watershed_sheet()
    genomes = sorted(map(str, ref))
    if a.limit:
        genomes = genomes[: a.limit]

    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():                                  # resumable
        with out.open() as fh:
            done = {r["genome_id"] for r in csv.DictReader(fh, delimiter="\t")
                    if r.get("status") in ("ok", "no_genome")}
    todo = [g for g in genomes if g not in done]
    print(f"2026p build: {len(todo)} to do ({len(done)} already done) x {a.workers} workers")

    fresh = not out.exists()
    with out.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        if fresh:
            w.writeheader()
        n_ok = n_fail = 0
        with ProcessPoolExecutor(max_workers=a.workers, initializer=_init) as ex:
            futs = {ex.submit(_build, g): g for g in todo}
            for i, f in enumerate(as_completed(futs), 1):
                r = f.result()
                w.writerow(r)
                fh.flush()
                n_ok += r["status"] in ("ok", "no_genome")
                n_fail += r["status"] == "failed"
                if i % 25 == 0 or r["status"] == "failed":
                    print(f"  [{i}/{len(todo)}] {r['genome_id']} {r['status']} "
                          f"{r['seconds']}s {r['error'][:80]}")
    print(f"\ndone: {n_ok} ok, {n_fail} failed -> {out}")


if __name__ == "__main__":
    main()
