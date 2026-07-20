#!/usr/bin/env python
"""H100 JOB 2 -- gapfill all three model sets in narrative 265353 to .MMGF.

Every model in the narrative (519 x {_2022p, _2025k, _2026p}) is gapfilled on
Carbon-Pyruvic-Acid -- Clayton's 2022 media, recovered from his own notebook -- via
the OFFICIAL kb_gapfill_metabolic_models, and saved back as <name>.MMGF. Uniform
treatment across all three sets is what makes the 2->3->4 comparison interpretable;
this is why the build job deliberately does not gapfill.

atp_safe=True (Chris's call). Note this ATP-corrects the 2022p models, which Clayton
never did -- so the resimulated-2022p column may replicate his original calls less
closely than the ungapfilled path did. That divergence is a REPORTABLE RESULT, not a
failure: it measures what ATP correction does to the 2022 predictions.

The rxn05759 / rxn03978 directionality fix is applied BEFORE gapfilling and baked
into the saved model. It has to happen first: those reactions leak free electrons
(rxn05759 reverse, charge-imbalanced by 2) and run NH3->NO2- oxidation (rxn03978
forward), and gapfilling will happily lean on that reducing power and add a
different reaction set. ATP correction does NOT catch this -- the leak is electrons,
not ATP, and mass balance passes it because only the charge is wrong.

    python scripts/gapfill_all.py --workers 50
    python scripts/gapfill_all.py --limit 2 --workers 1     # smoke test FIRST
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

NB = Path(__file__).resolve().parent.parent / "notebooks/PRJ-watershed_phenotype_replication"
sys.path.insert(0, str(NB))

TARGET_WS = 265353
SETS = ("_2022p", "_2025k", "_2026p")
GF_SUFFIX = ".MMGF"
FIELDS = ["model_name", "set", "status", "output_name", "directionality_fixed",
          "seconds", "error"]

_KBU = None


def _init():
    global _KBU
    import warnings
    warnings.filterwarnings("ignore")
    from util import ensure_kbase_token
    ensure_kbase_token()
    from kbutillib import KBUtilLib
    last = None
    for attempt in range(5):
        try:
            _KBU = KBUtilLib()
            # kb_gapfill_metabolic_models resolves the genome itself when genome_objs
            # is not supplied, via get_msgenome_from_ontology(native_python_api=
            # self.native_ontology). That defaults to False -> the SDK callback path
            # -> "Either set callback URL", which is fatal outside a container. There
            # is no callback URL here and never will be, so force the native path.
            #
            # Assign on _delegate, NOT the Impl wrapper: MSReconstructionUtilsImpl
            # proxies __getattr__ but not __setattr__, so `_KBU.recon.native_ontology
            # = True` sets the attribute on the wrapper while the delegate keeps
            # False -- it looks like it works and does nothing.
            _KBU.recon._delegate.native_ontology = True
            return
        except Exception as e:                        # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"worker init failed after 5 tries: {last}")


def _gapfill(name):
    from util import (fix_leaky_directionality, GAPFILL_MEDIA, KBASEMEDIA_WS,
                      GROWTH_THRESHOLD)
    t0 = time.time()
    tag = next((s for s in SETS if name.endswith(s)), "?")
    try:
        fba = _KBU.fba
        model = fba.get_model(f"{TARGET_WS}/{name}")           # fresh MSModelUtil per call
        n_fixed = fix_leaky_directionality(model)             # BEFORE gapfilling
        _KBU.recon.kb_gapfill_metabolic_models(
            workspace=str(TARGET_WS),
            model_objs=[model],
            # Must be workspace-qualified: process_media_list qualifies a BARE name
            # with the `workspace` arg (265353), and the media lives in KBaseMedia --
            # 265353 holds no Carbon-Pyruvic-Acid, so the bare name 404s.
            media_list=[f"{KBASEMEDIA_WS}/{GAPFILL_MEDIA}"],  # Carbon-Pyruvic-Acid
            atp_safe=True,
            suffix=GF_SUFFIX,
            save_models_to_kbase=True,
            save_report_to_kbase=False,
            # save_report_to_kbase=False does NOT stop the report being BUILT: the
            # HTML generation is gated on `internal_call`, not on that flag, and it
            # dereferences self.working_dir -- which only KBCallbackUtilsImpl sets,
            # so on MSReconstructionUtils it raises AttributeError AFTER the model
            # has already been saved. That is what makes it worth suppressing rather
            # than tolerating: the save succeeds, the call still raises, and the row
            # gets logged as failed.
            #
            # internal_call=True is the library's own lever for this -- it is what
            # kb_build_metabolic_models passes when it calls gapfilling internally
            # (ms_reconstruction_utils.py:653). It skips only the result-table concat
            # and the report block; neither is used here.
            internal_call=True,
        )
        # Regression guard: "gapfill returned" does NOT mean "model grows".
        # kb_gapfill_metabolic_models returns without error even when it found no
        # solution and added nothing -- which is exactly how all 519 2026p models
        # shipped non-growing and were only caught at the sweep (root cause: the
        # build dropped SK_/DM_ sinks; fixed in build_2026p.py). Verify the SAVED
        # model actually grows on the gapfill media, and mark it failed if not, so
        # the TSV never again reports ok for empty science. Same simulate settings
        # the sweep uses, so "ok here" means "will be scored as growing there".
        saved = fba.get_model(f"{TARGET_WS}/{name}{GF_SUFFIX}")
        gf_media = fba.get_media(GAPFILL_MEDIA, KBASEMEDIA_WS)
        res = fba.simulate_growth_phenotypes(
            saved, {"gapfill": gf_media}, add_missing_exchanges=True, biomass="bio1")
        flux = res["details"][0]["simulated_flux"]
        flux = 0.0 if flux is None or flux != flux else float(flux)
        if flux <= GROWTH_THRESHOLD:
            return {"model_name": name, "set": tag, "status": "failed", "output_name": "",
                    "directionality_fixed": n_fixed, "seconds": round(time.time() - t0, 1),
                    "error": f"gapfill saved a NON-GROWING model "
                             f"(biomass={flux:.4g} on {GAPFILL_MEDIA}; likely no gapfill "
                             f"solution)"}
        return {"model_name": name, "set": tag, "status": "ok",
                "output_name": f"{name}{GF_SUFFIX}", "directionality_fixed": n_fixed,
                "seconds": round(time.time() - t0, 1), "error": ""}
    except Exception as e:                            # noqa: BLE001
        return {"model_name": name, "set": tag, "status": "failed", "output_name": "",
                "directionality_fixed": "", "seconds": round(time.time() - t0, 1),
                "error": f"{type(e).__name__}: {e}"[:300]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=str(NB / "NBOutput/gapfill_all.tsv"))
    ap.add_argument("--workers", type=int, default=50)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--sets", default=",".join(SETS))
    a = ap.parse_args()

    import warnings
    warnings.filterwarnings("ignore")
    from util import get_msfba

    fba = get_msfba()
    objs = fba.list_ws_objects(TARGET_WS, type="KBaseFBA.FBAModel")
    want = tuple(a.sets.split(","))
    models = sorted(k for k in objs
                    if k.endswith(want) and not k.endswith(GF_SUFFIX))
    if a.limit:
        models = models[: a.limit]

    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():
        with out.open() as fh:
            done = {r["model_name"] for r in csv.DictReader(fh, delimiter="\t")
                    if r.get("status") == "ok"}
    todo = [m for m in models if m not in done]
    print(f"gapfill: {len(todo)} models to do ({len(done)} done) x {a.workers} workers")
    for s in want:
        print(f"    {s}: {sum(1 for m in todo if m.endswith(s))}")

    fresh = not out.exists()
    with out.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        if fresh:
            w.writeheader()
        n_ok = n_fail = 0
        with ProcessPoolExecutor(max_workers=a.workers, initializer=_init) as ex:
            futs = {ex.submit(_gapfill, m): m for m in todo}
            for i, f in enumerate(as_completed(futs), 1):
                r = f.result()
                w.writerow(r)
                fh.flush()
                n_ok += r["status"] == "ok"
                n_fail += r["status"] == "failed"
                if i % 25 == 0 or r["status"] == "failed":
                    print(f"  [{i}/{len(todo)}] {r['model_name']} {r['status']} "
                          f"{r['seconds']}s {r['error'][:80]}")
    print(f"\ndone: {n_ok} ok, {n_fail} failed -> {out}")


if __name__ == "__main__":
    main()
