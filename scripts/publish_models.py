#!/usr/bin/env python
"""Copy the 2022p and 2025k model sets into the deliverable narrative (265353).

Why this exists: narrative 119455 holds BOTH model generations, and for 18 genomes
the 2022 model was OVERWRITTEN in place by a 2025 rebuild -- same object, v1=2022,
v2=2025. Resolving by head version silently yields the 2025 model and those genomes
get dropped from the 2022p arm entirely (the old sweep ran 501, not 519).

So we resolve every genome's 2022 model by walking the object HISTORY and taking the
newest version saved before OLD_MODEL_CUTOFF, then copy that exact VERSIONED ref.
That recovers all 18 and gives a narrative holding a complete 519 x 2 corpus:

    <genome_id>_2022p    PATRIC annotations, KBase 2022 build   (Clayton's originals)
    <genome_id>_2025k    KBase re-annotation, KBase 2025 build

The 2026p set is built separately (H100 job 1, kb_build_metabolic_models), and all
three sets are gapfilled to .MMGF separately (H100 job 2).

    python scripts/publish_models.py --dry-run      # resolve + report, copy nothing
    python scripts/publish_models.py                # do it
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

NB = Path(__file__).resolve().parent.parent / "notebooks/PRJ-watershed_phenotype_replication"
sys.path.insert(0, str(NB))

TARGET_WS = 265353
OLD_CUTOFF = "2023-01-01"


def _hist(ws, wsid, objid):
    return ws.get_object_history({"wsid": wsid, "objid": objid})


def resolve(fba, ws, ref_genomes):
    """-> {genome_id: {"2022p": versioned_ref | None, "2025k": versioned_ref | None}}"""
    from util import partition_models, NARRATIVE_WS

    parts, _ = partition_models(ref_genomes, fba=fba)
    objs = fba.list_ws_objects(NARRATIVE_WS, type="KBaseFBA.FBAModel")
    out, recovered = {}, []

    for gid in sorted(map(str, ref_genomes)):
        entry = {"2022p": None, "2025k": None}

        # ---- 2025k: partition's "new" head is correct as-is -------------------
        new = parts["new"].get(gid)
        if new:
            entry["2025k"] = new[-1] if isinstance(new, list) else new

        # ---- 2022p: walk history, newest version saved BEFORE the cutoff ------
        old = parts["old"].get(gid)
        if old:
            entry["2022p"] = old[-1] if isinstance(old, list) else old
        else:
            # partition dropped it -> the plain object's head is a 2025 overwrite.
            info = objs.get(f"{gid}.fbamodel")
            if info:
                objid = info[0]
                cand = [h for h in _hist(ws, NARRATIVE_WS, objid)
                        if str(h[3])[:10] < OLD_CUTOFF]
                if cand:
                    v = max(cand, key=lambda h: str(h[3]))[4]
                    entry["2022p"] = f"{NARRATIVE_WS}/{objid}/{v}"
                    recovered.append(gid)

        out[gid] = entry
    return out, recovered


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--target", type=int, default=TARGET_WS)
    a = ap.parse_args()

    import warnings
    warnings.filterwarnings("ignore")
    from util import get_msfba, load_watershed_sheet

    fba = get_msfba()
    ws = fba.ws_client() if callable(fba.ws_client) else fba.ws_client
    ref, _, _ = load_watershed_sheet()

    plan, recovered = resolve(fba, ws, ref)
    n22 = sum(1 for e in plan.values() if e["2022p"])
    n25 = sum(1 for e in plan.values() if e["2025k"])
    print(f"resolved {len(plan)} genomes: 2022p={n22}/{len(plan)}  2025k={n25}/{len(plan)}")
    print(f"recovered via versioned ref (would otherwise be DROPPED): {len(recovered)}")
    if recovered:
        print("  " + ", ".join(recovered[:6]) + (" ..." if len(recovered) > 6 else ""))

    if a.dry_run:
        for gid in list(plan)[:3]:
            print(f"  {gid}: 2022p={plan[gid]['2022p']}  2025k={plan[gid]['2025k']}")
        print("\nDRY RUN -- nothing copied.")
        return

    done = fail = 0
    for gid, e in plan.items():
        for tag, ref_str in (("2022p", e["2022p"]), ("2025k", e["2025k"])):
            if not ref_str:
                continue
            try:
                ws.copy_object({"from": {"ref": ref_str},
                                "to": {"wsid": a.target, "name": f"{gid}_{tag}"}})
                done += 1
            except Exception as ex:
                print(f"  FAIL {gid}_{tag} <- {ref_str}: {type(ex).__name__}: {ex}")
                fail += 1
    print(f"\ncopied {done} objects into ws {a.target} ({fail} failed)")


if __name__ == "__main__":
    main()
