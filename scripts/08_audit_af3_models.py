#!/usr/bin/env python3
"""
Step 4B: All-30-models AlphaFold 3 confidence, pLDDT profile, and Motif-IX
R <-> KYFLS aromatic lock audit across all 6 trans complexes (model_0..4)
plus the cis Arch 3 reference (AF-A0A6N9NLZ4-F1).

Runs 100% self-contained from repository data/structures/af3_models/.
"""
from __future__ import annotations

import json
import re
import numpy as np

import lib_c5mtase as lib


COMPLEX_ORDER = [
    "phagep1_arch1_dmt_trans_64aa_core",
    "phagep1_arch1_dmt_trans_157aa",
    "bact_arch1_asoli_trans",
    "arch_arch1_nanob_trans",
    "posctrl_arch2_eaei_trans",
    "negctrl_hhai_asoli_trans",
]


def main():
    sdir = lib.STRUCT_DIR
    af3_dir = sdir / "af3_models"

    print("=" * 112)
    print("PART 1: ALL-5-MODELS CONFIDENCE AUDIT (6 COMPLEXES x 5 DIFFUSION MODELS = 30 MODELS)")
    print("=" * 112)
    print(f"{'Complex':36s} | {'ipTM (mean±sd [min-max])':26s} | {'pTM':6s} | {'rank':6s} | {'PAE_min':11s} | {'clash':5s}")
    print("-" * 112)
    for name in COMPLEX_ORDER:
        iptms, ptms, rscores, pae_mins, clashes = [], [], [], [], []
        for k in range(5):
            d = json.loads((af3_dir / f"fold_{name}_summary_confidences_{k}.json").read_text())
            iptms.append(d["iptm"])
            ptms.append(d["ptm"])
            rscores.append(d["ranking_score"])
            clashes.append(d.get("has_clash", 0))
            cpm = d.get("chain_pair_pae_min", [[0, 99], [99, 0]])
            pae_mins.append(min(cpm[0][1], cpm[1][0]))
        print(
            f"{name:36s} | {np.mean(iptms):.3f}±{np.std(iptms):.3f} [{min(iptms):.2f}-{max(iptms):.2f}]     | "
            f"{np.mean(ptms):.3f}  | {np.mean(rscores):.3f}  | {np.mean(pae_mins):4.2f}±{np.std(pae_mins):.2f} A | {sum(clashes)}/5"
        )

    print("\n" + "=" * 112)
    print("PART 2: MOTIF-IX ARGININE <-> KYFLS AROMATIC LOCK & CHAIN-B pLDDT PROFILE")
    print("=" * 112)
    for name in COMPLEX_ORDER:
        ch0 = lib.parse_cif_atoms(af3_dir / f"fold_{name}_model_0.cif")
        idsA, seqA, caA0, _ = lib.chain_ca_sequence(ch0["A"])
        idsB, seqB, caB0, bfB0 = lib.chain_ca_sequence(ch0["B"])
        m9 = re.search(r"R[RK][LMIFV][TIVSLH]P.{2,4}R[LVIM][QMV]", seqA) or re.search(r"RKLHPRECAR", seqA)
        m_kyf = re.search(r"[KRE][YFS][FYW][LF][SN]", seqB)
        m_csl = re.search(r"[CE][STLV][LF]S[QDR][IVL][LVIM][EDK]", seqB)

        dists = []
        for k in range(5):
            chk = lib.parse_cif_atoms(af3_dir / f"fold_{name}_model_{k}.cif")
            _, _, caA, _ = lib.chain_ca_sequence(chk["A"])
            _, _, caB, _ = lib.chain_ca_sequence(chk["B"])
            if m9 and m_kyf:
                dists.append(float(np.linalg.norm(caA[m9.start()] - caB[m_kyf.start() + 2])))
            elif name == "posctrl_arch2_eaei_trans" and m_kyf:
                m9b = re.search(r"RRMLPSEYERLQ", seqB)
                dists.append(float(np.linalg.norm(caB[m9b.start()] - caB[m_kyf.start() + 2])))

        csl_plddt = float(bfB0[m_csl.start():min(len(bfB0), m_csl.end() + 15)].mean()) if m_csl else float("nan")
        split_pt = max(0, m_csl.start() - 25) if m_csl else 0
        n_plddt = float(bfB0[:split_pt].mean()) if split_pt > 5 else float("nan")
        print(
            f"  {name:36s} | R1-F/Y (5m): {[round(x, 2) for x in dists]} (mean={np.mean(dists):.2f}±{np.std(dists):.2f} A) | "
            f"pLDDT N-tail={n_plddt:4.1f}, CSL-KYF={csl_plddt:4.1f}"
        )

    # Cis reference check
    a0 = lib.parse_cif_atoms(sdir / "AF-A0A6N9NLZ4-F1-model_v6.cif")
    ids0, seq0, ca0, _ = lib.chain_ca_sequence(a0["A"])
    m9_0 = re.search(r"RRFTPLECERLM", seq0)
    kyf_0 = re.search(r"KFYLS", seq0)
    d1 = float(np.linalg.norm(ca0[m9_0.start()] - ca0[kyf_0.start() + 2]))
    d2 = float(np.linalg.norm(ca0[m9_0.start() + 1] - ca0[kyf_0.start() + 2]))
    print(f"\n  CIS REF AF-A0A6N9NLZ4-F1: R{ids0[m9_0.start()]} <-> Y{ids0[kyf_0.start()+2]} = {d1:.2f} A (R2={d2:.2f} A)")


if __name__ == "__main__":
    main()
