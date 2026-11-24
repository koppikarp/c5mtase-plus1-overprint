#!/usr/bin/env python3
"""
Step 4: Quantitative 3D Structural Audit across ALL 30 AlphaFold 3 diffusion models
(6 trans complexes x 5 models) vs. cis Arch 3 reference (AF-A0A6N9NLZ4-F1) and
DNA/SAH-bound M.HhaI ternary crystal structure (PDB 1MHT).

Runs 100% self-contained from repository data/structures/:
  - data/structures/AF-A0A6N9NLZ4-F1-model_v6.cif
  - data/structures/1MHT.cif
  - data/structures/af3_models/fold_*_model_{0..4}.cif
  - data/structures/af3_models/fold_*_summary_confidences_{0..4}.json

Key methodological fix (rigid-core superposition):
  Motif IV anchor is strictly restricted to the rigid pre-catalytic beta-strand
  M4[-4..+1] (EXCLUDING the mobile base-flipping catalytic loop M4[+2..+8], which
  displaces ~27 A upon DNA binding in 1MHT and in apo diffusion models!).

Exports data/structures/trans_vs_cis_3d_audit.json with 5-model mean +/- SD for
every metric (ZERO hardcoded reference or model values!).
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import lib_c5mtase as lib


DNA_COMPS = {"DC", "DG", "DA", "DT", "5FC", "C", "G", "A", "T"}
# Strictly rigid Rossmann core anchors (excludes mobile Motif IV base-flipping loop M4[+2..+8]!)
RIGID_ANCHORS = [
    ("M4", -4, 2),    # 6 CA: rigid strand [-4..+1] around P of PCQ
    ("M6", -3, 7),    # 10 CA: Motif VI ENV core
    ("M8", -3, 9),    # 12 CA: Motif VIII QxRxR core
    ("M9", 0, 10),    # 10 CA: Motif IX (skipped automatically when absent, e.g. M.EaeI alpha)
    ("M10", -1, 7),   # 8 CA: Motif X GNxxV helix
]

COMPLEX_ORDER = [
    "phagep1_arch1_dmt_trans_64aa_core",
    "phagep1_arch1_dmt_trans_157aa",
    "bact_arch1_asoli_trans",
    "arch_arch1_nanob_trans",
    "posctrl_arch2_eaei_trans",
    "negctrl_hhai_asoli_trans",
]


def find_landmarks_strict(seq: str) -> Dict[str, int]:
    lm: Dict[str, int] = {}
    m4 = re.search(r"[GASCVITLFY][FTSAVILYCMWP]PC[QSNTG]", seq)
    if m4:
        lm["M4"] = m4.start() + 2  # 0-based index of P in xxPCQ
    m6 = re.search(r"[IVLMFWYC]ENV", seq)
    if m6:
        lm["M6"] = m6.start() + 1
    m8 = re.search(r"[PQ][HKQR][NRK]R[EQK]R", seq) or re.search(r"[PQ][QKHR][NKHRQ]R[EKPQ]R[VILYFM]", seq)
    if m8:
        lm["M8"] = m8.start()
    m9 = (
        re.search(r"R[RK][LMIFV][TIVSLH]P.{2,4}R[LVIM][QMV]", seq)
        or re.search(r"[RK][LMIV][TSIVL][PRK].{1,3}[CEWYFL].[RK][LM][QEH]", seq)
        or re.search(r"RKLHPRECAR", seq)
    )
    if m9:
        lm["M9"] = m9.start()
    m10 = re.search(r"GN[SA][VWILMT][ATVIL]", seq)
    if m10:
        lm["M10"] = m10.start()
    return lm


def collect_rigid_anchors(
    seq1: str, c1: np.ndarray, lm1: Dict[str, int],
    seq2: str, c2: np.ndarray, lm2: Dict[str, int],
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    P, Q, used = [], [], []
    for key, lo, hi in RIGID_ANCHORS:
        if key in lm1 and key in lm2:
            i1, i2 = lm1[key], lm2[key]
            for off in range(lo, hi):
                if 0 <= i1 + off < len(seq1) and 0 <= i2 + off < len(seq2):
                    P.append(c1[i1 + off])
                    Q.append(c2[i2 + off])
            used.append(key)
    return np.asarray(P, dtype=float), np.asarray(Q, dtype=float), used


def chain_heavy_coords(atoms: List[dict], res_lo: int | None = None, res_hi: int | None = None) -> np.ndarray:
    pts = [
        a["xyz"]
        for a in atoms
        if a["group"] == "ATOM"
        and not a["atom"].startswith("H")
        and a["comp"] in lib.AA3TO1
        and (res_lo is None or a["resseq"] >= res_lo)
        and (res_hi is None or a["resseq"] <= res_hi)
    ]
    return np.asarray(pts, dtype=float) if pts else np.zeros((0, 3))


def min_dist(A: np.ndarray, B: np.ndarray) -> float:
    if len(A) == 0 or len(B) == 0:
        return float("nan")
    d2 = np.sum((A[:, None, :] - B[None, :, :]) ** 2, axis=2)
    return float(np.sqrt(np.min(d2)))


def main():
    sdir = lib.STRUCT_DIR
    af3_dir = sdir / "af3_models"

    # 1. Load cis reference AF-A0A6N9NLZ4-F1 (485 aa)
    a0_chains = lib.parse_cif_atoms(sdir / "AF-A0A6N9NLZ4-F1-model_v6.cif")
    a0_atoms = a0_chains["A"]
    a0_ids, a0_seq, a0_ca, a0_bf = lib.chain_ca_sequence(a0_atoms)
    a0_lm = find_landmarks_strict(a0_seq)
    m_cis_csl = re.search(r"DLKDIVE", a0_seq)
    m_cis_kyf = re.search(r"KFYLS", a0_seq)
    assert m_cis_csl and m_cis_kyf
    cis_csl_0 = m_cis_csl.start()  # 0-based index of D394 (1-based 394 -> 0-based 393)
    cis_core_lo = cis_csl_0 - 6    # 0-based 387 -> 1-based 388 (I388..E441, 54 aa)
    cis_core_hi = min(len(a0_seq), cis_csl_0 + 48)
    cis_core_ca = a0_ca[cis_core_lo:cis_core_hi]
    cis_core_com = cis_core_ca.mean(axis=0)
    cis_core_heavy = chain_heavy_coords(a0_atoms, a0_ids[cis_core_lo], a0_ids[cis_core_hi - 1])

    # Exact cis inter-domain contacts with explicit residue ranges!
    mt_ca = a0_ca[:365]  # 1..365 (full C5-MTase domain Motifs I-X)
    ct_core_ca = a0_ca[cis_core_lo:cis_core_hi]  # 388..441 (54-aa CSLSDILE/KYFLS core)
    ct_full_ca = a0_ca[365:]  # 366..485 (full C-terminal domain)
    ca_c_core = int(np.sum(np.linalg.norm(mt_ca[:, None, :] - ct_core_ca[None, :, :], axis=2) < 8.0))
    ca_c_full = int(np.sum(np.linalg.norm(mt_ca[:, None, :] - ct_full_ca[None, :, :], axis=2) < 8.0))
    mt_hv = chain_heavy_coords(a0_atoms, 1, 365)
    hv_c_core = int(np.sum(np.linalg.norm(mt_hv[:, None, :] - cis_core_heavy[None, :, :], axis=2) < 4.5))
    ct_full_hv = chain_heavy_coords(a0_atoms, 366, 485)
    hv_c_full = int(np.sum(np.linalg.norm(mt_hv[:, None, :] - ct_full_hv[None, :, :], axis=2) < 4.5))

    cis_r1_idx = a0_lm["M9"]
    cis_y_idx = m_cis_kyf.start() + 2
    cis_r1_y = float(np.linalg.norm(a0_ca[cis_r1_idx] - a0_ca[cis_y_idx]))
    cis_r2_y = float(np.linalg.norm(a0_ca[cis_r1_idx + 1] - a0_ca[cis_y_idx]))
    cis_cys_idx = a0_lm["M4"] + 1  # Cys79 (1-based: a0_ids[cis_cys_idx])
    cis_cys_resnum = a0_ids[cis_cys_idx]

    # 2. Superpose 1MHT (M.HhaI-DNA-SAH ternary crystal) onto AF-A0A6N9NLZ4-F1 using RIGID CORE ONLY!
    mht_chains = lib.parse_cif_atoms(sdir / "1MHT.cif")
    mht_ids, mht_seq, mht_ca, _ = lib.chain_ca_sequence(mht_chains["A"])
    mht_lm = find_landmarks_strict(mht_seq)
    P_mht, Q_mht, used_mht = collect_rigid_anchors(mht_seq, mht_ca, mht_lm, a0_seq, a0_ca, a0_lm)
    R_mht, t_mht, rmsd_mht = lib.kabsch_align(P_mht, Q_mht)

    dna_raw = [
        a["xyz"]
        for ch_atoms in mht_chains.values()
        for a in ch_atoms
        if a["comp"] in DNA_COMPS and not a["atom"].startswith("H")
    ]
    sah_raw = [
        a["xyz"]
        for ch_atoms in mht_chains.values()
        for a in ch_atoms
        if a["comp"] == "SAH" and not a["atom"].startswith("H")
    ]
    dna_xyz = np.asarray(dna_raw, dtype=float) @ R_mht.T + t_mht
    sah_xyz = np.asarray(sah_raw, dtype=float) @ R_mht.T + t_mht

    cis_dna_ca_min = min_dist(cis_core_ca, dna_xyz)
    cis_dna_hv_min = min_dist(cis_core_heavy, dna_xyz)
    cis_sah_ca_min = min_dist(cis_core_ca, sah_xyz)
    cis_sah_hv_min = min_dist(cis_core_heavy, sah_xyz)
    cis_cys_ca_min = float(np.min(np.linalg.norm(cis_core_ca - a0_ca[cis_cys_idx], axis=1)))
    mht_cys_resnum = mht_ids[mht_lm["M4"] + 1]  # Cys81 in 1MHT!

    print("=" * 110)
    print(f"RIGID-CORE 1MHT -> AF-A0A6N9NLZ4-F1 SUPERPOSITION ({used_mht}, {len(P_mht)} CA): RMSD = {rmsd_mht:.2f} A!")
    print(f"  1MHT active-site Cys = Cys{mht_cys_resnum} | AF-A0A6N9NLZ4-F1 active-site Cys = Cys{cis_cys_resnum}")
    print(f"  AF-A0A6N9NLZ4-F1 cis core [{a0_ids[cis_core_lo]}..{a0_ids[cis_core_hi-1]}]:")
    print(f"    Motif-IX R{a0_ids[cis_r1_idx]} <-> KFYLS Y{a0_ids[cis_y_idx]} = {cis_r1_y:.2f} A (R2={cis_r2_y:.2f} A)")
    print(f"    Contacts MTase[1..365] vs core[388..441]: {ca_c_core} CA(<8A), {hv_c_core} heavy(<4.5A) | vs full[366..485]: {ca_c_full} CA, {hv_c_full} heavy")
    print(f"    Dist to 1MHT DNA: CA min = {cis_dna_ca_min:.2f} A, heavy min = {cis_dna_hv_min:.2f} A | SAH: CA={cis_sah_ca_min:.2f} A, heavy={cis_sah_hv_min:.2f} A | Cys{cis_cys_resnum}={cis_cys_ca_min:.2f} A")
    print("=" * 110)

    # 3. Audit ALL 6 complexes x ALL 5 models (30 models total!)
    trans_records = []
    for name in COMPLEX_ORDER:
        iptms, ptms, rscores, pae_mins, clashes, fdis = [], [], [], [], [], []
        core_rmsds, com_dists, docked_rmsds, fold_rmsds = [], [], [], []
        r1_kyf_dists, dna_ca_mins, dna_hv_mins, sah_ca_mins, cys_ca_mins = [], [], [], [], []
        hv_contacts_list, ca_contacts_list = [], []
        plddt_b_m0 = []
        len_A = len_B = 0
        csl_pos_B = kyf_pos_B = None
        chainA_cys_resnum = None
        r1_label = kyf_label = None

        for k in range(5):
            sp = af3_dir / f"fold_{name}_summary_confidences_{k}.json"
            cp = af3_dir / f"fold_{name}_model_{k}.cif"
            sd = json.loads(sp.read_text())
            iptms.append(float(sd["iptm"]))
            ptms.append(float(sd["ptm"]))
            rscores.append(float(sd["ranking_score"]))
            clashes.append(int(sd.get("has_clash", 0)))
            fdis.append(float(sd.get("fraction_disordered", 0.0)))
            cpm = sd.get("chain_pair_pae_min", [[0, 99], [99, 0]])
            pae_mins.append(float(min(cpm[0][1], cpm[1][0])))

            ch = lib.parse_cif_atoms(cp)
            idsA, seqA, caA, bfA = lib.chain_ca_sequence(ch["A"])
            idsB, seqB, caB, bfB = lib.chain_ca_sequence(ch["B"])
            if k == 0:
                len_A, len_B = len(seqA), len(seqB)
                plddt_b_m0 = [round(float(x), 1) for x in bfB]

            lmA = find_landmarks_strict(seqA)
            if k == 0 and "M4" in lmA:
                chainA_cys_resnum = idsA[lmA["M4"] + 1]

            # Rigid-core Kabsch alignment of Chain A onto AF-A0A6N9NLZ4-F1
            PA, QA, usedA = collect_rigid_anchors(seqA, caA, lmA, a0_seq, a0_ca, a0_lm)
            RA, tA, rmsdA = lib.kabsch_align(PA, QA)
            core_rmsds.append(round(rmsdA, 3))

            # Transform Chain B into AF-A0A6N9NLZ4-F1 + 1MHT frame (NO Chain B re-alignment!)
            caB_tr = caB @ RA.T + tA
            hvA = chain_heavy_coords(ch["A"])
            hvB_raw = chain_heavy_coords(ch["B"])
            hvB_tr = hvB_raw @ RA.T + tA

            # Locate conserved CSLSDILE + KYFLS core in Chain B
            m_csl = re.search(r"[CE][STLV][LF]S[QDR][IVL][LVIM][EDK]", seqB)
            m_kyf = re.search(r"[KRE][YFS][FYW][LF][SN]", seqB)
            assert m_csl and m_kyf
            if k == 0:
                csl_pos_B = m_csl.start() + 1
                kyf_pos_B = m_kyf.start() + 1

            # Match residues [-6..+47] relative to CSLSDILE start (capped to chain B length)
            b_csl0 = m_csl.start()
            off_lo = max(-6, -b_csl0)
            off_hi = min(48, len(seqB) - b_csl0, len(a0_seq) - cis_csl_0)
            b_idx = [b_csl0 + o for o in range(off_lo, off_hi)]
            c_idx = [cis_csl_0 + o for o in range(off_lo, off_hi)]

            b_core_tr = caB_tr[b_idx]
            c_core_ref = a0_ca[c_idx]
            com_d = float(np.linalg.norm(b_core_tr.mean(axis=0) - c_core_ref.mean(axis=0)))
            dock_rmsd = float(np.sqrt(np.mean(np.sum((b_core_tr - c_core_ref) ** 2, axis=1))))
            _, _, fold_r = lib.kabsch_align(b_core_tr, c_core_ref)
            com_dists.append(round(com_d, 3))
            docked_rmsds.append(round(dock_rmsd, 3))
            fold_rmsds.append(round(fold_r, 3))

            # Distance from docked Chain B conserved core to superimposed 1MHT DNA, SAH, and Chain A active-site Cys
            b_core_hv_tr = chain_heavy_coords(ch["B"], idsB[b_idx[0]], idsB[b_idx[-1]]) @ RA.T + tA
            dna_ca_mins.append(round(min_dist(b_core_tr, dna_xyz), 3))
            dna_hv_mins.append(round(min_dist(b_core_hv_tr, dna_xyz), 3))
            sah_ca_mins.append(round(min_dist(b_core_tr, sah_xyz), 3))
            if "M4" in lmA:
                cys_ca_mins.append(round(float(np.min(np.linalg.norm(caB[b_idx] - caA[lmA["M4"] + 1], axis=1))), 3))

            # Inter-chain contacts
            hv_c = int(np.sum(np.linalg.norm(hvA[:, None, :] - hvB_raw[None, :, :], axis=2) < 4.5))
            ca_c = int(np.sum(np.linalg.norm(caA[:, None, :] - caB[None, :, :], axis=2) < 8.0))
            hv_contacts_list.append(hv_c)
            ca_contacts_list.append(ca_c)

            # Motif-IX R1 <-> KYFLS F/Y distance
            if "M9" in lmA:
                r1_idx = lmA["M9"]
                f_idx = m_kyf.start() + 2
                r1_kyf_dists.append(round(float(np.linalg.norm(caA[r1_idx] - caB[f_idx])), 3))
                if k == 0:
                    r1_label = f"A:{seqA[r1_idx]}{idsA[r1_idx]}"
                    kyf_label = f"B:{seqB[f_idx]}{idsB[f_idx]}"
            elif name == "posctrl_arch2_eaei_trans":
                m9b = re.search(r"RRMLPSEYERLQ", seqB)
                r1_idx = m9b.start()
                f_idx = m_kyf.start() + 2
                r1_kyf_dists.append(round(float(np.linalg.norm(caB[r1_idx] - caB[f_idx])), 3))
                if k == 0:
                    r1_label = f"B:{seqB[r1_idx]}{idsB[r1_idx]} (intra-beta)"
                    kyf_label = f"B:{seqB[f_idx]}{idsB[f_idx]}"

        rec = {
            "name": name,
            "len_A": len_A, "len_B": len_B,
            "chainA_active_cys_resnum": chainA_cys_resnum,
            "r1_partner": r1_label, "kyf_aromatic": kyf_label,
            "iptm_5": iptms, "iptm_mean": round(float(np.mean(iptms)), 4), "iptm_sd": round(float(np.std(iptms)), 4),
            "ptm_5": ptms, "ptm_mean": round(float(np.mean(ptms)), 4), "ptm_sd": round(float(np.std(ptms)), 4),
            "ranking_score_5": rscores, "ranking_score_mean": round(float(np.mean(rscores)), 4), "ranking_score_sd": round(float(np.std(rscores)), 4),
            "pae_min_5": pae_mins, "pae_min_mean": round(float(np.mean(pae_mins)), 3), "pae_min_sd": round(float(np.std(pae_mins)), 3),
            "clashes": int(sum(clashes)),
            "chainA_rigid_rmsd_5": core_rmsds,
            "chainA_rigid_rmsd_mean": round(float(np.mean(core_rmsds)), 3),
            "r1_kyf_dists_5": r1_kyf_dists,
            "r1_kyf_mean": round(float(np.mean(r1_kyf_dists)), 3), "r1_kyf_sd": round(float(np.std(r1_kyf_dists)), 3),
            "com_dist_vs_cis_5": com_dists,
            "com_dist_vs_cis": round(float(np.mean(com_dists)), 3),
            "com_dist_vs_cis_sd": round(float(np.std(com_dists)), 3),
            "com_dist_vs_cis_model0": com_dists[0],
            "docked_rmsd_vs_cis_5": docked_rmsds,
            "docked_rmsd_vs_cis": round(float(np.mean(docked_rmsds)), 3),
            "docked_rmsd_vs_cis_sd": round(float(np.std(docked_rmsds)), 3),
            "docked_rmsd_vs_cis_model0": docked_rmsds[0],
            "internal_fold_rmsd_5": fold_rmsds,
            "internal_fold_rmsd_mean": round(float(np.mean(fold_rmsds)), 3),
            "dna_min_dist_ca_5": dna_ca_mins,
            "dna_min_dist": round(float(np.mean(dna_ca_mins)), 3),
            "dna_min_dist_sd": round(float(np.std(dna_ca_mins)), 3),
            "dna_min_dist_heavy_5": dna_hv_mins,
            "dna_min_dist_heavy_mean": round(float(np.mean(dna_hv_mins)), 3),
            "dna_min_dist_heavy_sd": round(float(np.std(dna_hv_mins)), 3),
            "sah_min_dist_ca_mean": round(float(np.mean(sah_ca_mins)), 3),
            "cys_min_dist_ca_mean": round(float(np.mean(cys_ca_mins)), 3) if cys_ca_mins else None,
            "heavy_contacts_4_5A_5": hv_contacts_list,
            "heavy_contacts_4_5A": int(round(float(np.mean(hv_contacts_list)))),
            "ca_contacts_8A_mean": round(float(np.mean(ca_contacts_list)), 1),
            "csl_pos_B": csl_pos_B,
            "kyf_pos_B": kyf_pos_B,
            "plddt_B_model0": plddt_b_m0,
        }
        trans_records.append(rec)
        print(
            f"  {name:34s} | ipTM={rec['iptm_mean']:.3f}±{rec['iptm_sd']:.3f} | "
            f"R1-F={rec['r1_kyf_mean']:5.2f}±{rec['r1_kyf_sd']:.2f} | "
            f"COM={rec['com_dist_vs_cis']:5.2f}±{rec['com_dist_vs_cis_sd']:.2f} | "
            f"DockRMSD={rec['docked_rmsd_vs_cis']:5.2f}±{rec['docked_rmsd_vs_cis_sd']:.2f} | "
            f"DNA(CA)={rec['dna_min_dist']:5.2f}±{rec['dna_min_dist_sd']:.2f} "
            f"DNA(hv)={rec['dna_min_dist_heavy_mean']:5.2f} | A_RMSD={rec['chainA_rigid_rmsd_mean']:.2f}"
        )

    out = {
        "superposition_method": {
            "description": "Kabsch SVD over strictly rigid C5-MTase Rossmann core anchors, excluding the mobile Motif IV base-flipping loop M4[+2..+8]",
            "anchors": [{"motif": k, "rel_start": lo, "rel_end": hi} for k, lo, hi in RIGID_ANCHORS],
            "mht1_to_cis_rmsd_A": round(rmsd_mht, 3),
            "mht1_to_cis_n_ca": int(len(P_mht)),
            "mht1_active_cys_resnum": int(mht_cys_resnum),
        },
        "cis_reference": {
            "name": "AF-A0A6N9NLZ4-F1",
            "len": len(a0_seq),
            "active_cys_resnum": int(cis_cys_resnum),
            "mtase_domain_range": [1, 365],
            "cterm_core_range": [int(a0_ids[cis_core_lo]), int(a0_ids[cis_core_hi - 1])],
            "cterm_full_range": [366, 485],
            "r1_kyf_dist": round(cis_r1_y, 3),
            "r2_kyf_dist": round(cis_r2_y, 3),
            "ca_contacts_8A_core": ca_c_core,
            "ca_contacts_8A_full": ca_c_full,
            "heavy_contacts_4_5A_core": hv_c_core,
            "heavy_contacts_4_5A_full": hv_c_full,
            "dna_min_dist_ca": round(cis_dna_ca_min, 3),
            "dna_min_dist_heavy": round(cis_dna_hv_min, 3),
            "sah_min_dist_ca": round(cis_sah_ca_min, 3),
            "sah_min_dist_heavy": round(cis_sah_hv_min, 3),
            "cys_min_dist_ca": round(cis_cys_ca_min, 3),
        },
        "trans_complexes": trans_records,
    }
    out_path = sdir / "trans_vs_cis_3d_audit.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nExported complete 30-model rigid-core 3D audit to {out_path}")


if __name__ == "__main__":
    main()
