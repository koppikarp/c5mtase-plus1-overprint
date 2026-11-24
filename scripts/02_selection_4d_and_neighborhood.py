#!/usr/bin/env python3
"""
Step 2: Unified, leak-free evolutionary selection & stop-codon depletion audit
for the pan-prokaryotic C5-MTase +1 dual-coding system.

Runs 100% self-contained from repository data files:
  - data/gtdb2339_c5mtase_plus1_census.json
  - data/gtdb2339_ctrl_c5mtases.json

Computes and exports data/selection_4d_and_stop_stats.json with ZERO hardcoded P-values:
  1. Unified whole-gene A/G/P/T 4D-site expected identity (EXP_4D) + consistent 35-80%
     Frame-0 AA-identity filter applied identically to BOTH sliding-window and paired tests.
  2. Primary non-pseudoreplicated gene-level paired Wilcoxon signed-rank tests
     (inside [-30..+120] vs conserved C-term 75 aa) + genome-level cluster-bootstrap 95% CIs.
  3. Non-circular 5-frame stop-codon depletion on homology-selected low-GC (GC3_in <= 0.38)
     Arch 1 loci (ZERO stop-count filter!):
       (a) Paired within-gene Wilcoxon signed-rank (+1 vs +2, -0, -1, -2)
       (b) 1:1 GC3-matched homology-negative controls (ZERO stop-count filter!)
  4. Non-pseudoreplicated +1 amino-acid preservation at synonymous Frame-0 4D substitutions
     (gene-level Mann-Whitney U + disjoint-pairs check) + locked-3rd-base MSA audit.
"""
from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy import stats

import lib_c5mtase as lib


def _align_4d_region(
    h1: str, n1: str, h2: str, n2: str, exp_global: float
) -> Tuple[float, float, float, int, int, int]:
    """
    Align two Frame-0 segments via NW and compute:
      (aa_id, obs_4d, exc_4d, n4d, p1_pres_same, p1_pres_tot)
    using the unified whole-gene exp_global baseline (A/G/P/T only).
    """
    a1, a2 = lib.nw_align(h1, h2)
    k1 = k2 = 0
    aln_pos = id_aa = 0
    n4d = same4d = 0
    aligned_codons = []
    for c1, c2 in zip(a1, a2):
        if c1 != "-" and c2 != "-":
            aln_pos += 1
            cd1 = n1[k1 * 3:k1 * 3 + 3]
            cd2 = n2[k2 * 3:k2 * 3 + 3]
            aligned_codons.append((c1, c2, cd1, cd2))
            if c1 == c2:
                id_aa += 1
                if c1 in lib.FOURFOLD_NOSTOP:
                    if len(cd1) == 3 and len(cd2) == 3 and "N" not in cd1 and "N" not in cd2:
                        n4d += 1
                        if cd1[2] == cd2[2]:
                            same4d += 1
            k1 += 1
            k2 += 1
        else:
            aligned_codons.append(None)
            if c1 != "-":
                k1 += 1
            if c2 != "-":
                k2 += 1
    # Test 4: At contiguous aligned codons where Host Frame-0 AA DIFFERS (h1a != h2a),
    # how often does the overlapping +1 amino acid stay IDENTICAL?
    p1_same = p1_tot = 0
    for u in range(len(aligned_codons) - 1):
        cur = aligned_codons[u]
        nxt = aligned_codons[u + 1]
        if cur is None or nxt is None:
            continue
        h1a, h2a, c1a, c2a = cur
        _h1b, _h2b, c1b, c2b = nxt
        if len(c1a) < 3 or len(c2a) < 3 or len(c1b) < 3 or len(c2b) < 3:
            continue
        if h1a != h2a:
            aa_p1_1 = lib.CODON_TABLE.get(c1a[1:3] + c1b[0], "X")
            aa_p1_2 = lib.CODON_TABLE.get(c2a[1:3] + c2b[0], "X")
            if aa_p1_1 != "X" and aa_p1_2 != "X":
                p1_tot += 1
                if aa_p1_1 == aa_p1_2:
                    p1_same += 1
    if aln_pos == 0 or n4d < 5:
        return (0.0, 0.0, 0.0, 0, 0, 0)
    aa_id = id_aa / aln_pos
    obs = same4d / n4d
    return (aa_id, obs, obs - exp_global, n4d, p1_same, p1_tot)


def analyze_subgroup_4d(grp: List[dict], label: str, n_boot: int = 2000, seed: int = 42) -> dict:
    """
    Unified 4D-site analysis on a subgroup of independent C5-MTase loci:
      - Whole-gene A/G/P/T EXP_4D baseline used identically in paired & sliding-window tests
      - Strict 35-80% inside & 30-80% outside Frame-0 AA-ID filter applied to BOTH
      - Gene-level collapsed paired Wilcoxon (primary, independent per genome)
      - Genome-level cluster bootstrap 95% CIs on every 25-aa sliding-window bin
    """
    n_loci = len(grp)
    pair_records = []
    # Store per-pair bin counts for genome-level cluster bootstrap!
    pair_bins: Dict[Tuple[int, int], Dict[int, Tuple[int, int, float]]] = {}
    gene_in_vals: Dict[int, List[float]] = defaultdict(list)
    gene_out_vals: Dict[int, List[float]] = defaultdict(list)
    gene_in_obs: Dict[int, List[float]] = defaultdict(list)
    gene_out_obs: Dict[int, List[float]] = defaultdict(list)
    gene_exp: Dict[int, List[float]] = defaultdict(list)
    gene_p1_pres: Dict[int, List[float]] = defaultdict(list)

    for i in range(n_loci):
        for j in range(i + 1, n_loci):
            if grp[i]["acc"] == grp[j]["acc"]:
                continue
            h1, n1, p1 = grp[i]["hprot"], grp[i]["hnt"], grp[i]["pcq"]
            h2, n2, p2 = grp[j]["hprot"], grp[j]["hnt"], grp[j]["pcq"]
            exp_g = lib.whole_gene_4d_exp(h1, n1, h2, n2)

            lo1, hi1 = max(0, p1 - 30), min(len(h1), p1 + 120)
            lo2, hi2 = max(0, p2 - 30), min(len(h2), p2 + 120)
            aa_in, obs_in, exc_in, n4d_in, p1s, p1t = _align_4d_region(
                h1[lo1:hi1], n1[lo1 * 3:hi1 * 3], h2[lo2:hi2], n2[lo2 * 3:hi2 * 3], exp_g
            )
            aa_out, obs_out, exc_out, n4d_out, _, _ = _align_4d_region(
                h1[-75:], n1[(len(h1) - 75) * 3:], h2[-75:], n2[(len(h2) - 75) * 3:], exp_g
            )
            if not (n4d_in >= 5 and n4d_out >= 5 and 0.35 <= aa_in <= 0.80 and 0.30 <= aa_out <= 0.80):
                continue

            p1_frac = (p1s / p1t) if p1t >= 5 else None
            pair_records.append({
                "i": i, "j": j, "acc_i": grp[i]["acc"], "acc_j": grp[j]["acc"],
                "aa_id_in": round(aa_in, 4), "aa_id_out": round(aa_out, 4),
                "exp_4d": round(exp_g, 4),
                "obs_in": round(obs_in, 4), "exc_in": round(exc_in, 4),
                "obs_out": round(obs_out, 4), "exc_out": round(exc_out, 4),
                "p1_pres_same": p1s, "p1_pres_tot": p1t,
                "p1_pres_frac": round(p1_frac, 4) if p1_frac is not None else None,
            })
            for idx in (i, j):
                gene_in_vals[idx].append(exc_in)
                gene_out_vals[idx].append(exc_out)
                gene_in_obs[idx].append(obs_in)
                gene_out_obs[idx].append(obs_out)
                gene_exp[idx].append(exp_g)
                if p1_frac is not None:
                    gene_p1_pres[idx].append(p1_frac)

            # Full-length NW alignment for 25-aa sliding window (EXACT same filtered pairs & exp_g!)
            a1, a2 = lib.nw_align(h1, h2)
            k1 = k2 = 0
            b_tmp = defaultdict(lambda: [0, 0, 0.0])
            for c1, c2 in zip(a1, a2):
                if c1 != "-" and c2 != "-":
                    rel = k1 - p1
                    bi = (rel // 25) * 25
                    if -50 <= bi <= 200 and c1 == c2 and c1 in lib.FOURFOLD_NOSTOP:
                        cd1 = n1[k1 * 3:k1 * 3 + 3]
                        cd2 = n2[k2 * 3:k2 * 3 + 3]
                        if len(cd1) == 3 and len(cd2) == 3 and "N" not in cd1 and "N" not in cd2:
                            b_tmp[bi][1] += 1
                            b_tmp[bi][0] += int(cd1[2] == cd2[2])
                            b_tmp[bi][2] += exp_g
                    k1 += 1
                    k2 += 1
                else:
                    if c1 != "-":
                        k1 += 1
                    if c2 != "-":
                        k2 += 1
            pair_bins[(i, j)] = {bi: tuple(v) for bi, v in b_tmp.items()}

    # Gene-level collapsed vectors (strictly 1 observation per independent locus/genome!)
    active_genes = sorted(gene_in_vals.keys())
    g_in = [float(np.mean(gene_in_vals[idx])) for idx in active_genes]
    g_out = [float(np.mean(gene_out_vals[idx])) for idx in active_genes]
    g_in_obs = [float(np.mean(gene_in_obs[idx])) for idx in active_genes]
    g_out_obs = [float(np.mean(gene_out_obs[idx])) for idx in active_genes]
    g_exp = [float(np.mean(gene_exp[idx])) for idx in active_genes]
    g_accs = [grp[idx]["acc"] for idx in active_genes]

    if len(g_in) >= 5:
        w_stat, w_p = stats.wilcoxon(g_in, g_out, alternative="greater")
        w_stat_2s, w_p_2s = stats.wilcoxon(g_in, g_out, alternative="two-sided")
    else:
        w_stat = w_p = w_stat_2s = w_p_2s = float("nan")

    # Pair-level summary (for supplementary reference)
    p_in_exc = [r["exc_in"] for r in pair_records]
    p_out_exc = [r["exc_out"] for r in pair_records]
    p_in_obs = [r["obs_in"] for r in pair_records]
    p_out_obs = [r["obs_out"] for r in pair_records]
    p_exp = [r["exp_4d"] for r in pair_records]

    # Genome-level cluster bootstrap 95% CIs on sliding window bins!
    rng = np.random.default_rng(seed)
    boot_samples = defaultdict(list)
    if len(active_genes) >= 3:
        ag_arr = np.array(active_genes)
        for _ in range(n_boot):
            draw = rng.choice(ag_arr, size=len(ag_arr), replace=True)
            b_acc = defaultdict(lambda: [0, 0, 0.0])
            for a_pos in range(len(draw)):
                ia = int(draw[a_pos])
                for b_pos in range(a_pos + 1, len(draw)):
                    ib = int(draw[b_pos])
                    if ia == ib:
                        continue
                    key = (ia, ib) if ia < ib else (ib, ia)
                    pb = pair_bins.get(key)
                    if pb:
                        for bi, (s, t, se) in pb.items():
                            b_acc[bi][0] += s
                            b_acc[bi][1] += t
                            b_acc[bi][2] += se
            for bi in range(-50, 201, 25):
                s, t, se = b_acc[bi]
                if t >= 10:
                    boot_samples[bi].append((s / t) - (se / t))

    # Point estimates for sliding-window bins
    agg_bins = defaultdict(lambda: [0, 0, 0.0])
    for pb in pair_bins.values():
        for bi, (s, t, se) in pb.items():
            agg_bins[bi][0] += s
            agg_bins[bi][1] += t
            agg_bins[bi][2] += se

    win_rows = []
    for bi in range(-50, 201, 25):
        s, t, se = agg_bins[bi]
        if t < 15:
            continue
        obs = s / t
        exp = se / t
        exc = obs - exp
        bs = boot_samples.get(bi, [])
        if len(bs) >= 50:
            ci_lo = float(np.percentile(bs, 2.5))
            ci_hi = float(np.percentile(bs, 97.5))
            se_boot = float(np.std(bs))
        else:
            se_boot = math.sqrt(obs * (1.0 - obs) / t)
            ci_lo = exc - 1.96 * se_boot
            ci_hi = exc + 1.96 * se_boot
        win_rows.append({
            "bin_start": bi, "bin_end": bi + 24, "bin_mid": bi + 12.0,
            "n_sites": int(t), "obs_4d": round(obs, 4), "exp_4d": round(exp, 4),
            "excess_4d": round(exc, 4),
            "ci95_lo": round(ci_lo, 4), "ci95_hi": round(ci_hi, 4),
            "se_boot": round(se_boot, 4),
        })

    # Gene-level & maximum-cardinality disjoint-pair +1 AA preservation
    g_p1_vals = [float(np.mean(gene_p1_pres[idx])) for idx in sorted(gene_p1_pres.keys())]
    used_disjoint = set()
    disjoint_p1 = []
    for r in sorted(pair_records, key=lambda x: abs(x["aa_id_in"] - 0.55)):
        if r["p1_pres_frac"] is not None and r["i"] not in used_disjoint and r["j"] not in used_disjoint:
            used_disjoint.add(r["i"])
            used_disjoint.add(r["j"])
            disjoint_p1.append(r["p1_pres_frac"])

    print(
        f"  [{label:18s}] n_loci={n_loci:2d} (active={len(active_genes):2d}, pairs={len(pair_records):3d}) | "
        f"gene-level med EXC: IN={np.median(g_in):+.3f} vs OUT={np.median(g_out):+.3f} "
        f"({sum(a > b for a, b in zip(g_in, g_out))}/{len(g_in)} IN>OUT, Wilcoxon P={w_p:.2e}) | "
        f"gene-level +1 AA pres med={np.median(g_p1_vals) if g_p1_vals else float('nan'):.3f} (n={len(g_p1_vals)})"
    )
    return {
        "n_loci_input": n_loci,
        "n_genes": len(active_genes),
        "n_pairs": len(pair_records),
        "gene_accs": g_accs,
        "gene_in_excess": [round(x, 4) for x in g_in],
        "gene_out_excess": [round(x, 4) for x in g_out],
        "gene_in_obs": [round(x, 4) for x in g_in_obs],
        "gene_out_obs": [round(x, 4) for x in g_out_obs],
        "gene_exp_4d": [round(x, 4) for x in g_exp],
        "gene_in_gt_out_count": int(sum(a > b for a, b in zip(g_in, g_out))),
        "gene_med_in_excess": round(float(np.median(g_in)), 4) if g_in else None,
        "gene_med_out_excess": round(float(np.median(g_out)), 4) if g_out else None,
        "gene_med_in_obs": round(float(np.median(g_in_obs)), 4) if g_in_obs else None,
        "gene_med_out_obs": round(float(np.median(g_out_obs)), 4) if g_out_obs else None,
        "gene_med_exp_4d": round(float(np.median(g_exp)), 4) if g_exp else None,
        "gene_wilcoxon_p_greater": float(w_p) if not math.isnan(w_p) else None,
        "gene_wilcoxon_p_twosided": float(w_p_2s) if not math.isnan(w_p_2s) else None,
        "pair_in_excess": p_in_exc,
        "pair_out_excess": p_out_exc,
        "pair_in_obs": p_in_obs,
        "pair_out_obs": p_out_obs,
        "pair_exp_4d": p_exp,
        "pair_med_in_excess": round(float(np.median(p_in_exc)), 4) if p_in_exc else None,
        "pair_med_out_excess": round(float(np.median(p_out_exc)), 4) if p_out_exc else None,
        "gene_p1_aa_preservation": [round(x, 4) for x in g_p1_vals],
        "disjoint_pair_p1_aa_preservation": [round(x, 4) for x in disjoint_p1],
        "windows": win_rows,
    }


def compute_locked_3rd_bases(fpcq_low_cases: List[dict]) -> dict:
    """
    Exact site-by-site locked 3rd-base audit across all low-GC (GC3_in <= 0.38)
    homology-verified GTDB Arch 1 FPCQ loci (ZERO false-positive high-GC noise!).
    """
    ref = [c for c in fpcq_low_cases if c["acc"] == "GCA_002355375.1"][0]
    ref_h, ref_nt, ref_pcq = ref["hprot"], ref["hnt"], ref["pcq"]
    lo_r, hi_r = max(0, ref_pcq - 45), min(len(ref_h), ref_pcq + 130)
    L_win = hi_r - lo_r

    h_rows, p1_rows, b3_rows = [], [], []
    for c in fpcq_low_cases:
        hprot, hnt, pcq = c["hprot"], c["hnt"], c["pcq"]
        p1_full = lib.translate(hnt, offset=1)
        q_off = max(0, pcq - 55)
        a_ref, a_q = lib.nw_align(ref_h[lo_r:hi_r], hprot[q_off:min(len(hprot), pcq + 140)])
        h_col, p1_col, b3_col = [], [], []
        qi = 0
        for cr, cq in zip(a_ref, a_q):
            if cr != "-":
                if cq != "-":
                    pos_q = q_off + qi
                    h_col.append(hprot[pos_q])
                    b3_col.append(hnt[pos_q * 3 + 2] if pos_q * 3 + 2 < len(hnt) else "-")
                    p1_col.append(p1_full[pos_q] if pos_q < len(p1_full) else "-")
                    qi += 1
                else:
                    h_col.append("-")
                    b3_col.append("-")
                    p1_col.append("-")
            else:
                if cq != "-":
                    qi += 1
        h_rows.append("".join(h_col))
        p1_rows.append("".join(p1_col))
        b3_rows.append("".join(b3_col))

    locked = []
    n_seq = len(fpcq_low_cases)
    for k in range(L_win):
        b3s = [r[k] for r in b3_rows if r[k] in "ACGT"]
        hs = [r[k] for r in h_rows if r[k] != "-"]
        p1s = [r[k] for r in p1_rows if r[k] != "-"]
        if len(b3s) < n_seq - 1:
            continue
        tb, cb = Counter(b3s).most_common(1)[0]
        th, ch = Counter(hs).most_common(1)[0]
        tp, cp = Counter(p1s).most_common(1)[0]
        if cb / len(b3s) >= 0.90:
            locked.append({
                "rel_pcq": k - 45,
                "pos3_base": tb,
                "pos3_count": cb,
                "pos3_total": len(b3s),
                "host_aa_top": th,
                "host_aa_counts": dict(Counter(hs)),
                "plus1_aa_top": tp,
                "plus1_aa_counts": dict(Counter(p1s)),
                "is_4d_host": th in lib.FOURFOLD_NOSTOP or th == "V",
            })
    return {
        "n_loci": n_seq,
        "loci_accs": [c["acc"] for c in fpcq_low_cases],
        "reference_acc": ref["acc"],
        "n_sites_ge90pct": len(locked),
        "n_sites_100pct": sum(1 for x in locked if x["pos3_count"] == x["pos3_total"]),
        "locked_sites": locked,
    }


def main():
    gtdb = json.loads((lib.DATA_DIR / "gtdb2339_c5mtase_plus1_census.json").read_text())
    ctrl_data = json.loads((lib.DATA_DIR / "gtdb2339_ctrl_c5mtases.json").read_text())
    ov = gtdb["overprinted"]
    ctrls_all = ctrl_data["controls"]

    # Subgroups for Test 1 & Test 2 (4D-site synonymous conservation):
    # Strictly intact Arch 1 loci (1 per unique genome where possible; GC3_in filtered for unsaturation)
    fpcq_low = [
        c for c in ov
        if c["motif4"] == "FPCQ" and c["gc3_in"] <= 0.45 and c.get("refined_arch") == "Arch1_intact"
    ]
    tspcq_all = [
        c for c in ov
        if c["motif4"] in ("TPCQ", "SPCQ") and c["gc3_in"] <= 0.60 and c.get("refined_arch") == "Arch1_intact"
    ]
    # Deduplicate by acc (keep first) for strict genome independence
    def uniq_acc(lst):
        seen, out = set(), []
        for x in lst:
            if x["acc"] not in seen:
                seen.add(x["acc"])
                out.append(x)
        return out

    fpcq_low = uniq_acc(fpcq_low)
    tspcq_all = uniq_acc(tspcq_all)
    ctrl_4d = uniq_acc([c for c in ctrls_all if c["gc3_in"] <= 0.45])

    print("=" * 105)
    print("TEST 1 & 2: UNIFIED WHOLE-GENE EXP_4D & GENE-LEVEL PAIRED WILCOXON (35-80% AA-ID)")
    print("=" * 105)
    sub_fpcq = analyze_subgroup_4d(fpcq_low, "Arch1 FPCQ")
    sub_tspcq = analyze_subgroup_4d(tspcq_all, "Arch1 TPCQ/SPCQ")
    sub_ctrl = analyze_subgroup_4d(ctrl_4d, "Ctrl +1-less")

    # Between-group gene-level Mann-Whitney U (Arch1 FPCQ vs Ctrl +1-less)
    _, fpcq_vs_ctrl_in_p = stats.mannwhitneyu(
        sub_fpcq["gene_in_excess"], sub_ctrl["gene_in_excess"], alternative="greater"
    )
    fpcq_delta = [a - b for a, b in zip(sub_fpcq["gene_in_excess"], sub_fpcq["gene_out_excess"])]
    ctrl_delta = [a - b for a, b in zip(sub_ctrl["gene_in_excess"], sub_ctrl["gene_out_excess"])]
    _, fpcq_vs_ctrl_delta_p = stats.mannwhitneyu(fpcq_delta, ctrl_delta, alternative="greater")
    tspcq_delta = [a - b for a, b in zip(sub_tspcq["gene_in_excess"], sub_tspcq["gene_out_excess"])]
    _, tspcq_vs_ctrl_in_p = stats.mannwhitneyu(
        sub_tspcq["gene_in_excess"], sub_ctrl["gene_in_excess"], alternative="greater"
    )
    _, tspcq_vs_ctrl_delta_p = stats.mannwhitneyu(tspcq_delta, ctrl_delta, alternative="greater")
    print(f"  Between-group gene-level MWU (Arch1 FPCQ vs Ctrl): IN excess P={fpcq_vs_ctrl_in_p:.2e} | (IN-OUT) delta P={fpcq_vs_ctrl_delta_p:.2e}")
    print(f"  Between-group gene-level MWU (Arch1 T/SPCQ vs Ctrl): IN excess P={tspcq_vs_ctrl_in_p:.2e} | (IN-OUT) delta P={tspcq_vs_ctrl_delta_p:.2e}")

    # ------------------------------------------------------------------
    # TEST 3: NON-CIRCULAR 5-FRAME STOP-CODON DEPLETION (ZERO STOP FILTER!)
    # ------------------------------------------------------------------
    print("\n" + "=" * 105)
    print("TEST 3: NON-CIRCULAR 5-FRAME STOP-CODON DEPLETION (HOMOLOGY-SELECTED, GC3_in <= 0.38)")
    print("=" * 105)
    # Select Arch 1 cases strictly by homology + low GC3_in <= 0.38 (ZERO stop-count filter!)
    low_cases = uniq_acc([
        c for c in ov
        if c.get("refined_arch") == "Arch1_intact" and c["gc3_in"] <= 0.38
    ])
    # Match 1:1 to nearest-GC3 homology-negative controls (ZERO stop-count filter!)
    pool_ctrls = [c for c in ctrls_all if c["gc3_in"] <= 0.45]
    avail = list(range(len(pool_ctrls)))
    matched_ctrls = []
    for cs in sorted(low_cases, key=lambda x: x["gc3_in"]):
        best_idx = min(avail, key=lambda k: abs(pool_ctrls[k]["gc3_in"] - cs["gc3_in"]))
        avail.remove(best_idx)
        matched_ctrls.append(pool_ctrls[best_idx])

    gc3_cases = [c["gc3_in"] for c in low_cases]
    gc3_mctrls = [c["gc3_in"] for c in matched_ctrls]
    _, gc3_mwu_p = stats.mannwhitneyu(gc3_cases, gc3_mctrls, alternative="two-sided")
    print(f"  Cases (homology-selected Arch1, GC3<=0.38): n={len(low_cases)}, med GC3={np.median(gc3_cases):.3f} [{min(gc3_cases):.3f}-{max(gc3_cases):.3f}]")
    print(f"  1:1 GC3-matched homology-negative controls: n={len(matched_ctrls)}, med GC3={np.median(gc3_mctrls):.3f} [{min(gc3_mctrls):.3f}-{max(gc3_mctrls):.3f}] (GC3 MWU P={gc3_mwu_p:.3f})")

    frames = ("+1", "+2", "-0", "-1", "-2")
    case_stops = {fr: [] for fr in frames}
    ctrl_stops = {fr: [] for fr in frames}
    for c in low_cases:
        d = lib.count_5frames_stops(c["hnt"], c["pcq"], rel_lo=-25, rel_hi=110)
        for fr in frames:
            case_stops[fr].append(d[fr])
    for c in matched_ctrls:
        tbl = lib.CODON_TABLE.copy()
        if c.get("use_t4"):
            tbl["TGA"] = "W"
        d = lib.count_5frames_stops(c["hnt"], c["pcq"], rel_lo=-25, rel_hi=110, table=tbl)
        for fr in frames:
            ctrl_stops[fr].append(d[fr])

    within_gene_wilcoxon = {}
    case_vs_ctrl_mwu = {}
    for fr in frames:
        _, p_cc = stats.mannwhitneyu(case_stops[fr], ctrl_stops[fr], alternative="two-sided")
        _, p_cc_less = stats.mannwhitneyu(case_stops[fr], ctrl_stops[fr], alternative="less")
        case_vs_ctrl_mwu[fr] = {
            "case_median": float(np.median(case_stops[fr])),
            "case_mean": round(float(np.mean(case_stops[fr])), 2),
            "ctrl_median": float(np.median(ctrl_stops[fr])),
            "ctrl_mean": round(float(np.mean(ctrl_stops[fr])), 2),
            "mwu_p_twosided": float(p_cc),
            "mwu_p_less": float(p_cc_less),
        }
        if fr != "+1":
            _, p_wg = stats.wilcoxon(case_stops["+1"], case_stops[fr], alternative="less")
            _, p_wg_2s = stats.wilcoxon(case_stops["+1"], case_stops[fr], alternative="two-sided")
            within_gene_wilcoxon[fr] = {
                "p_less": float(p_wg),
                "p_twosided": float(p_wg_2s),
                "n_pairs": len(low_cases),
                "n_plus1_lt_other": int(sum(a < b for a, b in zip(case_stops["+1"], case_stops[fr]))),
            }
        print(
            f"    Frame {fr:2s}: CASE med={np.median(case_stops[fr]):4.1f} (mean={np.mean(case_stops[fr]):4.1f}) | "
            f"CTRL med={np.median(ctrl_stops[fr]):4.1f} (mean={np.mean(ctrl_stops[fr]):4.1f}) | "
            f"Case-vs-Ctrl MWU P(2s)={p_cc:.2e}"
            + (f" | Within-gene (+1 < {fr}) Wilcoxon P={within_gene_wilcoxon[fr]['p_less']:.2e}" if fr != "+1" else "")
        )

    # ------------------------------------------------------------------
    # TEST 4: +1 AA PRESERVATION (GENE-LEVEL & DISJOINT PAIRS MWU!) + LOCKED 3RD BASES
    # ------------------------------------------------------------------
    print("\n" + "=" * 105)
    print("TEST 4: +1 AA PRESERVATION WHEN HOST FRAME-0 AA CHANGES (NON-PSEUDOREPLICATED MWU)")
    print("=" * 105)
    g_p1_case = sub_fpcq["gene_p1_aa_preservation"]
    g_p1_ctrl = sub_ctrl["gene_p1_aa_preservation"]
    _, p1_mwu_gene_gt = stats.mannwhitneyu(g_p1_case, g_p1_ctrl, alternative="greater")
    _, p1_mwu_gene_2s = stats.mannwhitneyu(g_p1_case, g_p1_ctrl, alternative="two-sided")

    dj_case = sub_fpcq["disjoint_pair_p1_aa_preservation"]
    dj_ctrl = sub_ctrl["disjoint_pair_p1_aa_preservation"]
    _, p1_mwu_dj_gt = stats.mannwhitneyu(dj_case, dj_ctrl, alternative="greater")
    _, p1_mwu_dj_2s = stats.mannwhitneyu(dj_case, dj_ctrl, alternative="two-sided")

    print(
        f"  Gene-level MWU: CASE med={np.median(g_p1_case)*100:.1f}% (mean={np.mean(g_p1_case)*100:.1f}%, n={len(g_p1_case)}) "
        f"vs CTRL med={np.median(g_p1_ctrl)*100:.1f}% (mean={np.mean(g_p1_ctrl)*100:.1f}%, n={len(g_p1_ctrl)}) | "
        f"MWU P(greater)={p1_mwu_gene_gt:.2e}, P(2s)={p1_mwu_gene_2s:.2e}"
    )
    print(
        f"  Disjoint-pair MWU: CASE med={np.median(dj_case)*100:.1f}% (n={len(dj_case)}) "
        f"vs CTRL med={np.median(dj_ctrl)*100:.1f}% (n={len(dj_ctrl)}) | "
        f"MWU P(greater)={p1_mwu_dj_gt:.2e}, P(2s)={p1_mwu_dj_2s:.2e}"
    )

    fpcq_low38 = [c for c in low_cases if c["motif4"] == "FPCQ"]
    locked_audit = compute_locked_3rd_bases(fpcq_low38)
    print(f"\n  Locked 3rd-base audit across {locked_audit['n_loci']} homology-verified low-GC (<=0.38) FPCQ Arch 1 loci:")
    for s in locked_audit["locked_sites"]:
        if 60 <= s["rel_pcq"] <= 115:
            print(
                f"    rel={s['rel_pcq']:+4d} | pos3={s['pos3_base']} ({s['pos3_count']}/{s['pos3_total']}) | "
                f"H0={s['host_aa_top']} {s['host_aa_counts']} | +1={s['plus1_aa_top']} {s['plus1_aa_counts']}"
            )

    out = {
        "methodology": {
            "fourfold_amino_acids": sorted(lib.FOURFOLD_NOSTOP),
            "exp_4d_definition": "Whole-gene A/G/P/T third-position base-composition baseline applied identically to paired and sliding-window comparisons",
            "aa_id_filter": {"inside_min": 0.35, "inside_max": 0.80, "outside_min": 0.30, "outside_max": 0.80},
            "stop_window_rel_pcq": [-25, 110],
            "stop_selection_rule": "Strictly homology-selected Arch1_intact loci (GC3_in <= 0.38) and 1:1 nearest-GC3 homology-negative controls; ZERO stop-count filter applied.",
            "bootstrap_replicates": 2000,
        },
        "fpcq_overprinted": sub_fpcq,
        "tspcq_overprinted": sub_tspcq,
        "ctrl_no_plus1": sub_ctrl,
        "stop_purge_5frames": {
            "window_rel_pcq": [-25, 110],
            "n_cases": len(low_cases),
            "n_ctrls": len(matched_ctrls),
            "case_accs": [c["acc"] for c in low_cases],
            "ctrl_accs": [c["acc"] for c in matched_ctrls],
            "case_gc3_median": round(float(np.median(gc3_cases)), 4),
            "ctrl_gc3_median": round(float(np.median(gc3_mctrls)), 4),
            "gc3_balance_mwu_p": round(float(gc3_mwu_p), 4),
            "cases": case_stops,
            "ctrls": ctrl_stops,
            "within_gene_wilcoxon_plus1_lt_other": within_gene_wilcoxon,
            "case_vs_ctrl_mwu": case_vs_ctrl_mwu,
        },
        "plus1_aa_preservation": {
            "test_name": "Mann-Whitney U test (gene-level independent means + disjoint non-overlapping pairs sensitivity check)",
            "gene_level": {
                "n_case_genes": len(g_p1_case),
                "n_ctrl_genes": len(g_p1_ctrl),
                "case_median": round(float(np.median(g_p1_case)), 4),
                "case_mean": round(float(np.mean(g_p1_case)), 4),
                "ctrl_median": round(float(np.median(g_p1_ctrl)), 4),
                "ctrl_mean": round(float(np.mean(g_p1_ctrl)), 4),
                "mwu_p_greater": float(p1_mwu_gene_gt),
                "mwu_p_twosided": float(p1_mwu_gene_2s),
            },
            "disjoint_pairs": {
                "n_case_pairs": len(dj_case),
                "n_ctrl_pairs": len(dj_ctrl),
                "case_median": round(float(np.median(dj_case)), 4),
                "ctrl_median": round(float(np.median(dj_ctrl)), 4),
                "mwu_p_greater": float(p1_mwu_dj_gt),
                "mwu_p_twosided": float(p1_mwu_dj_2s),
            },
        },
        "locked_3rd_bases_lowgc_fpcq": locked_audit,
        "p_values": {
            "fpcq_gene_wilcoxon_in_vs_out": sub_fpcq["gene_wilcoxon_p_greater"],
            "tspcq_gene_wilcoxon_in_vs_out": sub_tspcq["gene_wilcoxon_p_greater"],
            "ctrl_gene_wilcoxon_in_vs_out": sub_ctrl["gene_wilcoxon_p_greater"],
            "fpcq_vs_ctrl_in_excess_mwu_p": float(fpcq_vs_ctrl_in_p),
            "fpcq_vs_ctrl_delta_mwu_p": float(fpcq_vs_ctrl_delta_p),
            "tspcq_vs_ctrl_in_excess_mwu_p": float(tspcq_vs_ctrl_in_p),
            "tspcq_vs_ctrl_delta_mwu_p": float(tspcq_vs_ctrl_delta_p),
            "stop_purge_within_gene_max_p": max(v["p_less"] for v in within_gene_wilcoxon.values()),
            "stop_purge_plus1_case_vs_ctrl_mwu_p": case_vs_ctrl_mwu["+1"]["mwu_p_twosided"],
            "plus1_aa_preservation_gene_mwu_p": float(p1_mwu_gene_gt),
            "plus1_aa_preservation_gene_mwu_p_2s": float(p1_mwu_gene_2s),
            "plus1_aa_preservation_disjoint_mwu_p": float(p1_mwu_dj_gt),
            "plus1_aa_preservation_case_med": round(float(np.median(g_p1_case)), 4),
            "plus1_aa_preservation_ctrl_med": round(float(np.median(g_p1_ctrl)), 4),
        },
    }
    out_path = lib.DATA_DIR / "selection_4d_and_stop_stats.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nExported unified, zero-hardcoded selection stats to {out_path}")


if __name__ == "__main__":
    main()
