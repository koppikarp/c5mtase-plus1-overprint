#!/usr/bin/env python3
"""
Step 3B: Programmatic ML tree analysis, bootstrap support verification,
midpoint/outgroup rooting sensitivity, long-branch sensitivity check,
Fitch parsimony, and patristic congruence (full + non-overlapping Arch 3
and Arch 1-vs-Arch 3 Mantel controls).

Runs 100% self-contained from data/phylogeny/ and exports
data/phylogeny/ml_bootstrap_summary.json with ZERO hardcoded support numbers!
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import numpy as np

import lib_c5mtase as lib


def load_boot_splits(boot_path: Path) -> List[set]:
    lines = [l.strip() for l in boot_path.read_text().splitlines() if l.strip()]
    return [lib.get_unrooted_splits(lib.parse_newick(l)) for l in lines]


def find_lca_node(root: lib.Node, target_leaves: Set[str]) -> lib.Node | None:
    """Return the deepest node in a rooted tree whose leaf_names() is a superset of target_leaves."""
    best = None

    def walk(nd: lib.Node):
        nonlocal best
        lvs = nd.leaf_names()
        if target_leaves.issubset(lvs):
            if best is None or len(lvs) < len(best.leaf_names()):
                best = nd
            for c in nd.children:
                walk(c)

    walk(root)
    return best


def clade_metrics(
    rooted_tree: lib.Node,
    all_leaves: Set[str],
    boot_splits: List[set],
    target_leaves: Set[str],
) -> dict:
    """Compute exact LCA size, monophyly status, SH-like support, and ML bootstrap %."""
    lca = find_lca_node(rooted_tree, target_leaves)
    if lca is None:
        return {"found": False}
    lca_lvs = lca.leaf_names()
    sp = lib.canonical_split(lca_lvs, all_leaves)
    ml_boot = 100.0 * sum(1 for bs in boot_splits if sp in bs) / max(1, len(boot_splits))
    sh = round(lca.support * 100.0, 1) if lca.support is not None else None
    comp = Counter(lib.classify_taxon(x) for x in lca_lvs)
    return {
        "found": True,
        "monophyletic": (lca_lvs == target_leaves),
        "n_target": len(target_leaves),
        "n_lca": len(lca_lvs),
        "ml_boot": round(ml_boot, 1),
        "sh": sh,
        "composition": dict(comp),
        "branch_length": round(float(lca.length), 4),
    }


def fitch_parsimony(root: lib.Node, state_fn) -> Tuple[int, Set[str]]:
    changes = 0

    def post(nd: lib.Node) -> Set[str]:
        nonlocal changes
        if nd.is_leaf():
            return {state_fn(nd.name)}
        sets = [post(c) for c in nd.children]
        inter = set.intersection(*sets)
        if inter:
            return inter
        changes += 1
        return set.union(*sets)

    rset = post(root)
    return changes, rset


def cross_block_mantel(
    m1: np.ndarray, m2: np.ndarray, idx_a: List[int], idx_b: List[int], n_perm: int = 5000, seed: int = 42
) -> dict:
    """Permutation Mantel test strictly on cross-block pairs (i in A, j in B; zero shared DNA on Arch 3!)."""
    sub1 = m1[np.ix_(idx_a, idx_b)].ravel()
    sub2 = m2[np.ix_(idx_a, idx_b)].ravel()
    r_obs = lib.pearson_r(sub1, sub2)
    rho_obs = lib.pearson_r(lib.rank_vec(sub1), lib.rank_vec(sub2))
    rng = np.random.default_rng(seed)
    ge_r = ge_rho = 0
    rv1 = np.asarray(lib.rank_vec(sub1), dtype=float)
    a_arr = np.array(idx_a)
    b_arr = np.array(idx_b)
    for _ in range(n_perm):
        pa = rng.permutation(a_arr)
        pb = rng.permutation(b_arr)
        vp = m2[np.ix_(pa, pb)].ravel()
        if lib.pearson_r(sub1, vp) >= r_obs:
            ge_r += 1
        if lib.pearson_r(rv1, lib.rank_vec(vp)) >= rho_obs:
            ge_rho += 1
    return {
        "n_group_a": len(idx_a),
        "n_group_b": len(idx_b),
        "n_pairs": int(len(sub1)),
        "pearson_r": round(r_obs, 4),
        "spearman_rho": round(rho_obs, 4),
        "mantel_p_pearson": round((ge_r + 1) / (n_perm + 1), 6),
        "mantel_p_spearman": round((ge_rho + 1) / (n_perm + 1), 6),
        "n_perm": int(n_perm),
    }


def main():
    pdir = lib.PHYLO_DIR
    t1_raw = lib.parse_newick((pdir / "tree1_frame0_cat_ml_lg.nwk").read_text())
    t1 = lib.reroot_on_clade(t1_raw, ["OUT_Human_DNMT1", "OUT_Human_TRDMT1"])
    t1_mid = lib.midpoint_reroot(lib.parse_newick((pdir / "tree1_frame0_cat_ml_lg.nwk").read_text()))
    t2 = lib.parse_newick((pdir / "tree2_cslsdile_domain_ml_lg.nwk").read_text())
    t1_nlb_raw = lib.parse_newick((pdir / "tree1_no_longbranch_ml_lg.nwk").read_text())
    t1_nlb = lib.reroot_on_clade(t1_nlb_raw, ["OUT_Human_DNMT1", "OUT_Human_TRDMT1"])

    b1 = load_boot_splits(pdir / "tree1_boot100.nwk")
    b2 = load_boot_splits(pdir / "tree2_boot100.nwk")
    b1_nlb = load_boot_splits(pdir / "tree1_no_longbranch_boot100.nwk")

    lvs1 = t1.leaf_names()
    lvs2 = t2.leaf_names()
    lvs1_nlb = t1_nlb.leaf_names()

    # Check extreme branch lengths in Tree 1
    leaf_bls = sorted([(l.length, l.name) for l in t1.leaves()], reverse=True)
    print("Top 5 terminal branch lengths in Tree 1:", [(round(bl, 3), nm) for bl, nm in leaf_bls[:5]])

    # Exact target taxon sets
    a3_set = {x for x in lvs1 if x.startswith("A3_")}
    a2_set = {x for x in lvs1 if x.startswith("A2_")}
    a1_tpcq_set = {x for x in lvs1 if x.startswith("A1_TPCQ_")}
    a1_fpcq_set = {x for x in lvs1 if x.startswith("A1_FPCQ_")}
    arc_set = {x for x in lvs1 if "_Arc_" in x}

    # Check Arch 2 + closest Arch 1 FPCQ sister clade
    lca_a2 = find_lca_node(t1, a2_set)
    a2_parent = lca_a2.parent if lca_a2 and lca_a2.parent else lca_a2
    a2_parent_set = a2_parent.leaf_names() if a2_parent else set()

    # Check TPCQ + Arch 3 superclade (LCA of all Arch 3 + all A1_TPCQ)
    tpcq_a3_lca = find_lca_node(t1, a3_set | a1_tpcq_set)
    tpcq_a3_set = tpcq_a3_lca.leaf_names() if tpcq_a3_lca else set()

    # Check Tree 2 Arch 3 monophyly (reroot Tree 2 on Arch 1 FPCQ or unrooted split check!)
    t2_sp_a3 = lib.canonical_split(a3_set, lvs2)
    t2_a3_in_main = t2_sp_a3 in lib.get_unrooted_splits(t2)
    t2_a3_boot = round(100.0 * sum(1 for bs in b2 if t2_sp_a3 in bs) / len(b2), 1)
    # Reroot Tree 2 on an Arch 1 FPCQ clade to inspect Arch 3 & Arch 1 FPCQ SH supports
    t2_r = lib.reroot_on_clade(lib.parse_newick((pdir / "tree2_cslsdile_domain_ml_lg.nwk").read_text()),
                               [x for x in lvs2 if x.startswith("A1_FPCQ_Arc_")][:1])
    a3_t2_metrics = clade_metrics(t2_r, lvs2, b2, a3_set)
    a3_t2_metrics["unrooted_split_in_main"] = t2_a3_in_main
    a3_t2_metrics["unrooted_ml_boot"] = t2_a3_boot

    # Check Arch 1 FPCQ clade in Tree 2
    t2_r_a3 = lib.reroot_on_clade(lib.parse_newick((pdir / "tree2_cslsdile_domain_ml_lg.nwk").read_text()), a3_set)
    a1_fpcq_t2_set = {x for x in lvs2 if x.startswith("A1_FPCQ_")}
    a1_fpcq_t2_metrics = clade_metrics(t2_r_a3, lvs2, b2, a1_fpcq_t2_set)

    # Check DPANN Archaea nested inside Bacteria in Tree 1
    lca_arc = find_lca_node(t1, arc_set)
    arc_metrics = clade_metrics(t1, lvs1, b1, arc_set)

    # Sensitivity check: Tree 1 WITHOUT extreme long branch A1_SPCT_Baci_10645
    t1_nlb_a3 = clade_metrics(t1_nlb, lvs1_nlb, b1_nlb, a3_set)
    t1_nlb_a2 = clade_metrics(t1_nlb, lvs1_nlb, b1_nlb, a2_set)
    t1_nlb_tpcq_a3 = clade_metrics(t1_nlb, lvs1_nlb, b1_nlb, a3_set | a1_tpcq_set)

    # Sensitivity check: Midpoint rooting of Tree 1
    mid_a3 = clade_metrics(t1_mid, lvs1, b1, a3_set)
    mid_a2 = clade_metrics(t1_mid, lvs1, b1, a2_set)
    mid_tpcq_a3 = clade_metrics(t1_mid, lvs1, b1, a3_set | a1_tpcq_set)

    # Fitch parsimony on DNMT-rooted Tree 1 (binary +1/domain present vs absent)
    def bin_state(nm):
        c = lib.classify_taxon(nm)
        return "HAS_CSLS" if c in ("ARCH1_FPCQ", "ARCH1_TPCQ", "ARCH1_SPCQ", "ARCH2", "ARCH2_SPLIT_TRD", "ARCH3") else "NO_CSLS"

    fitch_n, fitch_root = fitch_parsimony(t1, bin_state)
    # Inspect the prokaryote ingroup child of root
    prok_node = max(t1.children, key=lambda c: len(c.leaf_names()))
    fitch_prok_n, fitch_prok_root = fitch_parsimony(prok_node, bin_state)

    # Check Bacillota SPCQ/SPCT clade containing SplitTRD_SPCQ_Baci_34735 (GCA_035934735.1)
    split_trd_set = {"SplitTRD_SPCQ_Baci_34735", "A1_SPCQ_Baci_09825"}
    split_trd_clade = clade_metrics(t1, lvs1, b1, split_trd_set)

    key_supports = {
        "Arch3_monophyly_tree1": clade_metrics(t1, lvs1, b1, a3_set),
        "Arch3_monophyly_tree2": a3_t2_metrics,
        "TPCQ_Arch3_superclade_tree1": clade_metrics(t1, lvs1, b1, tpcq_a3_set),
        "Arch2_monophyly_tree1": clade_metrics(t1, lvs1, b1, a2_set),
        "Arch2_parent_clade_tree1": clade_metrics(t1, lvs1, b1, a2_parent_set),
        "Bacillota_SPCQ_SplitTRD_clade_tree1": split_trd_clade,
        "Arch1_FPCQ_tree2": a1_fpcq_t2_metrics,
        "DPANN_Archaea_in_Bacteria_tree1": arc_metrics,
    }
    sensitivity = {
        "drop_longbranch_A1_SPCT_Baci_10645": {
            "n_taxa": len(lvs1_nlb),
            "Arch3_monophyly": t1_nlb_a3,
            "Arch2_monophyly": t1_nlb_a2,
            "TPCQ_Arch3_superclade": t1_nlb_tpcq_a3,
        },
        "midpoint_root_tree1": {
            "Arch3_monophyly": mid_a3,
            "Arch2_monophyly": mid_a2,
            "TPCQ_Arch3_superclade": mid_tpcq_a3,
        },
    }

    print("\n=== PROGRAMMATICALLY COMPUTED KEY CLADE SUPPORTS (UNIFORM THOROUGH BOOTS) ===")
    for k, v in key_supports.items():
        print(f"  {k:34s}: mono={v.get('monophyletic')} n={v.get('n_target')}/{v.get('n_lca')} ML_boot={v.get('ml_boot')}% SH={v.get('sh')}% comp={v.get('composition')}")
    print("\n=== SENSITIVITY CHECK: DROP LONG BRANCH A1_SPCT_Baci_10645 (96 TAXA) ===")
    for k, v in sensitivity["drop_longbranch_A1_SPCT_Baci_10645"].items():
        if isinstance(v, dict):
            print(f"  {k:26s}: mono={v.get('monophyletic')} n={v.get('n_target')}/{v.get('n_lca')} ML_boot={v.get('ml_boot')}% SH={v.get('sh')}%")
    print(f"Fitch parsimony (full rooted): {fitch_n} changes, root={sorted(fitch_root)} | prok ingroup: {fitch_prok_n} changes, root={sorted(fitch_prok_root)}")

    # ------------------------------------------------------------------
    # PATRISTIC CONGRUENCE & NON-OVERLAPPING MANTEL CONTROLS (M2)
    # ------------------------------------------------------------------
    n1, mat1 = lib.patristic_matrix(t1)
    n2, mat2 = lib.patristic_matrix(t2)
    idx1 = {lib.norm_id(x): i for i, x in enumerate(n1)}
    idx2 = {lib.norm_id(x): i for i, x in enumerate(n2)}
    shared = sorted(set(idx1) & set(idx2))

    m1_sh = mat1[np.ix_([idx1[s] for s in shared], [idx1[s] for s in shared])]
    m2_sh = mat2[np.ix_([idx2[s] for s in shared], [idx2[s] for s in shared])]

    mantel_all = lib.mantel_test(m1_sh, m2_sh, n_perm=5000, seed=42)

    # Subset indices in `shared`
    sh_a3 = [i for i, s in enumerate(shared) if lib.classify_taxon(s) == "ARCH3"]
    sh_a1 = [i for i, s in enumerate(shared) if lib.classify_taxon(s).startswith("ARCH1")]
    sh_a2 = [i for i, s in enumerate(shared) if lib.classify_taxon(s) == "ARCH2"]

    mantel_a3_only = lib.mantel_test(m1_sh[np.ix_(sh_a3, sh_a3)], m2_sh[np.ix_(sh_a3, sh_a3)], n_perm=5000, seed=42)
    mantel_a1_vs_a3 = cross_block_mantel(m1_sh, m2_sh, sh_a1, sh_a3, n_perm=5000, seed=42)
    mantel_a1_only = lib.mantel_test(m1_sh[np.ix_(sh_a1, sh_a1)], m2_sh[np.ix_(sh_a1, sh_a1)], n_perm=5000, seed=42)
    sh_no_lb = [i for i, s in enumerate(shared) if s != "A1_SPCT_Baci_10645"]
    mantel_no_lb = lib.mantel_test(m1_sh[np.ix_(sh_no_lb, sh_no_lb)], m2_sh[np.ix_(sh_no_lb, sh_no_lb)], n_perm=5000, seed=42)

    print("\n=== PATRISTIC CONGRUENCE (MANTEL TESTS, 5,000 PERMUTATIONS) ===")
    print(f"  All shared ({mantel_all['n_taxa']} taxa, {mantel_all['n_pairs']} pairs): r={mantel_all['pearson_r']}, rho={mantel_all['spearman_rho']}, P={mantel_all['mantel_p_pearson']}")
    print(f"  Arch 3 ONLY (ZERO overlapping DNA! {mantel_a3_only['n_taxa']} taxa, {mantel_a3_only['n_pairs']} pairs): r={mantel_a3_only['pearson_r']}, rho={mantel_a3_only['spearman_rho']}, P={mantel_a3_only['mantel_p_pearson']}")
    print(f"  Arch 1 vs Arch 3 (cross-block, {mantel_a1_vs_a3['n_pairs']} pairs): r={mantel_a1_vs_a3['pearson_r']}, rho={mantel_a1_vs_a3['spearman_rho']}, P={mantel_a1_vs_a3['mantel_p_pearson']}")
    print(f"  Arch 1 ONLY (caveat: shared DNA, {mantel_a1_only['n_taxa']} taxa, {mantel_a1_only['n_pairs']} pairs): r={mantel_a1_only['pearson_r']}, rho={mantel_a1_only['spearman_rho']}, P={mantel_a1_only['mantel_p_pearson']}")

    # Export all patristic pairs for Fig. 3c scatter plot
    pairs_data = []
    for i in range(len(shared)):
        for j in range(i + 1, len(shared)):
            a, b = shared[i], shared[j]
            c1 = lib.classify_taxon(a)
            c2 = lib.classify_taxon(b)
            if c1 == "ARCH3" and c2 == "ARCH3":
                grp = "Arch3-Arch3"
            elif "ARCH3" in (c1, c2):
                grp = "Arch1-Arch3"
            elif c1.startswith("ARCH2") or c2.startswith("ARCH2"):
                grp = "Arch2-any"
            else:
                grp = "Arch1-Arch1"
            pairs_data.append({
                "a": a, "b": b,
                "d_f0": round(float(m1_sh[i, j]), 4),
                "d_p1": round(float(m2_sh[i, j]), 4),
                "pair_type": grp,
            })

    summary = {
        "n_taxa_tree1": len(n1),
        "n_taxa_tree2": len(n2),
        "n_shared": len(shared),
        "n_pairs": len(pairs_data),
        "tree1_composition": dict(Counter(lib.classify_taxon(x) for x in n1)),
        "tree2_composition": dict(Counter(lib.classify_taxon(x) for x in n2)),
        "pearson_r": mantel_all["pearson_r"],
        "spearman_rho": mantel_all["spearman_rho"],
        "mantel_p_5000perms": mantel_all["mantel_p_pearson"],
        "mantel_subchecks": {
            "all_shared": mantel_all,
            "all_shared_no_longbranch_outlier": mantel_no_lb,
            "arch3_only_zero_overlap": mantel_a3_only,
            "arch1_vs_arch3_cross": mantel_a1_vs_a3,
            "arch1_only_shared_dna": mantel_a1_only,
        },
        "fitch_parsimony": {
            "full_tree_changes": fitch_n,
            "full_root_states": sorted(fitch_root),
            "prokaryote_ingroup_changes": fitch_prok_n,
            "prokaryote_root_states": sorted(fitch_prok_root),
        },
        "key_clade_supports": key_supports,
        "sensitivity_checks": sensitivity,
        "patristic_pairs": pairs_data,
    }
    out_path = pdir / "ml_bootstrap_summary.json"
    out_path.write_text(json.dumps(summary, indent=2))
    print(f"\nExported 100% programmatically verified phylogeny summary to {out_path}")


if __name__ == "__main__":
    main()
