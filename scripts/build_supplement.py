#!/usr/bin/env python3
"""Deterministically build Supplementary Figure S1 (full uncollapsed ML trees) and
Supplementary Information LaTeX/PDF for the C5-MTase +1 overprint preprint directly
from verified repository datasets.
"""
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PREPRINT = ROOT / "preprint"
RENDERED = PREPRINT / "figures/rendered"

PAL = {
    "ARCH1_FPCQ": "#C85A3F",
    "ARCH1_TPCQ": "#B03A2E",
    "ARCH1_SPCQ": "#D97B48",
    "ARCH2": "#8E4B78",
    "ARCH2_SPLIT_TRD": "#A36490",
    "ARCH3": "#22859A",
    "NOp1": "#7E848C",
    "OUT": "#4A6FA5",
    "ink": "#111111",
    "muted": "#5C6168",
}


# ---------------------------------------------------------------------------
# Minimal Newick parser & rerooting
# ---------------------------------------------------------------------------
class Node:
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None

    def is_leaf(self):
        return len(self.children) == 0

    def leaves(self):
        if self.is_leaf():
            return [self]
        out = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def parse_newick(text):
    s = text.strip().rstrip(";")
    idx = [0]

    def parse_node():
        node = Node()
        if idx[0] < len(s) and s[idx[0]] == "(":
            idx[0] += 1
            while True:
                ch = parse_node()
                ch.parent = node
                node.children.append(ch)
                if idx[0] < len(s) and s[idx[0]] == ",":
                    idx[0] += 1
                elif idx[0] < len(s) and s[idx[0]] == ")":
                    idx[0] += 1
                    break
                else:
                    break
        # Read label and branch length
        tok = []
        while idx[0] < len(s) and s[idx[0]] not in ",)":
            tok.append(s[idx[0]])
            idx[0] += 1
        label_str = "".join(tok).strip()
        if ":" in label_str:
            nm, bl = label_str.rsplit(":", 1)
            node.name = nm.strip()
            try:
                node.length = max(0.0, float(bl.strip()))
            except ValueError:
                node.length = 0.0
        else:
            node.name = label_str
        return node

    return parse_node()


def reroot_on_outgroup(root, outgroup_pred):
    """Reroot an unrooted tree at the midpoint of the branch leading to the outgroup clade."""
    all_nodes = []

    def collect(n):
        all_nodes.append(n)
        for c in n.children:
            collect(c)

    collect(root)
    all_leaves = set(root.leaves())
    og_leaves = {l for l in all_leaves if outgroup_pred(l.name)}
    if not og_leaves or og_leaves == all_leaves:
        return root

    # Find node whose descendant leaf set equals og_leaves (or its complement)
    target = None
    for n in all_nodes:
        if n is root:
            continue
        ls = set(n.leaves())
        if ls == og_leaves or (all_leaves - ls) == og_leaves:
            target = n
            break
    if target is None:
        # Fallback: first outgroup leaf
        target = next(iter(og_leaves))

    # Build adjacency graph preserving internal edge support labels
    adj = {}

    def add_edge(u, v, w, lbl):
        adj.setdefault(u, []).append((v, w, lbl))
        adj.setdefault(v, []).append((u, w, lbl))

    for n in all_nodes:
        for c in n.children:
            add_edge(n, c, c.length, c.name if not c.is_leaf() else "")

    p = target.parent
    w = target.length
    el = target.name if not target.is_leaf() else ""
    adj[p] = [x for x in adj[p] if x[0] is not target]
    adj[target] = [x for x in adj[target] if x[0] is not p]

    new_root = Node("ROOT", 0.0)
    adj[new_root] = [(target, w / 2.0, el), (p, w / 2.0, el)]
    adj[target].append((new_root, w / 2.0, el))
    adj[p].append((new_root, w / 2.0, el))

    def rebuild(curr, prev, bl, elbl):
        nn = Node(curr.name if curr.is_leaf() else elbl, bl)
        for nxt, nw, nlbl in adj.get(curr, []):
            if nxt is prev:
                continue
            ch = rebuild(nxt, curr, nw, nlbl)
            ch.parent = nn
            nn.children.append(ch)
        return nn

    r = rebuild(new_root, None, 0.0, "")
    # Collapse degree-2 internal nodes (other than root)
    def collapse(n):
        new_ch = []
        for c in n.children:
            collapse(c)
            if len(c.children) == 1:
                gc = c.children[0]
                gc.length += c.length
                if not gc.name and c.name:
                    gc.name = c.name
                gc.parent = n
                new_ch.append(gc)
            else:
                new_ch.append(c)
        n.children = new_ch

    collapse(r)
    return r


def ladderize(node):
    if node.is_leaf():
        return 1
    sizes = [(ladderize(c), c) for c in node.children]
    sizes.sort(key=lambda x: x[0])
    node.children = [c for _, c in sizes]
    return sum(s for s, _ in sizes)


def tip_color(name):
    if name.startswith("OUT_"):
        return PAL["OUT"]
    if name.startswith("A3_") or "Arch3" in name:
        return PAL["ARCH3"]
    if name.startswith("A2_"):
        return PAL["ARCH2"]
    if name.startswith("SplitTRD"):
        return PAL["ARCH2_SPLIT_TRD"]
    if name.startswith("NOp1"):
        return PAL["NOp1"]
    if "TPCQ" in name:
        return PAL["ARCH1_TPCQ"]
    if "SPCQ" in name or "SPCT" in name:
        return PAL["ARCH1_SPCQ"]
    return PAL["ARCH1_FPCQ"]


def plot_tree(ax, root, title, max_x_pad=0.42, font_sz=4.3, cap_bl=None):
    ladderize(root)
    leaves = root.leaves()
    n_leaves = len(leaves)
    y_map = {lf: i for i, lf in enumerate(leaves)}
    x_map = {}

    def calc_pos(n, x):
        bl = min(n.length, cap_bl) if (cap_bl is not None and n.is_leaf()) else n.length
        xn = x + bl
        x_map[n] = xn
        if not n.is_leaf():
            for c in n.children:
                calc_pos(c, xn)
            y_map[n] = 0.5 * (y_map[n.children[0]] + y_map[n.children[-1]])

    calc_pos(root, 0.0)
    xmax = max(x_map.values())

    def draw_node(n):
        xn, yn = x_map[n], y_map[n]
        if n.is_leaf():
            col = tip_color(n.name)
            ax.scatter([xn], [yn], s=5.5, color=col, zorder=4, linewidths=0)
            lbl = n.name
            if cap_bl is not None and n.length > cap_bl:
                lbl += f" [bl={n.length:.2f}]"
            ax.text(xn + xmax * 0.015, yn, lbl, fontsize=font_sz, va="center", ha="left",
                    color=col if col != PAL["NOp1"] else "#555B62", family="Arial")
            return
        ys = [y_map[c] for c in n.children]
        ax.plot([xn, xn], [min(ys), max(ys)], color="#222222", lw=0.45, zorder=2)
        for c in n.children:
            xc, yc = x_map[c], y_map[c]
            col = tip_color(c.name) if c.is_leaf() else "#222222"
            ls = "--" if (cap_bl is not None and c.is_leaf() and c.length > cap_bl) else "-"
            ax.plot([xn, xc], [yc, yc], color=col if c.is_leaf() else "#222222", lw=0.48, ls=ls, zorder=2)
            draw_node(c)
        # IQ-TREE 2 SH-aLRT / UFBoot2 support label for well-supported internal nodes
        if n is not root and n.name:
            if "/" in n.name:
                try:
                    sh_v, uf_v = [float(x) for x in n.name.split("/", 1)]
                    if uf_v >= 70 or sh_v >= 80:
                        ax.text(xn - xmax * 0.008, yn + 0.16, f"{int(round(sh_v))}/{int(round(uf_v))}",
                                fontsize=3.25, ha="right", va="bottom", color="#555555", family="Arial")
                except ValueError:
                    pass
            else:
                try:
                    sh = float(n.name)
                    if sh >= 0.70:
                        ax.text(xn - xmax * 0.008, yn + 0.18, f"{int(round(sh * 100))}",
                                fontsize=3.5, ha="right", va="bottom", color="#555555", family="Arial")
                except ValueError:
                    pass

    draw_node(root)
    # Scale bar
    sb = 0.5
    ax.plot([0.0, sb], [-2.2, -2.2], color="#111111", lw=0.9)
    ax.text(sb / 2.0, -3.2, f"{sb:.1f} sub/site", fontsize=4.8, ha="center", va="top", color="#111111", family="Arial")
    ax.set_xlim(-xmax * 0.03, xmax * (1.0 + max_x_pad))
    ax.set_ylim(-4.5, n_leaves + 0.8)
    ax.axis("off")
    ax.set_title(title, fontsize=7.5, fontweight="bold", loc="left", pad=4, family="Arial")


def build_supp_fig1():
    t1 = parse_newick((DATA / "phylogeny/tree1_frame0_cat_iqtree.treefile").read_text())
    t1 = reroot_on_outgroup(t1, lambda nm: nm.startswith("OUT_Human_"))
    t2 = parse_newick((DATA / "phylogeny/tree2_cslsdile_domain_iqtree.treefile").read_text())
    t2 = reroot_on_outgroup(t2, lambda nm: nm.startswith("A3_") or nm == "A2_UP_A0A1I2_Paeni")

    plt.rcParams["font.family"] = "Arial"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 9.8), gridspec_kw={"width_ratios": [1.14, 1.0]})
    fig.subplots_adjust(left=0.02, right=0.99, top=0.94, bottom=0.05, wspace=0.06)

    plot_tree(ax1, t1, "a   Frame-0 C5-MTase IQ-TREE 2 (Tree 1, n = 97, Q.pfam+R6)",
              max_x_pad=0.58, font_sz=4.15, cap_bl=2.2)
    plot_tree(ax2, t2, "b   CSLSDILE/KYFLS IQ-TREE 2 (Tree 2, n = 67, VT+G4)",
              max_x_pad=0.62, font_sz=4.35, cap_bl=1.8)

    # Shared top legend
    handles = [
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["ARCH1_FPCQ"], label="Arch 1 FPCQ"),
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["ARCH1_TPCQ"], label="Arch 1 TPCQ"),
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["ARCH1_SPCQ"], label="Arch 1 SPCQ/SPCT"),
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["ARCH2"], label="Arch 2 Split-IX"),
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["ARCH2_SPLIT_TRD"], label="Split-TRD/IX/X"),
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["ARCH3"], label="Arch 3 C-term fusion"),
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["NOp1"], label="+1-less control"),
        plt.Line2D([], [], marker="o", ls="none", ms=4.5, color=PAL["OUT"], label="Euk. outgroup"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.992), ncol=4,
               frameon=False, fontsize=5.8, handletextpad=0.25, columnspacing=1.0)

    RENDERED.mkdir(parents=True, exist_ok=True)
    fig.savefig(RENDERED / "supp_fig1_full_trees.pdf")
    fig.savefig(RENDERED / "supp_fig1_full_trees.png", dpi=220)
    plt.close(fig)
    print("Wrote supp_fig1_full_trees.pdf/.png")


# ---------------------------------------------------------------------------
# LaTeX helpers
# ---------------------------------------------------------------------------
def esc(s):
    if s is None:
        return "--"
    s = str(s)
    for k, v in [("\\", r"\textbackslash{}"), ("_", r"\_"), ("%", r"\%"),
                 ("&", r"\&"), ("#", r"\#"), ("$", r"\$"), (">", r"$>$"), ("<", r"$<$")]:
        s = s.replace(k, v)
    return s


def fmt_p(p):
    if p is None:
        return "--"
    if p < 1e-3:
        m, e = f"{p:.2e}".split("e")
        return f"${m} \\times 10^{{{int(e)}}}$"
    return f"{p:.4f}"


def clean_tax(s):
    if not s:
        return "--"
    return re.sub(r"^[dpocfgs]__", "", str(s))


def arch_label(a):
    return {
        "Arch1_intact": "Arch 1 (intact)",
        "Arch1_interrupted": "Arch 1 (interrupted)",
        "Arch2_split_IX": "Arch 2 (split IX)",
        "Arch2_split_TRD_IX_X": "Split TRD/IX/X",
        "Arch3_inframe_cterm": "Arch 3 (C-term fusion)",
    }.get(a, a)


def build_supplement_tex():
    gtdb = json.loads((DATA / "gtdb2339_c5mtase_plus1_census.json").read_text())
    phage = json.loads((DATA / "inphared_phage_c5mtase_census.json").read_text())
    uniref = json.loads((DATA / "uniref100_deletion_audit.json").read_text())
    tpcq = json.loads((DATA / "uniprot_tpcq_dna_audit.json").read_text())
    sel = json.loads((DATA / "selection_4d_and_stop_stats.json").read_text())
    phy = json.loads((DATA / "phylogeny/ml_bootstrap_summary.json").read_text())
    af3 = json.loads((DATA / "structures/trans_vs_cis_3d_audit.json").read_text())

    # Assemble all 45 GTDB loci
    gtdb_rows = []
    for r in gtdb["overprinted"]:
        if r.get("h0_IX"):
            m9_loc = "Frame 0"
        elif r.get("p1_IX") or r["refined_arch"] == "Arch2_split_IX":
            m9_loc = "Frame +1"
        elif r["refined_arch"] == "Arch2_split_TRD_IX_X":
            m9_loc = "Frame +1 (TRD--X)"
        else:
            m9_loc = "Frame 0 (div.)"
        gtdb_rows.append((
            r["acc"], clean_tax(r["domain"]), clean_tax(r["phylum"]), clean_tax(r["order"]),
            r.get("motif4", "--"), arch_label(r["refined_arch"]),
            f"{r['hlen']}", f"{r['p1_len']}", f"{r['gc']*100:.1f}",
            m9_loc,
        ))
    for r in gtdb["other"]:
        gtdb_rows.append((
            r["acc"], clean_tax(r["domain"]), clean_tax(r["phylum"]), clean_tax(r["order"]),
            r.get("motif4", "--"), arch_label(r["refined_arch"]),
            "--" if "interrupted" in r["refined_arch"] else "Full", f"{r['len']}", f"{r['gc']*100:.1f}",
            "Frame +1" if r["refined_arch"] == "Arch2_split_IX" else "Frame 0",
        ))
    for r in gtdb["inframe_cterm"]:
        gtdb_rows.append((
            r["acc"], clean_tax(r["domain"]), clean_tax(r["phylum"]), clean_tax(r["order"]),
            "TPCQ/fusion", arch_label(r["refined_arch"]),
            f"{r['len']}", "In-frame", f"{r['gc']*100:.1f}", "Frame 0",
        ))
    assert len(gtdb_rows) == 45, f"Expected 45 GTDB loci, got {len(gtdb_rows)}"
    gtdb_rows.sort(key=lambda x: (x[5], x[2], x[0]))

    # Assemble all 39 phage loci
    phage_rows = []
    for r in phage["overprinted"]:
        phage_rows.append((
            r["acc"], r.get("motif4", "--"), arch_label(r["refined_arch"]),
            f"{r['hlen']}", f"{r['p1_len']}",
            "Yes" if r.get("h0_IX") else "Div/Split",
            "Yes" if r.get("h0_X") else "No",
            f"{r['gc']*100:.1f}",
            f"{r.get('pid', 0):.1f}",
        ))
    assert len(phage_rows) == 39, f"Expected 39 phage loci, got {len(phage_rows)}"
    phage_rows.sort(key=lambda x: (x[1], x[2], x[0]))

    # Table S1 body
    t1_lines = "\n".join(
        " & ".join(esc(v) for v in row) + r" \\"
        for row in gtdb_rows
    )
    # Table S2 body
    t2_lines = "\n".join(
        " & ".join(esc(v) for v in row) + r" \\"
        for row in phage_rows
    )
    # Table S3a body (33 UniProtKB audited)
    t3a_lines = "\n".join(
        f"{e['mmseqs2_rank']} & {esc(e['accession'])} & {esc(e['status'])} & "
        f"{e['identity']*100:.1f}\\% & {e['q_cov']:.1f}\\% & {e['e_value']:.1e} & "
        f"{esc(e['organism']) if e['organism'] else '--'} & 0 / 0 \\\\"
        for e in sorted(uniref["entries"], key=lambda x: x["mmseqs2_rank"])
    )
    # Table S3b body (6 DNA-verified UniProt entries)
    t3b_lines = "\n".join(
        f"{esc(x['acc'])} & {esc(x['entry'])} & {esc(x['org'])} & {esc(x['ena_pid'])} & "
        f"{x['f0_len']} & {esc(x['motif4'])} & {'Yes' if x['has_plus1'] else 'No (beta)'} & {x['p1_stops_win']} \\\\"
        for x in tpcq
    )

    # Table S4: selection + phylogeny audit
    pv = sel["p_values"]
    dj = sel["disjoint_4d_sensitivity"]
    aa = sel["plus1_aa_preservation"]
    st = sel["stop_purge_5frames"]
    st_low = sel["stop_purge_5frames_low_allgc3"]
    ks = phy["key_clade_supports"]
    ms = phy["mantel_subchecks"]
    sc = phy["sensitivity_checks"]

    t4_rows = [
        ("4D non-stop-risk excess (Arch 1 FPCQ, $n=13$ loci / 51 pairs)", "Inside vs.\\ outside (paired 1-sided Wilcoxon)", "Med. $+0.258$ vs.\\ $+0.015$ (13/13)", fmt_p(pv["fpcq_gene_wilcoxon_in_vs_out"])),
        ("4D non-stop-risk excess (Arch 1 TPCQ/SPCQ, $n=7$ loci / 12 pairs)", "Inside vs.\\ outside (paired 1-sided Wilcoxon)", "Med. $+0.225$ vs.\\ $+0.035$ (6/7)", fmt_p(pv["tspcq_gene_wilcoxon_in_vs_out"])),
        ("4D inside excess: FPCQ vs.\\ +1-less controls ($n=13$ vs.\\ 26)", "Between-group 1-sided Mann--Whitney $U$", "Med. $+0.258$ vs.\\ $+0.030$", fmt_p(pv["fpcq_vs_ctrl_in_excess_mwu_p"])),
        ("4D inside$-$outside $\\Delta$: FPCQ vs.\\ controls ($n=13$ vs.\\ 26)", "Between-group 1-sided Mann--Whitney $U$", "Delta enrichment", fmt_p(pv["fpcq_vs_ctrl_delta_mwu_p"])),
        ("4D inside excess: TPCQ/SPCQ vs.\\ controls ($n=7$ vs.\\ 26)", "Between-group 1-sided Mann--Whitney $U$", "Med. $+0.225$ vs.\\ $+0.030$", fmt_p(pv["tspcq_vs_ctrl_in_excess_mwu_p"])),
        ("4D disjoint non-overlapping pairs (200 random matchings)", "Zero genome reuse (paired 1-sided Wilcoxon)", f"Med. frac in$>$out {dj['random_matchings_median_frac_in_gt_out']*100:.1f}\\%", fmt_p(dj["random_matchings_median_p_greater"])),
        ("4D greedy disjoint FPCQ vs.\\ disjoint controls (6 vs.\\ 10 pairs)", "Zero genome reuse (1-sided MWU)", f"Med. $+{dj['greedy_med_in_excess']:.3f}$ inside", fmt_p(dj["greedy_vs_ctrl_in_mwu_p_greater"])),
        ("5-frame stop purge: Frame +1 vs.\\ flanking null ($n=30$ genomes)", "Exact 1-sided paired Wilcoxon (1,000 null reps)", "Obs. med. 0.0 vs.\\ null 4.7 (30/30)", fmt_p(pv["stop_purge_plus1_vs_flanking_null_p_1s"])),
        ("5-frame stop purge: low all-codon $\\text{GC3} \\le 0.45$ ($n=14$)", "Exact 1-sided paired Wilcoxon vs.\\ flanking null", "Obs. med. 0.0 vs.\\ null 6.8 (14/14)", fmt_p(st_low["obs_vs_flanking_null_wilcoxon"]["+1"]["wilcoxon_p_less"])),
        ("5-frame stop purge: Frame +1 vs.\\ all 4 other frames ($n=30$)", "Within-gene max 1-sided paired Wilcoxon", "Max $P$ across $+2, -0, -1, -2$", fmt_p(pv["stop_purge_within_gene_max_p_1s"])),
        ("Frame $-0$ shadowing: \\texttt{TTA}/\\texttt{CTA} vs.\\ \\texttt{TCA} ($n=30$)", "Exact paired Wilcoxon vs.\\ flanking null", "TTA/CTA 1.20 vs.\\ 3.24; TCA 1.20 vs.\\ 1.19", f"{fmt_p(pv['stop_purge_minus0_tta_cta_p_1s'])} / {fmt_p(pv['stop_purge_minus0_tca_p_2s'])}"),
        ("+1 AA preservation at Frame-0-divergent sites ($n=13$ vs.\\ 26)", "Gene-level 1-sided MWU (disjoint 6 vs.\\ 10)", f"{aa['gene_level']['case_median']*100:.1f}\\% vs.\\ {aa['gene_level']['ctrl_median']*100:.1f}\\% ({aa['disjoint_pairs']['case_median']*100:.1f}\\% vs.\\ {aa['disjoint_pairs']['ctrl_median']*100:.1f}\\%)", f"{fmt_p(pv['plus1_aa_preservation_gene_mwu_p'])} ({fmt_p(pv['plus1_aa_preservation_disjoint_mwu_p'])})"),
        ("Phylogeny Tree 1: Arch 3 monophyly ($n=7$ taxa, 4 phyla)", "IQ-TREE 2 Q.pfam+R6 (SH-aLRT / UFBoot2)", f"SH {ks['Arch3_monophyly_tree1']['sh_alrt']:.1f}\\%, UFBoot {ks['Arch3_monophyly_tree1']['ufboot']:.0f}\\%", "Drop-LB: 98.7\\% / 100\\%"),
        ("Phylogeny Tree 1: TPCQ + Arch 3 superclade ($n=29$ taxa)", "IQ-TREE 2 Q.pfam+R6 (SH-aLRT / UFBoot2)", f"SH {ks['TPCQ_Arch3_superclade_tree1']['sh_alrt']:.1f}\\%, UFBoot {ks['TPCQ_Arch3_superclade_tree1']['ufboot']:.0f}\\%", "Drop-LB: 99.1\\% / 100\\%"),
        ("Phylogeny Tree 1: Arch 2 monophyly / parent clade ($n=3$ / 5)", "IQ-TREE 2 Q.pfam+R6 (SH-aLRT / UFBoot2)", f"SH {ks['Arch2_monophyly_tree1']['sh_alrt']:.1f}\\% / {ks['Arch2_parent_clade_tree1']['sh_alrt']:.1f}\\%, UF {ks['Arch2_monophyly_tree1']['ufboot']:.0f}\\% / {ks['Arch2_parent_clade_tree1']['ufboot']:.0f}\\%", "Drop-LB: 63.8\\% / 67\\%"),
        ("Phylogeny Tree 1: Bacillota Split-TRD/IX/X in SPCQ ($n=6$)", "IQ-TREE 2 Q.pfam+R6 (SH-aLRT / UFBoot2)", f"SH {ks['Bacillota_SPCQ_SplitTRD_clade_tree1']['sh_alrt']:.1f}\\%, UFBoot {ks['Bacillota_SPCQ_SplitTRD_clade_tree1']['ufboot']:.0f}\\%", "Drop-LB: 77.2\\% / 76\\%"),
        ("Patristic congruence Tree 1 vs.\\ Tree 2 (all 66 shared taxa)", "Pearson $r$ / Spearman $\\rho$ (5,000 Mantel perms)", f"$r = {ms['all_shared']['pearson_r']:.3f}$, $\\rho = {ms['all_shared']['spearman_rho']:.3f}$", fmt_p(ms["all_shared"]["mantel_p_pearson"])),
        ("Patristic congruence excl.\\ long-branch outlier (65 taxa)", "Pearson $r$ / Spearman $\\rho$ (5,000 Mantel perms)", f"$r = {ms['all_shared_no_longbranch_outlier']['pearson_r']:.3f}$, $\\rho = {ms['all_shared_no_longbranch_outlier']['spearman_rho']:.3f}$", fmt_p(ms["all_shared_no_longbranch_outlier"]["mantel_p_pearson"])),
        ("Patristic congruence Arch 3 only (zero DNA overlap, $n=7$)", "Pearson $r$ / Spearman $\\rho$ (5,000 Mantel perms)", f"$r = {ms['arch3_only_zero_overlap']['pearson_r']:.3f}$, $\\rho = {ms['arch3_only_zero_overlap']['spearman_rho']:.3f}$", fmt_p(ms["arch3_only_zero_overlap"]["mantel_p_pearson"])),
        ("Patristic congruence Arch 1 vs.\\ Arch 3 cross-pairs ($n=385$)", "Pearson $r$ / Spearman $\\rho$ (5,000 Mantel perms)", f"$r = {ms['arch1_vs_arch3_cross']['pearson_r']:.3f}$, $\\rho = {ms['arch1_vs_arch3_cross']['spearman_rho']:.3f}$", fmt_p(ms["arch1_vs_arch3_cross"]["mantel_p_pearson"])),
    ]
    t4_lines = "\n".join(f"{a} & {b} & {c} & {d} \\\\" for a, b, c, d in t4_rows)
    au_rows = []
    for h in phy["au_topology_tests"]["hypotheses"]:
        verdict = "\\textbf{Rejected}" if h["rejected"] else "\\textbf{Supported}"
        au_rows.append(
            f"{esc(h['label'])} & {h['logl']:.2f} & {h['delta_logl']:.2f} & {h['bp_rell']:.4f} & "
            f"{h['p_kh']:.4f} & {h['p_sh']:.4f} & {fmt_p(h['c_elw'])} & {fmt_p(h['p_au'])} & {verdict} \\\\"
        )
    au_lines = "\n".join(au_rows)

    # Table S5: Per-model AF3 benchmark (all 30 models)
    nice_af3 = {
        "posctrl_arch2_eaei_trans": "M.EaeI (Arch 2, pos. ctrl)",
        "phagep1_arch1_dmt_trans_64aa_core": "Phage P1 dmt + 64aa (A1)",
        "phagep1_arch1_dmt_trans_157aa": "Phage P1 dmt + 157aa (A1)",
        "bact_arch1_asoli_trans": "A. soli (Arch 1, Bact.)",
        "arch_arch1_nanob_trans": "Nanobdellota (Arch 1, Arch.)",
        "negctrl_hhai_asoli_trans": "M.HhaI + A. soli +1 (neg. ctrl)",
    }
    t5_rows = []
    for rec in af3["trans_complexes"]:
        nm = nice_af3.get(rec["name"], rec["name"])
        for m in range(5):
            t5_rows.append(
                f"{esc(nm)} & {m} & {rec['iptm_5'][m]:.2f} & {rec['ptm_5'][m]:.2f} & "
                f"{rec['ranking_score_5'][m]:.2f} & {rec['pae_min_5'][m]:.2f} & "
                f"{rec['r1_kyf_dists_5'][m]:.2f} & {rec['com_dist_vs_cis_5'][m]:.2f} & "
                f"{rec['docked_rmsd_vs_cis_5'][m]:.2f} & {rec['heavy_contacts_4_5A_5'][m]} & {(rec['clashes'][m] if isinstance(rec['clashes'], list) else rec['clashes'])} \\\\"
            )
    t5_lines = "\n".join(t5_rows)

    derep = phage["dereplicated_prevalence"]
    derep_note = (
        f"Dereplicated phage prevalence across all 230 C5-MTases: "
        f"$<95\\%$ identity: {derep['lt_95pct']['n_overprinted_clusters']}/{derep['lt_95pct']['n_clusters']} clusters "
        f"({derep['lt_95pct']['pct_overprinted']:.1f}\\%); "
        f"$<90\\%$ identity: {derep['lt_90pct']['n_overprinted_clusters']}/{derep['lt_90pct']['n_clusters']} clusters "
        f"({derep['lt_90pct']['pct_overprinted']:.1f}\\%, TPCQ 5/6 [83.3\\%], FPCQ 10/47 [21.3\\%]); "
        f"$<80\\%$ identity: {derep['lt_80pct']['n_overprinted_clusters']}/{derep['lt_80pct']['n_clusters']} clusters "
        f"({derep['lt_80pct']['pct_overprinted']:.1f}\\%, TPCQ 5/6 [83.3\\%])."
    )

    tex = rf"""\documentclass[10pt,letterpaper]{{article}}

\usepackage[letterpaper,left=0.65in,right=0.65in,top=0.72in,bottom=0.72in,headheight=19pt]{{geometry}}
\usepackage{{fontspec}}
\usepackage{{amsmath,amssymb}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{longtable}}
\usepackage{{array}}
\usepackage{{xcolor}}
\usepackage{{titlesec}}
\usepackage{{fancyhdr}}
\usepackage{{lastpage}}
\usepackage[hidelinks]{{hyperref}}
\urlstyle{{same}}

\setmainfont{{Arial}}[
  Scale=0.92,
  BoldFont={{Arial Bold}},
  ItalicFont={{Arial Italic}},
  BoldItalicFont={{Arial Bold Italic}}
]
\setsansfont{{Arial}}[Scale=0.92]
\setmonofont{{Arial}}[Scale=0.92]

\linespread{{1.08}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{4pt}}
\definecolor{{RuleGray}}{{HTML}}{{BDBDBD}}
\color{{black}}

\titleformat{{\section}}
  {{\fontsize{{11.5}}{{14}}\selectfont\bfseries\color{{black}}}}
  {{}}{{0pt}}{{}}
  [\vspace{{-3pt}}\textcolor{{RuleGray}}{{\rule{{\textwidth}}{{0.5pt}}}}]
\titlespacing*{{\section}}{{0pt}}{{12pt}}{{6pt}}

\pagestyle{{fancy}}
\fancyhf{{}}
\renewcommand{{\headrulewidth}}{{0.4pt}}
\renewcommand{{\headrule}}{{\hbox to\headwidth{{\color{{RuleGray}}\leaders\hrule height \headrulewidth\hfill}}}}
\fancyhead[L]{{\footnotesize\color{{black}}Koppikar (2026) \textbar\ Supplementary Information: Pan-prokaryotic +1 overprinted C5-MTase family}}
\fancyhead[R]{{\footnotesize\color{{black}}Supplementary Information}}
\fancyfoot[C]{{\footnotesize\color{{black}}S\thepage\ / S\pageref{{LastPage}}}}

\newcommand{{\gene}}[1]{{\textit{{#1}}}}
\newcommand{{\spname}}[1]{{\textit{{#1}}}}

\begin{{document}}

\noindent{{\fontsize{{13.5}}{{17}}\selectfont\bfseries\color{{black}}\raggedright
Supplementary Information\\[3pt]
A Pan-Prokaryotic and Viral +1 Overprinted Gene Family Inside Full-Length C5-Cytosine DNA Methyltransferases Explains the Origin of Split-Subunit Methyltransferases\par}}
\vspace{{4pt}}
\noindent{{\small Pratik Koppikar\textsuperscript{{1,*}} \;\textbar\; \textsuperscript{{1}}University of Cambridge, Cambridge, United Kingdom \;\textbar\; \href{{mailto:psk34@cam.ac.uk}}{{psk34@cam.ac.uk}}\\
Repository and raw datasets: \url{{https://github.com/koppikarp/c5mtase-plus1-overprint}}\par}}
\vspace{{3pt}}
\noindent{{\color{{RuleGray}}\rule{{\textwidth}}{{0.7pt}}}}

\section*{{Contents}}
\begin{{itemize}}
  \item \textbf{{Supplementary Figure S1:}} Full uncollapsed IQ-TREE 2 maximum-likelihood phylogenies of Frame-0 C5-MTases (Tree 1, $n=97$ taxa, \texttt{{Q.pfam+R6}}) and the CSLSDILE/KYFLS domain (Tree 2, $n=67$ taxa, \texttt{{VT+G4}}).
  \item \textbf{{Supplementary Table S1:}} Complete census of all 45 verified GTDB R232 loci across 13 prokaryotic phyla.
  \item \textbf{{Supplementary Table S2:}} Complete census of all 39 INPHARED bacteriophage C5-MTase +1 loci and dereplicated prevalence.
  \item \textbf{{Supplementary Table S3:}} UniRef100/UniProtKB database erasure audit ($n=33$) and DNA-level reading-frame audit of UniProt TPCQ entries ($n=6$).
  \item \textbf{{Supplementary Table S4:}} Complete summary of evolutionary selection, stop-purging, IQ-TREE 2 phylogenetic support, AU topology tests, Mantel congruence, and sensitivity checks.
  \item \textbf{{Supplementary Table S5:}} Per-model AlphaFold 3 benchmark across all 30 diffusion models (6 complexes $\times$ 5 models).
\end{{itemize}}

\newpage
\section*{{Supplementary Figure S1}}
\begin{{center}}
  \includegraphics[width=\textwidth,height=0.75\textheight,keepaspectratio]{{figures/rendered/supp_fig1_full_trees.pdf}}
\end{{center}}
\vspace{{2pt}}
\noindent\begin{{minipage}}{{\textwidth}}\footnotesize\raggedright
\textbf{{Supplementary Figure S1. Full uncollapsed IQ-TREE 2 maximum-likelihood phylogenies supporting Figure 3.}} \textbf{{a}}, Complete DNMT1/TRDMT1-rooted \texttt{{Q.pfam+R6}} ML phylogeny of the Frame-0 C5-MTase catalytic core (Tree 1; 97 unique non-redundant taxa, 185 MAFFT L-INS-i columns). The extreme long-branch taxon \texttt{{A1\_SPCT\_Baci\_10645}} (branch length 3.46) is truncated visually with a dashed branch and annotated with its true length; removing this taxon in sensitivity re-inference (96 taxa, \texttt{{Q.pfam+R6}}) yields Arch 3 monophyly SH-aLRT 98.7\% / UFBoot2 100\% and TPCQ+Arch 3 superclade SH-aLRT 99.1\% / UFBoot2 100\%. \textbf{{b}}, Complete harmonized \texttt{{VT+G4}} ML phylogeny of the CSLSDILE/KYFLS domain (Tree 2; 67 taxa, 77 trimmed columns), rooted on the Arch 3 clade for display (note that the single non-shared long-branch taxon \texttt{{A2\_UP\_A0A1I2\_Paeni}}, bl 2.16, clusters with Arch 3 via long-branch attraction in the 77-aa alignment; 8-taxon clade SH-aLRT 96.5\% / UFBoot2 90\%). Internal node labels report SH-aLRT / UFBoot2 (\texttt{{-bnni}}) support values where SH-aLRT $\ge 80\%$ or UFBoot2 $\ge 70\%$.
\end{{minipage}}

\newpage
\section*{{Supplementary Table S1: Complete Census of All 45 Verified GTDB R232 Loci}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.2pt}}\renewcommand{{\arraystretch}}{{1.10}}
\begin{{longtable}}{{@{{}}llllllllll@{{}}}}
\toprule
\textbf{{Genome}} & \textbf{{Domain}} & \textbf{{Phylum}} & \textbf{{Order}} & \textbf{{Motif IV}} & \textbf{{Architecture}} & \textbf{{F0 aa}} & \textbf{{+1 aa}} & \textbf{{GC\%}} & \textbf{{Motif IX}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Genome}} & \textbf{{Domain}} & \textbf{{Phylum}} & \textbf{{Order}} & \textbf{{Motif IV}} & \textbf{{Architecture}} & \textbf{{F0 aa}} & \textbf{{+1 aa}} & \textbf{{GC\%}} & \textbf{{Motif IX}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t1_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S2: Complete Census of All 39 INPHARED Bacteriophage +1 Loci}}
\noindent{{\footnotesize {derep_note}\par}}
\vspace{{4pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{4.2pt}}\renewcommand{{\arraystretch}}{{1.10}}
\begin{{longtable}}{{@{{}}lllllllll@{{}}}}
\toprule
\textbf{{GenBank Acc.}} & \textbf{{Motif IV}} & \textbf{{Architecture}} & \textbf{{Frame-0 aa}} & \textbf{{+1 ORF aa}} & \textbf{{Frame-0 IX}} & \textbf{{Frame-0 X}} & \textbf{{GC\%}} & \textbf{{Seed ID\%}} \\
\midrule
\endfirsthead
\toprule
\textbf{{GenBank Acc.}} & \textbf{{Motif IV}} & \textbf{{Architecture}} & \textbf{{Frame-0 aa}} & \textbf{{+1 ORF aa}} & \textbf{{Frame-0 IX}} & \textbf{{Frame-0 X}} & \textbf{{GC\%}} & \textbf{{Seed ID\%}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t2_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S3: Database Erasure Audit and DNA-Level Verification}}
\noindent{{\footnotesize\textbf{{Table S3a.}} All 33 audited non-UniParc UniProtKB matches among top MMseqs2 UniRef100 hits (query: \spname{{Aneurinibacillales}} +1 polypeptide, GCA\_002355375.1): 29/33 (87.9\%) marked deleted/obsolete; 4/33 active, all with 0 Pfam and 0 InterPro hits.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.5pt}}\renewcommand{{\arraystretch}}{{1.08}}
\begin{{longtable}}{{@{{}}rlllrrll@{{}}}}
\toprule
\textbf{{Rank}} & \textbf{{Accession}} & \textbf{{UniProtKB Status}} & \textbf{{ID\%}} & \textbf{{Q.Cov.}} & \textbf{{E-value}} & \textbf{{Organism (if active)}} & \textbf{{Pfam / InterPro}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Rank}} & \textbf{{Accession}} & \textbf{{UniProtKB Status}} & \textbf{{ID\%}} & \textbf{{Q.Cov.}} & \textbf{{E-value}} & \textbf{{Organism (if active)}} & \textbf{{Pfam / InterPro}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t3a_lines}
\end{{longtable}}
}}

\vspace{{8pt}}
\noindent{{\footnotesize\textbf{{Table S3b.}} DNA-level reading-frame audit of five UniProt TPCQ entries initially appearing +1-less in protein records (5/5 carry intact unannotated +1 CSLSDILE/KYFLS ORFs with 0 internal stops) plus the Arch 2 $\beta$ reference \texttt{{A0A1I2CIH5}}.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.2pt}}\renewcommand{{\arraystretch}}{{1.10}}
\begin{{longtable}}{{@{{}}llllrlll@{{}}}}
\toprule
\textbf{{UniProt}} & \textbf{{Entry}} & \textbf{{Organism}} & \textbf{{ENA/GenBank}} & \textbf{{F0 aa}} & \textbf{{Motif IV}} & \textbf{{Intact +1?}} & \textbf{{+1 Stops}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t3b_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S4: Evolutionary Selection, Phylogenetic Support, and AU Topology Tests}}
\noindent{{\footnotesize\textbf{{Table S4a.}} Summary of 4D synonymous conservation, 5-frame stop purging, IQ-TREE 2 clade support (1,000 SH-aLRT / 1,000 UFBoot2 with \texttt{{-bnni}}), Mantel congruence, and drop-long-branch sensitivity checks.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.5pt}}\renewcommand{{\arraystretch}}{{1.12}}
\begin{{longtable}}{{@{{}}p{{0.31\textwidth}}p{{0.25\textwidth}}p{{0.24\textwidth}}p{{0.15\textwidth}}@{{}}}}
\toprule
\textbf{{Analysis / Cohort}} & \textbf{{Statistical Test / Model}} & \textbf{{Observed Effect}} & \textbf{{$P$-value / Check}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Analysis / Cohort}} & \textbf{{Statistical Test / Model}} & \textbf{{Observed Effect}} & \textbf{{$P$-value / Check}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t4_lines}
\end{{longtable}}
}}

\vspace{{6pt}}
\noindent{{\footnotesize\textbf{{Table S4b.}} Approximately Unbiased (AU) topology tests on Tree 1 (97 taxa, \texttt{{Q.pfam+R6}}, 10,000 RELL resamplings in IQ-TREE 2) comparing the unconstrained ML tree against constrained alternative evolutionary hypotheses.\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{2.8pt}}\renewcommand{{\arraystretch}}{{1.12}}
\begin{{longtable}}{{@{{}}lrrrrrrrr@{{}}}}
\toprule
\textbf{{Topology Hypothesis}} & \textbf{{$\ln L$}} & \textbf{{$\Delta\ln L$}} & \textbf{{bp-RELL}} & \textbf{{$p$-KH}} & \textbf{{$p$-SH}} & \textbf{{c-ELW}} & \textbf{{$p$-AU}} & \textbf{{Result}} \\
\midrule
\endhead
\bottomrule
\endfoot
{au_lines}
\end{{longtable}}
}}

\newpage
\section*{{Supplementary Table S5: Per-Model AlphaFold 3 Benchmark Across All 30 Diffusion Models}}
\noindent{{\footnotesize All six two-chain complexes (five diffusion models each, seed 42; 0/30 steric clashes) evaluated programmatically against \textit{{cis}} \texttt{{AF-A0A6N9NLZ4-F1}} (Motif IX R290--Y411 = 4.02\,\AA) and rigid-core superposed M.HhaI--DNA (PDB 1MHT; 0.84\,\AA\ core RMSD).\par}}
\vspace{{3pt}}
{{\scriptsize\setlength{{\tabcolsep}}{{3.5pt}}\renewcommand{{\arraystretch}}{{1.08}}
\begin{{longtable}}{{@{{}}lrrrrrrrrrr@{{}}}}
\toprule
\textbf{{Complex}} & \textbf{{Mod.}} & \textbf{{ipTM}} & \textbf{{pTM}} & \textbf{{Rank}} & \textbf{{MinPAE}} & \textbf{{R1--F/Y}} & \textbf{{COM}} & \textbf{{DockRMSD}} & \textbf{{Contacts}} & \textbf{{Clash}} \\
\midrule
\endfirsthead
\toprule
\textbf{{Complex}} & \textbf{{Mod.}} & \textbf{{ipTM}} & \textbf{{pTM}} & \textbf{{Rank}} & \textbf{{MinPAE}} & \textbf{{R1--F/Y}} & \textbf{{COM}} & \textbf{{DockRMSD}} & \textbf{{Contacts}} & \textbf{{Clash}} \\
\midrule
\endhead
\bottomrule
\endfoot
{t5_lines}
\end{{longtable}}
}}

\end{{document}}
"""
    out_tex = PREPRINT / "c5mtase_plus1_supplement.tex"
    out_tex.write_text(tex)
    print(f"Wrote {out_tex} ({len(gtdb_rows)} GTDB, {len(phage_rows)} phage, {len(uniref['entries'])} UniRef, 30 AF3 models)")


if __name__ == "__main__":
    build_supp_fig1()
    build_supplement_tex()
