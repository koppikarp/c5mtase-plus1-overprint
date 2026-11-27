"""Figure 3, Panel a: Conceptual competing evolutionary hypotheses tested by IQ-TREE 2 ML phylogenies & AU tests."""
import json
from matplotlib.patches import FancyBboxPatch
from style import DATA_DIR, PAL


def _pill(ax, x, y, w, h, fc, txt, tc="white", fs=5.7):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.12",
                                facecolor=fc, edgecolor="none", zorder=3))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, fontweight="bold", color=tc, zorder=4)


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "phylogeny/ml_bootstrap_summary.json").read_text())
    ks = d["key_clade_supports"]
    b_a2 = int(round(ks["Arch2_monophyly_tree1"]["ufboot"]))
    b_a3 = int(round(ks["Arch3_monophyly_tree1"]["ufboot"]))
    n_t1 = d["n_taxa_tree1"]
    au = {h["id"]: h["p_au"] for h in d["au_topology_tests"]["hypotheses"]}

    ax.set_xlim(0, 10.0)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    ax.text(5.0, 5.92, "Phylogenetic & AU test of evolutionary direction", ha="center", va="center",
            fontsize=6.5, fontweight="bold", color=PAL["ink"])

    ax.add_patch(FancyBboxPatch((0.20, 0.35), 4.45, 5.10, boxstyle="round,pad=0,rounding_size=0.22",
                                facecolor="white", ec=PAL["ink"], lw=0.85, zorder=1))
    ax.text(2.42, 5.08, "Alternative models", ha="center", va="center", fontsize=6.3, fontweight="bold", color=PAL["muted"])
    ax.text(2.42, 4.62, "Split (Arch 2) or fusion\n(Arch 3) is basal", ha="center", va="center", fontsize=5.8, color=PAL["muted"])

    ax.plot([0.85, 1.55], [2.85, 2.85], lw=0.9, color=PAL["muted"])
    ax.plot([1.55, 1.55], [2.10, 3.60], lw=0.9, color=PAL["muted"])
    ax.plot([1.55, 2.35], [3.60, 3.60], lw=0.9, color=PAL["arch2"])
    ax.plot([1.55, 2.35], [2.10, 2.10], lw=0.9, color=PAL["plus1"])
    _pill(ax, 2.42, 3.30, 1.95, 0.60, PAL["arch2"], "Arch 2 / 3", fs=5.6)
    _pill(ax, 2.42, 1.80, 1.95, 0.60, PAL["plus1"], "Arch 1 (+1)", fs=5.6)
    ax.text(2.42, 0.90, f"Rejected (AU test):\nArch 3 P = {au['H1_arch3_basal']:.3f}; Arch 2 P = {au['H2_arch2_basal']:.3f}",
            ha="center", va="center", fontsize=5.4, style="italic", color=PAL["muted"])

    ax.add_patch(FancyBboxPatch((5.15, 0.35), 4.65, 5.10, boxstyle="round,pad=0,rounding_size=0.22",
                                facecolor="white", ec=PAL["plus1"], lw=0.95, zorder=1))
    ax.text(7.48, 5.08, "Supported model", ha="center", va="center", fontsize=6.3, fontweight="bold", color=PAL["plus1"])
    ax.text(7.48, 4.62, "Arch 1 (+1) arose first;\nArch 2 & 3 nested within", ha="center", va="center", fontsize=5.8, color=PAL["ink"])

    ax.plot([5.55, 6.10], [2.85, 2.85], lw=1.0, color=PAL["plus1"])
    ax.plot([6.10, 6.10], [1.95, 3.75], lw=1.0, color=PAL["plus1"])
    ax.plot([6.10, 7.20], [3.75, 3.75], lw=1.0, color=PAL["plus1"])
    ax.plot([6.65, 6.65], [3.15, 3.75], lw=1.0, color=PAL["arch2"])
    ax.plot([6.65, 7.20], [3.15, 3.15], lw=1.0, color=PAL["arch2"])
    ax.plot([6.10, 7.20], [1.95, 1.95], lw=1.0, color=PAL["plus1"])
    ax.plot([6.65, 6.65], [1.95, 2.52], lw=1.0, color=PAL["arch3"])
    ax.plot([6.65, 7.20], [2.52, 2.52], lw=1.0, color=PAL["arch3"])

    _pill(ax, 7.25, 3.52, 2.30, 0.46, PAL["plus1"], "Arch 1 FPCQ", fs=5.4)
    _pill(ax, 7.25, 2.92, 2.30, 0.46, PAL["arch2"], f"Arch 2 (UF {b_a2}%)", fs=5.3)
    _pill(ax, 7.25, 2.29, 2.30, 0.46, PAL["arch3"], f"Arch 3 (UF {b_a3}%)", fs=5.3)
    _pill(ax, 7.25, 1.72, 2.30, 0.46, PAL["plus1"], "Arch 1 TPCQ", fs=5.4)

    ax.text(7.48, 0.90, f"Rooted on Human DNMT1 +\nTRDMT1 ({n_t1} taxa, Q.pfam+R6)", ha="center", va="center",
            fontsize=5.5, color=PAL["ink"])
