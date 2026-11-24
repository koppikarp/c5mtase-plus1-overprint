"""Figure 3, Panel b: DNMT-rooted LG+Gamma ML Phylogeny (Tree 1, 97 unique taxa collapsed to exact clades)."""
import json
from matplotlib.patches import Polygon
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "phylogeny/ml_bootstrap_summary.json").read_text())
    ks = d["key_clade_supports"]
    n_t1 = d["n_taxa_tree1"]
    a3_b, a3_sh = int(round(ks["Arch3_monophyly_tree1"]["ml_boot"])), ks["Arch3_monophyly_tree1"]["sh"]
    tp_b, tp_sh = int(round(ks["TPCQ_Arch3_superclade_tree1"]["ml_boot"])), ks["TPCQ_Arch3_superclade_tree1"]["sh"]
    a2_b, a2_sh = int(round(ks["Arch2_monophyly_tree1"]["ml_boot"])), ks["Arch2_monophyly_tree1"]["sh"]
    a2p_b, a2p_sh = int(round(ks["Arch2_parent_clade_tree1"]["ml_boot"])), ks["Arch2_parent_clade_tree1"]["sh"]

    ax.set_xlim(-0.3, 14.8)
    ax.set_ylim(-0.6, 12.8)
    ax.axis("off")

    ax.text(6.8, 12.52, f"Frame-0 C5-MTase ML phylogeny ({n_t1} unique taxa, LG+Γ, 100 thorough ML bootstraps)",
            ha="center", va="center", fontsize=6.5, fontweight="bold", color=PAL["ink"])

    # Exact DNMT-rooted partition of all 97 unique taxa (sum: 2 + 11 + 3 + 5 + 8 + 21 + 3 + 3 + 5 + 6 + 1 + 22 + 7 = 97!):
    def branch(x0, y0, x1, y1, col=PAL["ink"], lw=0.9):
        ax.plot([x0, x0, x1], [y0, y1, y1], color=col, lw=lw, solid_capstyle="round", zorder=2)

    def wedge(x, y, dx, dy, fc, label, sub="", tc=PAL["ink"], bold=False):
        tri = Polygon([(x, y), (x + dx, y + dy), (x + dx, y - dy)],
                      closed=True, facecolor=fc, edgecolor=PAL["ink"], linewidth=0.5, zorder=3)
        ax.add_patch(tri)
        ax.text(x + dx + 0.16, y + (0.12 if sub else 0.0), label,
                ha="left", va="center", fontsize=5.8, fontweight="bold" if bold else "normal", color=tc)
        if sub:
            ax.text(x + dx + 0.16, y - 0.23, sub, ha="left", va="center", fontsize=5.15, color=PAL["muted"])

    def badge(x, y, txt, col=PAL["ink"]):
        ax.text(x, y, txt, ha="right", va="bottom", fontsize=5.0, fontweight="bold", color=col,
                bbox=dict(boxstyle="round,pad=0.11", fc="white", ec="none", alpha=0.88), zorder=5)

    # Root at x=0.2, y=3.6
    branch(0.2, 3.6, 4.6, -0.1, PAL["ctrl"], lw=0.95)
    wedge(4.6, -0.1, 0.75, 0.20, PAL["ctrl_lt"], "Outgroup: Human DNMT1, TRDMT1 (n = 2)", tc=PAL["muted"])

    branch(0.2, 3.6, 0.9, 5.1, PAL["ink"], lw=0.95)
    badge(0.88, 5.20, "ML 95%")

    # Basal +1-less phages (n = 11; ML 95% / SH 100%)
    branch(0.9, 5.1, 4.6, 0.75, PAL["ctrl"], lw=0.9)
    badge(3.2, 0.85, "ML 95% / SH 100%", PAL["ctrl"])
    wedge(4.6, 0.75, 0.95, 0.28, PAL["ctrl"], "+1-less phage outgroup (n = 11)", tc=PAL["muted"])

    # Main radiation (n = 84) -> Clade B (n = 47, upper) & Clade A (n = 37, lower)
    branch(0.9, 5.1, 1.5, 6.2, PAL["ink"], lw=0.95)

    # === UPPER CLADE (n = 47, SH 99.1%) ===
    branch(1.5, 6.2, 2.25, 8.75, PAL["ink"], lw=1.0)
    ax.text(2.15, 8.62, "SH 99.1%", ha="right", va="top", fontsize=4.9, fontweight="bold", color=PAL["ink"],
            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.90), zorder=5)

    # Arch 1 FPCQ + Arch 2 (n = 6, SH 90.5%; inner 5-taxon Arch 2 + 2 Arch 1: ML 77% / SH 91.9%)
    branch(2.25, 8.75, 3.45, 6.25, PAL["plus1"], lw=0.95)
    ax.text(3.38, 6.13, f"SH 90.5%\n(inner ML {a2p_b}%)", ha="right", va="top", fontsize=4.8, fontweight="bold",
            color=PAL["arch2"], bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.90), zorder=5)
    branch(3.45, 6.25, 4.65, 6.78, PAL["arch2"], lw=1.05)
    badge(4.62, 6.88, f"ML {a2_b}% / SH {a2_sh:.1f}%", PAL["arch2"])
    wedge(4.65, 6.78, 0.85, 0.27, PAL["arch2"], "Arch 2: Split Motif IX (n = 3)",
          "M.EaeI, Thiotrichales, Cyanobacteriota", tc=PAL["arch2"], bold=True)
    branch(3.45, 6.25, 4.65, 5.75, PAL["plus1"], lw=0.95)
    wedge(4.65, 5.75, 0.85, 0.25, PAL["plus1"], "Arch 1 FPCQ sisters (n = 3)",
          "Chloroflexota, JACPUC01, phage AY129338", tc=PAL["ink"])

    # +1-less phage FPCQ (n = 5, ML 76% / SH 98.1%) vs Core SPCQ/TPCQ/Arch3 (n = 36)
    branch(2.25, 8.75, 3.05, 9.65, PAL["ink"], lw=0.95)
    branch(3.05, 9.65, 4.65, 7.82, PAL["ctrl"], lw=0.85)
    ax.text(4.62, 7.91, "ML 76% / SH 98.1%", ha="right", va="bottom", fontsize=4.8, fontweight="bold",
            color=PAL["ctrl"], bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.90), zorder=5)
    wedge(4.65, 7.82, 0.85, 0.24, PAL["ctrl"], "+1-less phage FPCQ clade (n = 5)", tc=PAL["muted"])

    branch(3.05, 9.65, 3.65, 10.25, PAL["plus1"], lw=1.0)
    badge(3.62, 10.36, "SH 82.9%")
    # Arch 1 SPCQ Bacillota (n = 6, ML 71% / SH 89.4%)
    branch(3.65, 10.25, 4.65, 8.72, PAL["plus1"], lw=0.9)
    ax.text(4.62, 8.81, "ML 71% / SH 89.4%", ha="right", va="bottom", fontsize=4.8, fontweight="bold",
            color=PAL["plus1"], bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.90), zorder=5)
    wedge(4.65, 8.72, 0.85, 0.25, PAL["plus1"], "Arch 1 SPCQ/SPCT: Bacillota (n = 6)", tc=PAL["ink"])

    # Long-branch outlier leaf A1_SPCT_Baci_10645 (n = 1)
    branch(3.65, 10.25, 4.10, 10.85, PAL["plus1"], lw=1.0)
    branch(4.10, 10.85, 4.95, 9.48, PAL["plus1"], lw=0.75)
    ax.scatter([4.95], [9.48], s=14, color=PAL["plus1"], zorder=4)
    ax.text(5.12, 9.48, "Arch 1 SPCT long-branch taxon (n = 1; bl = 3.77)",
            ha="left", va="center", fontsize=5.3, color=PAL["muted"])

    # 100%-overprinted TPCQ + Arch 3 Superclade (n = 29, ML 95% / SH 89.8%)
    branch(4.10, 10.85, 4.60, 11.00, PAL["plus1"], lw=1.15)
    badge(4.58, 11.12, f"ML {tp_b}% / SH {tp_sh:.1f}%", PAL["plus1"])
    branch(4.60, 11.00, 5.25, 10.35, PAL["plus1"], lw=1.0)
    wedge(5.25, 10.35, 0.95, 0.33, PAL["plus1"], "Arch 1 TPCQ superclade (n = 22)",
          "17 TPCQ (incl. P1 dmt, 5/5 UniProt) + 3 FPCQ + 2 SPCQ", tc=PAL["plus1"], bold=True)
    # DERIVED ARCH 3 nested inside Arch 1 TPCQ (n = 7, ML 97% / SH 99.5%)
    branch(4.60, 11.00, 5.25, 11.58, PAL["arch3"], lw=1.15)
    badge(5.23, 11.70, f"ML {a3_b}% / SH {a3_sh:.1f}%", PAL["arch3"])
    wedge(5.25, 11.58, 0.90, 0.30, PAL["arch3"], "Arch 3: Cis C-term fusions (n = 7)",
          "Single origin across 4 phyla (A0A6N9NLZ4)", tc=PAL["arch3"], bold=True)

    # === LOWER CLADE (n = 37, SH 95.5%) ===
    branch(1.5, 6.2, 2.25, 3.1, PAL["ink"], lw=1.0)
    badge(2.22, 3.20, "SH 95.5%")

    # Basal Arch 1 FPCQ Patescibacteriota (n = 3, ML 94% / SH 97.4%)
    branch(2.25, 3.1, 4.65, 1.58, PAL["plus1"], lw=0.9)
    badge(4.62, 1.68, "ML 94% / SH 97.4%", PAL["plus1"])
    wedge(4.65, 1.58, 0.85, 0.23, PAL["plus1"], "Arch 1 FPCQ: Patescibacteriota (n = 3)", tc=PAL["ink"])

    branch(2.25, 3.1, 3.10, 3.55, PAL["ink"], lw=0.95)
    badge(3.08, 3.65, "SH 89.2%")
    branch(3.10, 3.55, 4.65, 2.42, PAL["ctrl"], lw=0.85)
    badge(4.62, 2.52, "SH 95.6%", PAL["ctrl"])
    wedge(4.65, 2.42, 0.85, 0.23, PAL["ctrl"], "+1-less: M.SssI & phage clade (n = 5)", tc=PAL["muted"])

    branch(3.10, 3.55, 3.85, 4.02, PAL["ink"], lw=0.95)
    branch(3.85, 4.02, 4.65, 3.32, PAL["ctrl"], lw=0.9)
    badge(4.62, 3.42, "SH 91.7%", PAL["ctrl"])
    wedge(4.65, 3.32, 0.90, 0.28, PAL["ctrl"], "+1-less: classic MTases (n = 8)",
          "M.HhaI, E. coli Dcm, M.EcoRII, Arch. M.MthTI", tc=PAL["muted"])

    branch(3.85, 4.02, 4.65, 4.52, PAL["plus1"], lw=1.05)
    badge(4.62, 4.62, "SH 77.8%", PAL["plus1"])
    wedge(4.65, 4.52, 1.00, 0.36, PAL["plus1"], "Arch 1 FPCQ: Bacteria, phages & Archaea (n = 21)",
          "DPANN Archaea nested in Bacteria (subclade SH 95.0%)", tc=PAL["plus1"], bold=True)
