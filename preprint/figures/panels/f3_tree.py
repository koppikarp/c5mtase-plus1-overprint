"""Figure 3, Panel b: DNMT-rooted IQ-TREE 2 (Q.pfam+R6) ML Phylogeny (Tree 1, 97 unique taxa collapsed to exact clades)."""
import json
from matplotlib.patches import Polygon
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "phylogeny/ml_bootstrap_summary.json").read_text())
    ks = d["key_clade_supports"]
    n_t1 = d["n_taxa_tree1"]
    a3_sh, a3_uf = ks["Arch3_monophyly_tree1"]["sh_alrt"], int(round(ks["Arch3_monophyly_tree1"]["ufboot"]))
    tp_sh, tp_uf = ks["TPCQ_Arch3_superclade_tree1"]["sh_alrt"], int(round(ks["TPCQ_Arch3_superclade_tree1"]["ufboot"]))
    a2_sh, a2_uf = ks["Arch2_monophyly_tree1"]["sh_alrt"], int(round(ks["Arch2_monophyly_tree1"]["ufboot"]))
    a2p_sh, a2p_uf = ks["Arch2_parent_clade_tree1"]["sh_alrt"], int(round(ks["Arch2_parent_clade_tree1"]["ufboot"]))
    a2o_sh, a2o_uf = ks["Arch2_outer_clade_tree1"]["sh_alrt"], int(round(ks["Arch2_outer_clade_tree1"]["ufboot"]))
    sp_sh, sp_uf = ks["Bacillota_SPCQ_SplitTRD_clade_tree1"]["sh_alrt"], int(round(ks["Bacillota_SPCQ_SplitTRD_clade_tree1"]["ufboot"]))

    ax.set_xlim(-0.3, 14.9)
    ax.set_ylim(-0.6, 12.8)
    ax.axis("off")

    ax.text(6.8, 12.52, f"Frame-0 C5-MTase IQ-TREE 2 phylogeny ({n_t1} taxa, Q.pfam+R6, SH-aLRT / UFBoot2)",
            ha="center", va="center", fontsize=6.4, fontweight="bold", color=PAL["ink"])

    # Exact DNMT-rooted IQ-TREE 2 partition of all 97 unique taxa:
    # Outgroup (2) + Basal [11 phage + 1 SPCT LB] (12) + Lower [3 Pate + 8 classic + 25 FPCQ/DPANN] (36)
    # + Upper [2 outer + 2 A1 sisters + 3 Arch2 + 5 phage + 6 Bacillota + 22 TPCQ + 7 Arch3] (47) = 97!
    def branch(x0, y0, x1, y1, col=PAL["ink"], lw=0.9):
        ax.plot([x0, x0, x1], [y0, y1, y1], color=col, lw=lw, solid_capstyle="round", zorder=2)

    def wedge(x, y, dx, dy, fc, label, sub="", tc=PAL["ink"], bold=False):
        tri = Polygon([(x, y), (x + dx, y + dy), (x + dx, y - dy)],
                      closed=True, facecolor=fc, edgecolor=PAL["ink"], linewidth=0.5, zorder=3)
        ax.add_patch(tri)
        ax.text(x + dx + 0.16, y + (0.12 if sub else 0.0), label,
                ha="left", va="center", fontsize=5.75, fontweight="bold" if bold else "normal", color=tc)
        if sub:
            ax.text(x + dx + 0.16, y - 0.23, sub, ha="left", va="center", fontsize=5.1, color=PAL["muted"])

    def badge(x, y, txt, col=PAL["ink"], va="bottom"):
        ax.text(x, y, txt, ha="right", va=va, fontsize=4.85, fontweight="bold", color=col,
                bbox=dict(boxstyle="round,pad=0.09", fc="white", ec="none", alpha=0.88), zorder=5)

    # Root at x=0.2, y=3.6
    branch(0.2, 3.6, 4.6, -0.15, PAL["ctrl"], lw=0.95)
    wedge(4.6, -0.15, 0.75, 0.20, PAL["ctrl_lt"], "Outgroup: Human DNMT1, TRDMT1 (n = 2)", tc=PAL["muted"])

    # Ingroup (n=95, 90.8/91)
    branch(0.2, 3.6, 0.9, 5.1, PAL["ink"], lw=0.95)
    badge(0.83, 5.20, "91 / 91")

    # Basal 12-taxon clade (72.1/73): 11 +1-less phages (90.6/81) + 1 long-branch SPCT (bl=3.46)
    branch(0.9, 5.1, 1.9, 1.05, PAL["ctrl"], lw=0.9)
    badge(1.84, 1.13, "72 / 73", PAL["ctrl"])
    branch(1.9, 1.05, 4.6, 0.65, PAL["ctrl"], lw=0.85)
    badge(4.56, 0.73, "91 / 81", PAL["ctrl"])
    wedge(4.6, 0.65, 0.90, 0.25, PAL["ctrl"], "+1-less phage outgroup (n = 11)", tc=PAL["muted"])
    branch(1.9, 1.05, 4.95, 1.45, PAL["plus1"], lw=0.75)
    ax.scatter([4.95], [1.45], s=13, color=PAL["plus1"], zorder=4)
    ax.text(5.12, 1.45, "Arch 1 SPCT long-branch taxon (n = 1; bl = 3.46)",
            ha="left", va="center", fontsize=5.2, color=PAL["muted"])

    # Main radiation (n=83) -> Upper Clade (n=47) & Lower Clade (n=36)
    branch(0.9, 5.1, 1.5, 6.2, PAL["ink"], lw=0.95)

    # === LOWER CLADE (n = 36, SH 95.2 / UF 68) ===
    branch(1.5, 6.2, 2.25, 3.45, PAL["ink"], lw=1.0)
    badge(2.18, 3.33, "95 / 68", va="top")
    # Basal Arch 1 FPCQ Patescibacteriota (n = 3, 98.6/98)
    branch(2.25, 3.45, 4.65, 2.30, PAL["plus1"], lw=0.9)
    badge(4.62, 2.39, "99 / 98", PAL["plus1"])
    wedge(4.65, 2.30, 0.85, 0.22, PAL["plus1"], "Arch 1 FPCQ: Patescibacteriota (n = 3)", tc=PAL["ink"])

    # Inner 33-taxon clade (77.6/48) -> 8 +1-less (92.1/69) & 25 Arch 1 FPCQ + DPANN (89.0/66)
    branch(2.25, 3.45, 3.15, 3.92, PAL["ink"], lw=0.95)
    branch(3.15, 3.92, 4.65, 3.22, PAL["ctrl"], lw=0.85)
    badge(4.62, 3.31, "92 / 69", PAL["ctrl"])
    wedge(4.65, 3.22, 0.90, 0.26, PAL["ctrl"], "+1-less: M.SssI, M.MthTI & phage (n = 8)", tc=PAL["muted"])

    branch(3.15, 3.92, 4.65, 4.38, PAL["plus1"], lw=1.05)
    badge(4.62, 4.48, "89 / 66", PAL["plus1"])
    wedge(4.65, 4.38, 1.00, 0.36, PAL["plus1"], "Arch 1 FPCQ: Bacteria, phages & Archaea (n = 25)",
          "21 Arch 1 (DPANN in Bacteria: 81/72) + 4 +1-less (HhaI/Dcm)", tc=PAL["plus1"], bold=True)

    # === UPPER CLADE (n = 47, SH 100 / UF 93) ===
    branch(1.5, 6.2, 2.25, 8.85, PAL["ink"], lw=1.0)
    badge(2.18, 8.96, "100 / 93", va="bottom")

    # Arch 2 + Arch 1 FPCQ clade (n = 7, SH 92.2 / UF 85)
    branch(2.25, 8.85, 3.10, 6.50, PAL["plus1"], lw=0.95)
    badge(3.06, 6.40, f"{a2o_sh:.0f} / {a2o_uf}", PAL["arch2"], va="top")
    # Outer sister pair (n = 2: phage AY129338 + phage JN564907)
    branch(3.10, 6.50, 4.65, 5.55, PAL["ctrl"], lw=0.8)
    wedge(4.65, 5.55, 0.80, 0.21, PAL["ctrl_lt"], "Phage AY129338 (A1) + JN564907 (n = 2)", tc=PAL["muted"])
    # Inner 5-taxon clade (88.3 / 69) -> Arch 2 (3) + Arch 1 FPCQ sisters (2)
    branch(3.10, 6.50, 3.75, 6.95, PAL["arch2"], lw=1.0)
    badge(3.72, 7.05, f"{a2p_sh:.0f} / {a2p_uf}", PAL["arch2"])
    branch(3.75, 6.95, 4.65, 6.42, PAL["plus1"], lw=0.9)
    badge(4.62, 6.50, "96 / 90", PAL["plus1"])
    wedge(4.65, 6.42, 0.85, 0.22, PAL["plus1"], "Arch 1 FPCQ sisters: Chlo., JACPUC01 (n = 2)", tc=PAL["ink"])
    branch(3.75, 6.95, 4.65, 7.45, PAL["arch2"], lw=1.05)
    badge(4.62, 7.54, f"{a2_sh:.0f} / {a2_uf}", PAL["arch2"])
    wedge(4.65, 7.45, 0.85, 0.27, PAL["arch2"], "Arch 2: Split Motif IX (n = 3)",
          "M.EaeI, Thiotrichales, Cyanobacteriota (AU basal P = 0.022)", tc=PAL["arch2"], bold=True)

    # SPCQ / TPCQ / Arch 3 arm (n = 40) -> 5 +1-less phages (96.4/99) & 35 +1-carrying (75.1/60)
    branch(2.25, 8.85, 3.05, 9.75, PAL["ink"], lw=0.95)
    branch(3.05, 9.75, 4.65, 8.48, PAL["ctrl"], lw=0.85)
    badge(4.62, 8.56, "96 / 99", PAL["ctrl"])
    wedge(4.65, 8.48, 0.85, 0.23, PAL["ctrl"], "+1-less phage FPCQ clade (n = 5)", tc=PAL["muted"])

    branch(3.05, 9.75, 3.70, 10.35, PAL["plus1"], lw=1.0)
    badge(3.67, 10.45, "75 / 60", PAL["plus1"])

    # Bacillota SPCQ/SPCT + Split-TRD/IX/X (n = 6, 78.7/75)
    branch(3.70, 10.35, 4.65, 9.42, PAL["plus1"], lw=0.95)
    badge(4.62, 9.51, f"{sp_sh:.0f} / {sp_uf}", PAL["plus1"])
    wedge(4.65, 9.42, 0.85, 0.26, PAL["plus1"], "Arch 1 SPCQ/SPCT: Bacillota (n = 6)",
          "5 intact Arch 1 + 1 split-TRD/IX/X (GCA_035934735.1)", tc=PAL["ink"])

    # 100%-overprinted 29-taxon TPCQ + Arch 3 Superclade (SH 99.5 / UF 100!)
    branch(3.70, 10.35, 4.55, 11.00, PAL["plus1"], lw=1.15)
    badge(4.52, 11.12, f"{tp_sh:.1f} / {tp_uf}", PAL["plus1"])
    branch(4.55, 11.00, 5.25, 10.42, PAL["plus1"], lw=1.0)
    wedge(5.25, 10.42, 0.95, 0.31, PAL["plus1"], "Arch 1 TPCQ superclade (n = 22)",
          "17 TPCQ (incl. P1 dmt, 5/5 UniProt) + 3 FPCQ + 2 SPCQ", tc=PAL["plus1"], bold=True)

    # DERIVED ARCH 3 nested inside Arch 1 TPCQ (n = 7, SH 93.5 / UF 98; drop-LB 98.7 / 100!)
    branch(4.55, 11.00, 5.25, 11.58, PAL["arch3"], lw=1.15)
    badge(5.22, 11.70, f"{a3_sh:.1f} / {a3_uf}", PAL["arch3"])
    wedge(5.25, 11.58, 0.90, 0.30, PAL["arch3"], "Arch 3: Cis C-term fusions (n = 7)",
          "Single origin across 4 phyla (AU basal P = 0.005)", tc=PAL["arch3"], bold=True)
