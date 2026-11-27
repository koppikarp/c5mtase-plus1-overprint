"""Graphical Abstract panel: evolutionary and structural unification of the three C5-MTase architectures."""
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from style import PAL


def _box(ax, x, y, w, h, fc, ec="none", lw=0.7, r=0.18, alpha=1.0, zorder=2):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={r}",
        facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(p)
    return p


def draw(fig):
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 18.3)
    ax.set_ylim(0, 7.8)
    ax.axis("off")

    # Clean white rounded cards with crisp black borders
    _box(ax, 0.20, 0.30, 7.65, 7.20, fc="white", ec=PAL["ink"], lw=0.85, r=0.28, zorder=0)
    _box(ax, 9.85, 4.10, 8.25, 3.40, fc="white", ec=PAL["ink"], lw=0.85, r=0.28, zorder=0)
    _box(ax, 9.85, 0.30, 8.25, 3.40, fc="white", ec=PAL["ink"], lw=0.85, r=0.28, zorder=0)

    # LEFT CARD: Ancestral intact dual-coding (Arch 1)
    ax.text(4.02, 7.12, "Ancestral intact dual-coding (Arch 1)", ha="center", va="center",
            fontsize=8.0, fontweight="bold", color=PAL["ink"])
    ax.text(4.02, 6.68, "13 prokaryotic phyla (Bacteria, Archaea) + 17% of phage C5-MTases",
            ha="center", va="center", fontsize=6.1, color=PAL["muted"])

    # Frame 0 gene bar
    _box(ax, 0.58, 5.52, 6.90, 0.60, fc=PAL["frame0"], ec=PAL["ink"], lw=0.55, r=0.10, zorder=3)
    _box(ax, 1.80, 5.52, 0.55, 0.60, fc=PAL["motif4"], ec="none", r=0.04, zorder=4)
    _box(ax, 5.38, 5.52, 0.72, 0.60, fc=PAL["motif9"], ec="none", r=0.04, zorder=4)
    ax.text(1.18, 5.82, "Frame 0", ha="center", va="center", fontsize=6.0, fontweight="bold", color="white", zorder=5)
    ax.text(2.07, 5.82, "IV", ha="center", va="center", fontsize=6.0, fontweight="bold", color="white", zorder=5)
    ax.text(3.75, 5.82, "Catalytic core (I–VIII)", ha="center", va="center", fontsize=6.0, color="white", zorder=5)
    ax.text(5.74, 5.82, "IX", ha="center", va="center", fontsize=6.2, fontweight="bold", color="white", zorder=5)
    ax.text(6.82, 5.82, "X", ha="center", va="center", fontsize=6.0, fontweight="bold", color="white", zorder=5)

    # +1 gene bar
    _box(ax, 1.90, 4.72, 2.55, 0.54, fc=PAL["plus1"], ec=PAL["ink"], lw=0.55, r=0.10, zorder=3)
    ax.text(3.18, 4.99, "+1  CSLSDILE / KYFLS", ha="center", va="center", fontsize=6.0, fontweight="bold", color="white", zorder=5)
    ax.text(5.74, 4.99, "Stops before TRD & IX", ha="center", va="center", fontsize=5.8, color=PAL["muted"])

    # Arrow down to trans assembly
    ax.annotate("", xy=(4.02, 3.88), xytext=(4.02, 4.56),
                arrowprops=dict(arrowstyle="-|>", color=PAL["ink"], lw=0.85, mutation_scale=8))
    ax.text(4.20, 4.22, "trans assembly", ha="left", va="center", fontsize=6.0, style="italic", color=PAL["muted"])

    # Trans heterodimer cartoon
    _box(ax, 0.95, 1.25, 3.30, 2.18, fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=1.1, r=0.42, zorder=2)
    ax.text(2.28, 2.75, "C5-MTase", ha="center", va="center", fontsize=7.0, fontweight="bold", color=PAL["frame0"])
    ax.text(2.28, 2.34, "(intact Motifs I–X)", ha="center", va="center", fontsize=5.8, color=PAL["muted"])
    _box(ax, 1.22, 1.54, 0.96, 0.42, fc=PAL["motif4"], ec="none", r=0.14, zorder=3)
    ax.text(1.70, 1.75, "Cys (IV)", ha="center", va="center", fontsize=5.6, fontweight="bold", color="white", zorder=4)

    _box(ax, 3.48, 1.58, 1.05, 1.38, fc=PAL["motif9"], ec=PAL["ink"], lw=0.6, r=0.22, zorder=4)
    ax.text(4.00, 2.44, "Motif IX", ha="center", va="center", fontsize=5.9, fontweight="bold", color="white", zorder=5)
    ax.text(4.00, 1.96, "R1", ha="center", va="center", fontsize=6.4, fontweight="bold", color="white", zorder=5)

    _box(ax, 5.15, 1.58, 2.02, 1.38, fc=PAL["plus1"], ec=PAL["ink"], lw=0.6, r=0.28, zorder=4)
    ax.text(6.16, 2.46, "+1 Microprotein", ha="center", va="center", fontsize=6.0, fontweight="bold", color="white", zorder=5)
    ax.text(6.16, 1.98, "KYFLS (F/Y)", ha="center", va="center", fontsize=6.0, color="white", zorder=5)

    ax.plot([4.53, 5.15], [2.15, 2.15], ls="--", lw=1.2, color=PAL["ink"], zorder=6)
    ax.scatter([4.53, 5.15], [2.15, 2.15], s=14, color=PAL["ink"], zorder=7)
    ax.text(4.84, 3.22, "3.7–4.0 Å", ha="center", va="center", fontsize=6.1, fontweight="bold", color=PAL["ink"])
    ax.text(4.02, 0.72, "Phage P1 dmt, Bacteria & Archaea (AF3 ipTM 0.63–0.83)",
            ha="center", va="center", fontsize=6.0, color=PAL["ink"])

    # CENTER EVOLUTIONARY TRANSITION ARROWS (clean clearance!)
    arr_up = FancyArrowPatch((7.96, 4.75), (9.74, 5.75), connectionstyle="arc3,rad=-0.16",
                             arrowstyle="-|>", mutation_scale=11, lw=1.35, color=PAL["arch2"])
    arr_dn = FancyArrowPatch((7.96, 3.05), (9.74, 2.05), connectionstyle="arc3,rad=0.16",
                             arrowstyle="-|>", mutation_scale=11, lw=1.35, color=PAL["arch3"])
    ax.add_patch(arr_up)
    ax.add_patch(arr_dn)
    ax.text(8.78, 6.22, "Motif IX\nhandoff", ha="center", va="center", fontsize=6.2, fontweight="bold", color=PAL["arch2"])
    ax.text(8.88, 4.70, "1 event\n(FPCQ)", ha="center", va="center", fontsize=5.7, color=PAL["muted"])
    ax.text(8.88, 3.08, "1 event\n(TPCQ)", ha="center", va="center", fontsize=5.7, color=PAL["muted"])
    ax.text(8.78, 1.56, "C-term\nfusion", ha="center", va="center", fontsize=6.2, fontweight="bold", color=PAL["arch3"])

    # TOP-RIGHT CARD: Derived split-subunit (Arch 2)
    ax.text(13.98, 7.12, "Derived split-subunit (Arch 2: M.EcoHK31I, M.EaeI)",
            ha="center", va="center", fontsize=7.3, fontweight="bold", color=PAL["ink"])
    _box(ax, 10.22, 6.12, 4.05, 0.50, fc=PAL["frame0"], ec=PAL["ink"], lw=0.5, r=0.09, zorder=3)
    ax.text(12.24, 6.37, "α subunit (Motifs I–VIII, X; no IX)", ha="center", va="center", fontsize=5.7, color="white", zorder=4)
    _box(ax, 11.45, 5.44, 2.20, 0.48, fc=PAL["plus1"], ec=PAL["ink"], lw=0.5, r=0.09, zorder=3)
    _box(ax, 13.60, 5.44, 0.85, 0.48, fc=PAL["motif9"], ec=PAL["ink"], lw=0.5, r=0.09, zorder=4)
    ax.text(12.52, 5.68, "+1 β (CSLSDILE)", ha="center", va="center", fontsize=5.6, fontweight="bold", color="white", zorder=5)
    ax.text(14.02, 5.68, "IX", ha="center", va="center", fontsize=5.9, fontweight="bold", color="white", zorder=5)
    ax.text(12.35, 4.66, "Trans lock becomes intra-β", ha="center", va="center", fontsize=6.0, style="italic", color=PAL["arch2"])

    _box(ax, 15.05, 4.50, 2.75, 2.15, fc=PAL["arch2_lt"], ec=PAL["arch2"], lw=0.95, r=0.26, zorder=2)
    _box(ax, 15.30, 5.50, 1.00, 0.82, fc=PAL["plus1"], ec="none", r=0.14, zorder=3)
    _box(ax, 16.55, 5.50, 1.00, 0.82, fc=PAL["motif9"], ec="none", r=0.14, zorder=3)
    ax.plot([16.30, 16.55], [5.91, 5.91], lw=1.2, color=PAL["ink"], zorder=4)
    ax.text(15.80, 5.91, "F78", ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=4)
    ax.text(17.05, 5.91, "R152", ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=4)
    ax.text(16.42, 4.94, "Intra-β lock\n3.76 Å (ipTM 0.84)", ha="center", va="center", fontsize=5.8, color=PAL["ink"], zorder=4)

    # BOTTOM-RIGHT CARD: Derived C-terminal fusion (Arch 3)
    ax.text(13.98, 3.32, "Derived C-terminal fusion (Arch 3: de-overprinted)",
            ha="center", va="center", fontsize=7.3, fontweight="bold", color=PAL["ink"])
    _box(ax, 10.22, 2.12, 2.85, 0.54, fc=PAL["frame0"], ec=PAL["ink"], lw=0.5, r=0.09, zorder=3)
    _box(ax, 12.32, 2.12, 0.52, 0.54, fc=PAL["motif9"], ec="none", r=0.04, zorder=4)
    _box(ax, 13.07, 2.12, 1.50, 0.54, fc=PAL["arch3"], ec=PAL["ink"], lw=0.5, r=0.09, zorder=3)
    ax.text(11.27, 2.39, "Single frame (I–X)", ha="center", va="center", fontsize=5.7, color="white", zorder=5)
    ax.text(12.58, 2.39, "IX", ha="center", va="center", fontsize=5.7, fontweight="bold", color="white", zorder=5)
    ax.text(13.82, 2.39, "CSLS / KYFLS", ha="center", va="center", fontsize=5.6, fontweight="bold", color="white", zorder=5)
    ax.text(12.35, 1.24, "Dual-coding conflict resolved;\ntrans lock becomes cis (4 phyla)",
            ha="center", va="center", fontsize=5.9, style="italic", color=PAL["arch3"])

    _box(ax, 15.05, 0.70, 2.75, 2.15, fc=PAL["arch3_lt"], ec=PAL["arch3"], lw=0.95, r=0.26, zorder=2)
    _box(ax, 15.30, 1.68, 1.00, 0.82, fc=PAL["motif9"], ec="none", r=0.14, zorder=3)
    _box(ax, 16.55, 1.68, 1.00, 0.82, fc=PAL["arch3"], ec="none", r=0.14, zorder=3)
    ax.plot([16.30, 16.55], [2.09, 2.09], lw=1.2, color=PAL["ink"], zorder=4)
    ax.text(15.80, 2.09, "R290", ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=4)
    ax.text(17.05, 2.09, "Y411", ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=4)
    ax.text(16.42, 1.12, "Cis monomer lock\n4.02 Å (UFBoot 98%)", ha="center", va="center", fontsize=5.8, color=PAL["ink"], zorder=4)
