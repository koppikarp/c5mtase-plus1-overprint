"""Figure 4, Panel a: Conceptual schematic of the trans-vs-cis 3D clamp structure prediction test."""
from matplotlib.patches import FancyBboxPatch
from style import PAL


def _box(ax, x, y, w, h, fc, ec="none", lw=0.65, r=0.16, zorder=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder))


def draw(fig):
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 10.0)
    ax.set_ylim(0, 6.2)
    ax.axis("off")

    ax.text(5.0, 5.92, "Trans structure prediction test", ha="center", va="center",
            fontsize=6.7, fontweight="bold", color=PAL["ink"])

    # Step 1: Unconnected inputs
    _box(ax, 0.20, 3.65, 3.85, 1.85, fc="white", ec=PAL["ink"], lw=0.8)
    ax.text(2.12, 5.20, "1. Unlinked chains", ha="center", va="center", fontsize=6.1, fontweight="bold", color=PAL["ink"])
    _box(ax, 0.48, 4.42, 3.30, 0.44, fc=PAL["frame0"], ec="none", r=0.09)
    _box(ax, 2.75, 4.42, 0.55, 0.44, fc=PAL["motif9"], ec="none", r=0.05, zorder=3)
    ax.text(1.60, 4.64, "Chain A (Motifs I–X)", ha="center", va="center", fontsize=5.4, color="white", zorder=4)
    ax.text(3.02, 4.64, "IX", ha="center", va="center", fontsize=5.5, fontweight="bold", color="white", zorder=4)
    _box(ax, 0.48, 3.84, 3.30, 0.42, fc=PAL["plus1"], ec="none", r=0.09)
    ax.text(2.13, 4.05, "Chain B (+1 CSLS / KYFLS)", ha="center", va="center", fontsize=5.4, fontweight="bold", color="white", zorder=4)

    # Arrow right
    ax.annotate("", xy=(5.00, 4.55), xytext=(4.15, 4.55),
                arrowprops=dict(arrowstyle="-|>", color=PAL["ink"], lw=0.9, mutation_scale=8))
    ax.text(4.58, 4.92, "Co-fold\n5 models", ha="center", va="center", fontsize=5.4, color=PAL["muted"])

    # Step 2: Compare trans vs cis + DNA
    _box(ax, 5.10, 3.65, 4.70, 1.85, fc="white", ec=PAL["ink"], lw=0.8)
    ax.text(7.45, 5.20, "2. Superpose onto cis & DNA", ha="center", va="center", fontsize=6.1, fontweight="bold", color=PAL["ink"])
    _box(ax, 5.42, 3.88, 1.45, 0.98, fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=0.8, r=0.22)
    ax.text(6.14, 4.37, "MTase\n(Kabsch)", ha="center", va="center", fontsize=5.5, color=PAL["frame0"], fontweight="bold")
    _box(ax, 7.05, 3.98, 0.72, 0.78, fc=PAL["motif9"], ec="none", r=0.14)
    ax.text(7.41, 4.37, "R1", ha="center", va="center", fontsize=5.9, fontweight="bold", color="white")
    _box(ax, 8.18, 3.98, 1.32, 0.78, fc=PAL["plus1"], ec=PAL["arch3"], lw=1.0, r=0.16)
    ax.text(8.84, 4.37, "trans vs\ncis F/Y", ha="center", va="center", fontsize=5.5, fontweight="bold", color="white")
    ax.plot([7.77, 8.18], [4.37, 4.37], ls="--", lw=1.0, color=PAL["ink"])

    # Bottom panel: Control design strip
    _box(ax, 0.20, 0.35, 9.60, 2.85, fc="white", ec=PAL["ink"], lw=0.8)
    ax.text(5.0, 2.88, "3. Built-in controls across 6 complexes (30 models)", ha="center", va="center",
            fontsize=6.1, fontweight="bold", color=PAL["ink"])

    items = [
        (0.45, PAL["arch2"], "Positive control", "M.EaeI α + β\n(known dimer)"),
        (3.60, PAL["plus1"], "Arch 1 tests (×4)", "Phage P1, Bacteria,\nDPANN Archaea"),
        (6.85, PAL["ctrl"], "Negative control", "+1-less M.HhaI\n+ A. soli +1"),
    ]
    for ix, col, h1, h2 in items:
        _box(ax, ix, 0.62, 2.72, 1.92, fc="white", ec=col, lw=0.95, r=0.16, zorder=3)
        ax.text(ix + 1.36, 2.15, h1, ha="center", va="center", fontsize=5.9, fontweight="bold", color=col, zorder=4)
        ax.text(ix + 1.36, 1.28, h2, ha="center", va="center", fontsize=5.6, color=PAL["ink"], linespacing=1.3, zorder=4)
