"""Figure 5, Panel a: Bacteriophage P1 dmt locus map and synonymous-F0 / nonsense-+1 validation mutations."""
from matplotlib.patches import FancyBboxPatch
from style import PAL


def _bar(ax, x, y, w, h, fc, ec=None, lw=0.55, r=0.08, zorder=3):
    ec = ec if ec is not None else PAL["ink"]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder))


def draw(fig):
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 17.6)
    ax.set_ylim(0, 5.8)
    ax.axis("off")

    ax.text(8.8, 5.55, "Bacteriophage P1 dmt locus (AF234172: 66,412–68,676, + strand, 754 aa)",
            ha="center", va="center", fontsize=6.8, fontweight="bold", color=PAL["ink"])

    # Scale mapping: 0..754 aa mapped onto x = 0.5..17.1 (dx = 16.6)
    def sx(aa):
        return 0.5 + (aa / 754.0) * 16.6

    # Frame 0 dmt bar: C5-MTase (1..473) + Dam N6-adenine MTase (485..754)
    y0 = 4.35
    _bar(ax, sx(0), y0, sx(473) - sx(0), 0.52, fc=PAL["frame0"], ec=PAL["ink"], lw=0.55)
    _bar(ax, sx(480), y0, sx(754) - sx(480), 0.52, fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=0.7)
    ax.plot([sx(473), sx(480)], [y0 + 0.26, y0 + 0.26], lw=1.0, color=PAL["ink"])

    # Conserved C5-MTase motifs on P1 dmt
    for s, e, lbl, col in [(75, 84, "IV", PAL["motif4"]), (123, 134, "VI", "#365685"),
                           (201, 218, "VIII", "#365685"), (385, 405, "IX", PAL["motif9"]),
                           (445, 463, "X", "#365685")]:
        _bar(ax, sx(s), y0, max(0.36, sx(e) - sx(s)), 0.52, fc=col, ec="none", r=0.04, zorder=4)
        ax.text(sx(s) + max(0.36, sx(e) - sx(s)) / 2, y0 + 0.26, lbl,
                ha="center", va="center", fontsize=5.4, fontweight="bold", color="white", zorder=5)

    ax.text(sx(300), y0 + 0.26, "TRD", ha="center", va="center", fontsize=5.6, color="white", zorder=5)
    ax.text(sx(617), y0 + 0.26, "C-terminal Dam (N6-adenine MTase, 485–754 aa)",
            ha="center", va="center", fontsize=5.8, fontweight="bold", color=PAL["frame0"], zorder=5)
    ax.text(sx(235), y0 + 0.72, "N-terminal C5-MTase domain (1–473 aa; Motif IX R385 intact)",
            ha="center", va="center", fontsize=5.8, color=PAL["ink"])

    # Embedded +1 ORF underneath (aa 100..263 in dmt coordinates: nt 299..790)
    y1 = 3.40
    # Full 157aa SD+ATG isoform (aa 106..263) in light terracotta, 64aa SD+ATG core (aa 199..263) in bold terracotta!
    _bar(ax, sx(106), y1, sx(263) - sx(106), 0.46, fc=PAL["plus1_lt"], ec=PAL["plus1"], lw=0.75, zorder=3)
    _bar(ax, sx(199), y1, sx(263) - sx(199), 0.46, fc=PAL["plus1"], ec=PAL["ink"], lw=0.6, zorder=4)
    ax.text(sx(152), y1 + 0.23, "+1 (157 aa)", ha="center", va="center", fontsize=5.4, color=PAL["plus1"], fontweight="bold", zorder=5)
    ax.text(sx(231), y1 + 0.23, "64aa core", ha="center", va="center", fontsize=5.4, color="white", fontweight="bold", zorder=5)

    # SD+ATG flags
    ax.annotate("GGAG + ATG\n(66,728; 157 aa)", xy=(sx(106), y1), xytext=(sx(45), 2.52),
                fontsize=5.4, color=PAL["plus1"], ha="center", va="top",
                arrowprops=dict(arrowstyle="-|>", color=PAL["plus1"], lw=0.65, mutation_scale=6))
    ax.annotate("GGAG + ATG\n(67,007; 64 aa, ipTM 0.83)", xy=(sx(199), y1), xytext=(sx(235), 2.52),
                fontsize=5.5, fontweight="bold", color=PAL["plus1"], ha="center", va="top",
                arrowprops=dict(arrowstyle="-|>", color=PAL["plus1"], lw=0.75, mutation_scale=6))

    # Zoomed codon-register box on bottom-right showing exact synonymous Frame-0 / nonsense +1 switches!
    _bar(ax, 7.15, 0.18, 10.20, 2.82, fc="white", ec=PAL["ink"], lw=0.8, r=0.18, zorder=1)
    ax.text(12.25, 2.72, "Synonymous Frame-0 / Nonsense +1 validation switches",
            ha="center", va="center", fontsize=6.0, fontweight="bold", color=PAL["ink"])

    muts = [
        ("c.642T>A\n(nt 67,053)", "F0: V214V\nGTT → GTA (silent)", "+1: L115*\nTTG → TAG (in CLLSQVLE)"),
        ("c.645C>A\n(nt 67,056)", "F0: V215V\nGTC → GTA (silent)", "+1: S116*\nTCG → TAG (in CLLSQVLE)"),
    ]
    for k, (m_id, f0_s, p1_s) in enumerate(muts):
        ry = 1.92 - k * 1.10
        _bar(ax, 7.40, ry - 0.42, 2.15, 0.84, fc="white", ec=PAL["ink"], lw=0.55, r=0.10, zorder=2)
        ax.text(8.47, ry, m_id, ha="center", va="center", fontsize=5.5, fontweight="bold", color=PAL["ink"], zorder=3, linespacing=1.25)
        _bar(ax, 9.75, ry - 0.42, 3.45, 0.84, fc=PAL["frame0_lt"], ec=PAL["frame0"], lw=0.55, r=0.10, zorder=2)
        ax.text(11.47, ry, f0_s, ha="center", va="center", fontsize=5.4, color=PAL["frame0"], fontweight="bold", zorder=3, linespacing=1.25)
        _bar(ax, 13.40, ry - 0.42, 3.70, 0.84, fc=PAL["plus1_lt"], ec=PAL["plus1"], lw=0.55, r=0.10, zorder=2)
        ax.text(15.25, ry, p1_s, ha="center", va="center", fontsize=5.4, color=PAL["plus1"], fontweight="bold", zorder=3, linespacing=1.25)

    # Dashed connector from 64aa core to zoom box
    ax.plot([sx(215), 7.15], [y1, 2.65], ls=":", lw=0.75, color=PAL["plus1"])
