"""Figure 2, Panel a: Conceptual schematic of the non-stop-risk 4D (A/G/P/T) dual-coding test."""
from matplotlib.patches import FancyBboxPatch
from style import PAL


def _box(ax, x, y, w, h, fc, ec="none", lw=0.6, r=0.10, zorder=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder))


def draw(fig):
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 10.0)
    ax.set_ylim(0, 6.4)
    ax.axis("off")

    ax.text(5.0, 6.12, "4D synonymous test (A, G, P, T only)", ha="center", va="center",
            fontsize=6.8, fontweight="bold", color=PAL["ink"])

    # Real DNA sequence: G C T | A G C | G G A  -> Frame 0: Ala(4D), Ser, Gly(4D)
    # Frame +1 starting at pos 1: C T A (Leu) | G C G (Ala) -> wobble base 'T' of Ala is pos 2 of +1 Leu!
    bases = ["G", "C", "T", "A", "G", "C", "G", "G", "A"]
    x0 = 1.05
    bw = 0.84
    for i, b in enumerate(bases):
        bx = x0 + i * (bw + 0.04)
        is_wobble = (i in (2, 8))
        fc = PAL["plus1_lt"] if is_wobble else "white"
        ec = PAL["plus1"] if is_wobble else PAL["ink"]
        lw = 1.05 if is_wobble else 0.65
        _box(ax, bx, 3.85, bw, 0.78, fc=fc, ec=ec, lw=lw, r=0.08, zorder=3)
        ax.text(bx + bw / 2, 4.24, b, ha="center", va="center",
                fontsize=7.2, fontweight="bold", color=PAL["plus1"] if is_wobble else PAL["ink"], zorder=4)

    # Frame 0 brackets above
    for k, (lbl, is4d) in enumerate([("Ala (4D)", True), ("Ser", False), ("Gly (4D)", True)]):
        lx = x0 + (k * 3) * (bw + 0.04)
        rx = lx + 3 * bw + 2 * 0.04
        ax.plot([lx + 0.06, lx + 0.06, rx - 0.06, rx - 0.06], [4.72, 4.92, 4.92, 4.72],
                lw=0.85, color=PAL["frame0"])
        ax.text((lx + rx) / 2, 5.22, lbl, ha="center", va="center",
                fontsize=6.0, fontweight="bold" if is4d else "normal", color=PAL["frame0"])
    ax.text(0.12, 5.05, "F0", ha="left", va="center", fontsize=6.5, fontweight="bold", color=PAL["frame0"])

    # Frame +1 brackets below (CTA = Leu, GCG = Ala)
    for k, lbl in enumerate(["+1 Leu (pos 2 = F0 4D)", "+1 Ala"]):
        lx = x0 + (1 + k * 3) * (bw + 0.04)
        rx = lx + 3 * bw + 2 * 0.04
        ax.plot([lx + 0.06, lx + 0.06, rx - 0.06, rx - 0.06], [3.76, 3.55, 3.55, 3.76],
                lw=0.85, color=PAL["plus1"])
        ax.text((lx + rx) / 2, 3.24, lbl, ha="center", va="center",
                fontsize=5.8, fontweight="bold", color=PAL["plus1"])
    ax.text(0.12, 3.45, "+1", ha="left", va="center", fontsize=6.5, fontweight="bold", color=PAL["plus1"])

    # Two concise contrast boxes at bottom (clean white background, crisp borders)
    _box(ax, 0.25, 0.25, 4.55, 2.35, fc="white", ec=PAL["plus1"], lw=0.9, r=0.16)
    ax.text(2.52, 2.22, "Inside +1 gene", ha="center", va="center", fontsize=6.3, fontweight="bold", color=PAL["plus1"])
    ax.text(2.52, 1.25, "F0 3rd base = +1 2nd base\n(missense in +1)\n4D bases locked across\nphyla; stops purged",
            ha="center", va="center", fontsize=5.8, color=PAL["ink"], linespacing=1.35)

    _box(ax, 5.20, 0.25, 4.55, 2.35, fc="white", ec=PAL["ink"], lw=0.85, r=0.16)
    ax.text(7.48, 2.22, "Outside +1 / controls", ha="center", va="center", fontsize=6.3, fontweight="bold", color=PAL["ink"])
    ax.text(7.48, 1.25, "F0 3rd base freely drifts\n(neutral baseline)\nExcess 4D ≈ 0;\n~8 stops in +1 frame",
            ha="center", va="center", fontsize=5.8, color=PAL["ink"], linespacing=1.35)
