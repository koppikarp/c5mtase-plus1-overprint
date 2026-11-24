"""Figure 2, Panel c: Independent gene-level paired INSIDE vs OUTSIDE excess 4D conservation."""
import json
import math
import numpy as np
from style import DATA_DIR, PAL


def _fmt_p(p):
    exp = int(math.floor(math.log10(p)))
    mant = p / (10 ** exp)
    return f"$P = {mant:.1f} \\times 10^{{{exp}}}$"


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "selection_4d_and_stop_stats.json").read_text())

    ax.axhline(0.0, ls="--", lw=0.65, color=PAL["muted"], zorder=0)

    groups = [
        ("fpcq_overprinted", 0.0, 1.0, PAL["plus1"], "FPCQ"),
        ("tspcq_overprinted", 2.3, 3.3, PAL["arch3"], "TPCQ/SPCQ"),
    ]
    rng = np.random.default_rng(7)
    for key, x0, x1, col, nm in groups:
        sub = d[key]
        yin = np.array(sub["gene_in_excess"])
        yout = np.array(sub["gene_out_excess"])
        n_g = sub["n_genes"]
        p_val = sub["gene_wilcoxon_p_greater"]
        title_lbl = f"{nm} (n = {n_g} loci)\n{_fmt_p(p_val)}"

        for a, b in zip(yin, yout):
            ax.plot([x0, x1], [a, b], lw=0.70, color=col, alpha=0.45, zorder=1)
        j0 = rng.normal(0, 0.035, size=len(yin))
        j1 = rng.normal(0, 0.035, size=len(yout))
        ax.scatter(x0 + j0, yin, s=16, color=col, alpha=0.88, edgecolors="white", linewidths=0.3, zorder=3)
        ax.scatter(x1 + j1, yout, s=16, color=PAL["ctrl"], alpha=0.88, edgecolors="white", linewidths=0.3, zorder=3)
        ax.plot([x0 - 0.18, x0 + 0.18], [np.median(yin)] * 2, lw=1.9, color=PAL["ink"], zorder=4)
        ax.plot([x1 - 0.18, x1 + 0.18], [np.median(yout)] * 2, lw=1.9, color=PAL["ink"], zorder=4)

        top_y = max(yin.max(), yout.max()) + 0.035
        ax.plot([x0, x0, x1, x1], [top_y - 0.015, top_y, top_y, top_y - 0.015], lw=0.7, color=PAL["ink"])
        ax.text((x0 + x1) / 2, top_y + 0.015, title_lbl, ha="center", va="bottom", fontsize=5.7, color=PAL["ink"])

    ax.set_xticks([0.0, 1.0, 2.3, 3.3])
    ax.set_xticklabels(["Inside\n+1", "Outside\n(IX–X)", "Inside\n+1", "Outside\n(IX–X)"], fontsize=5.9)
    ax.set_xlim(-0.55, 3.85)
    ax.set_ylim(-0.18, 0.48)
    ax.set_ylabel("Gene-level mean excess 4D", fontsize=6.5)
