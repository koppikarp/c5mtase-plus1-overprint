"""Figure 2, Panel d: Zero-confounder 5-frame stop-codon depletion vs. within-gene flanking-codon synonymous null."""
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
    d = json.loads((DATA_DIR / "selection_4d_and_stop_stats.json").read_text())["stop_purge_5frames"]

    frames = ["+1", "+2", "-0", "-1", "-2"]
    x = np.arange(len(frames))
    w = 0.34
    rng = np.random.default_rng(11)

    obs_meds = [np.median(d["cases"][f]) for f in frames]
    null_meds = [np.median(d["flanking_null"][f]) for f in frames]

    ax.bar(x - w / 2, obs_meds, width=w, color=PAL["plus1"], alpha=0.88,
           label=f"Arch 1 observed (n = {d['n_cases']}, no filter)")
    ax.bar(x + w / 2, null_meds, width=w, color=PAL["ctrl"], alpha=0.72,
           label="Within-gene flanking-codon null")

    for i, f in enumerate(frames):
        vc = np.array(d["cases"][f])
        vn = np.array(d["flanking_null"][f])
        ax.scatter(i - w / 2 + rng.normal(0, 0.04, len(vc)), vc, s=6.5, color=PAL["ink"], alpha=0.45, zorder=3)
        ax.scatter(i + w / 2 + rng.normal(0, 0.04, len(vn)), vn, s=6.5, color=PAL["ink"], alpha=0.28, zorder=3)

    r_p1 = d["obs_vs_flanking_null_wilcoxon"]["+1"]
    p_null_p1 = r_p1["wilcoxon_p_less"]
    p_wg_max = max(v["p_less"] for v in d["within_gene_wilcoxon_plus1_lt_other"].values())

    ax.plot([-w / 2, -w / 2, w / 2, w / 2], [11.6, 12.4, 12.4, 11.6], lw=0.75, color=PAL["ink"])
    ax.text(0.05, 12.9, f"{obs_meds[0]:.0f} vs. {null_meds[0]:.1f} stops\n{_fmt_p(p_null_p1)} ({r_p1['n_obs_lt_null']}/{d['n_cases']})",
            ha="center", va="bottom", fontsize=5.35, fontweight="bold", color=PAL["plus1"])
    ax.text(2.0, 10.6, "TTA/CTA\nshadow", ha="center", va="bottom", fontsize=4.85,
            color=PAL["muted"], style="italic")
    ax.text(2.70, 18.6, f"Within-gene +1 < all 4 frames\n(exact Wilcoxon max {_fmt_p(p_wg_max)})",
            ha="center", va="top", fontsize=5.15, color=PAL["muted"])

    ax.set_xticks(x)
    ax.set_xticklabels([f"Frame {f}" for f in frames], fontsize=6.0)
    ax.set_ylim(-0.6, 23.5)
    lo, hi = d["window_rel_pcq"]
    ax.set_ylabel(f"Stop codons in [{lo}..+{hi} aa]", fontsize=6.5)
    ax.legend(frameon=False, fontsize=5.1, loc="upper right", handlelength=1.0)
