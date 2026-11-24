"""Figure 2, Panel d: Non-circular homology-selected 5-frame stop-codon depletion."""
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

    case_meds = [np.median(d["cases"][f]) for f in frames]
    ctrl_meds = [np.median(d["ctrls"][f]) for f in frames]

    ax.bar(x - w / 2, case_meds, width=w, color=PAL["plus1"], alpha=0.85,
           label=f"Arch 1 +1 (n = {d['n_cases']}, no stop filter)")
    ax.bar(x + w / 2, ctrl_meds, width=w, color=PAL["ctrl"], alpha=0.75,
           label=f"1:1 GC3 ctrl (n = {d['n_ctrls']}, P = {d['gc3_balance_mwu_p']:.2f})")

    for i, f in enumerate(frames):
        vc = np.array(d["cases"][f])
        vt = np.array(d["ctrls"][f])
        ax.scatter(i - w / 2 + rng.normal(0, 0.04, len(vc)), vc, s=8, color=PAL["ink"], alpha=0.50, zorder=3)
        ax.scatter(i + w / 2 + rng.normal(0, 0.04, len(vt)), vt, s=8, color=PAL["ink"], alpha=0.35, zorder=3)

    p_cc = d["case_vs_ctrl_mwu"]["+1"]["mwu_p_twosided"]
    p_wg_max = max(v["p_less"] for v in d["within_gene_wilcoxon_plus1_lt_other"].values())

    ax.plot([-w / 2, -w / 2, w / 2, w / 2], [13.8, 14.6, 14.6, 13.8], lw=0.75, color=PAL["ink"])
    ax.text(0.0, 15.1, f"{case_meds[0]:.0f} vs. {ctrl_meds[0]:.0f} stops\n{_fmt_p(p_cc)}", ha="center", va="bottom",
            fontsize=5.6, fontweight="bold", color=PAL["plus1"])
    ax.text(2.55, 16.5, f"Within-gene +1 < all 4 frames\n(Wilcoxon max {_fmt_p(p_wg_max)})", ha="center", va="top",
            fontsize=5.3, color=PAL["muted"])

    ax.set_xticks(x)
    ax.set_xticklabels([f"Frame {f}" for f in frames], fontsize=6.0)
    ax.set_ylim(-0.6, 22.8)
    lo, hi = d["window_rel_pcq"]
    ax.set_ylabel(f"Stop codons in [{lo}..+{hi} aa]", fontsize=6.5)
    ax.legend(frameon=False, fontsize=5.2, loc="upper right", handlelength=1.0)
