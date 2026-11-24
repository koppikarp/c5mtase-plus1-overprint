"""Figure 2, Panel b: Sliding-window composition-corrected excess 4D conservation (whole-gene EXP_4D, ±1 genome-bootstrap s.e.)."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "selection_4d_and_stop_stats.json").read_text())

    ax.axvspan(-30, 120, color=PAL["plus1_lt"], alpha=0.55, zorder=0)
    ax.axvspan(75, 118, color=PAL["plus1"], alpha=0.16, zorder=1)
    ax.axhline(0.0, ls="--", lw=0.7, color=PAL["muted"], zorder=1)
    ax.axvline(0.0, ls=":", lw=0.7, color=PAL["motif4"], zorder=1)
    ax.axvline(120.0, ls=":", lw=0.7, color=PAL["plus1"], zorder=1)

    for key, prefix, col, mk, lw in [
        ("fpcq_overprinted", "FPCQ", PAL["plus1"], "o", 1.5),
        ("tspcq_overprinted", "TPCQ/SPCQ", PAL["arch3"], "s", 1.35),
        ("ctrl_no_plus1", "Controls", PAL["ctrl"], "^", 1.15),
    ]:
        sub = d[key]
        wins = sub["windows"]
        xs = np.array([w["bin_mid"] for w in wins])
        ys = np.array([w["excess_4d"] for w in wins])
        se = np.array([w["se_boot"] for w in wins])
        lbl = f"{prefix} ({sub['n_genes']} loci)"
        ax.fill_between(xs, ys - se, ys + se, color=col, alpha=0.14, lw=0, zorder=2)
        ax.plot(xs, ys, marker=mk, ms=3.4, lw=lw, color=col, label=lbl, zorder=3)

    ax.text(3, 0.51, "Motif IV\n(PCQ)", fontsize=5.5, color=PAL["motif4"], ha="left", va="top", fontweight="bold")
    ax.text(90, 0.53, "CSLSDILE\n/ KYFLS", fontsize=5.4, color=PAL["plus1"], ha="center", va="top", fontweight="bold")
    ax.text(124, 0.30, "+1 stop\ncliff", fontsize=5.4, color=PAL["plus1"], ha="left", va="top")
    ax.text(-26, -0.15, "+1 ORF span", fontsize=5.2, color=PAL["plus1"], ha="left", va="bottom", style="italic")

    ax.set_xlim(-55, 215)
    ax.set_ylim(-0.17, 0.58)
    ax.set_xlabel("Position relative to Motif IV PCQ (aa)", fontsize=6.5)
    ax.set_ylabel("Excess 4D identity (OBS − EXP)", fontsize=6.5)
    ax.legend(frameon=True, facecolor="white", framealpha=0.90, edgecolor="none",
              fontsize=5.3, loc="upper right", handlelength=1.3, labelspacing=0.22)
