"""Figure 5, Panel b: Chain B per-residue AF3 pLDDT aligned at CSLSDILE, showing compact ordered core."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "structures/trans_vs_cis_3d_audit.json").read_text())
    rec = {r["name"]: r for r in d["trans_complexes"]}

    # Align every Chain B residue coordinate relative to start of CSLSDILE = 0!
    ax.axvspan(0, 20, color=PAL["plus1_lt"], alpha=0.65, zorder=0)
    ax.axvline(0, ls=":", lw=0.65, color=PAL["plus1"], zorder=1)
    ax.axvline(-13, ls="--", lw=0.75, color=PAL["ink"], zorder=1)
    ax.axhline(70, ls=":", lw=0.6, color=PAL["muted"], zorder=0)

    specs = [
        ("phagep1_arch1_dmt_trans_64aa_core", "P1 64aa core (ipTM 0.83)", PAL["plus1"], 1.85, "-"),
        ("phagep1_arch1_dmt_trans_157aa", "P1 157aa full (ipTM 0.65)", PAL["frame0"], 1.25, "-"),
        ("bact_arch1_asoli_trans", "A. soli 153aa (ipTM 0.64)", PAL["motif9"], 1.15, "--"),
        ("negctrl_hhai_asoli_trans", "M.HhaI neg. ctrl (ipTM 0.16)", PAL["ctrl"], 1.10, ":"),
    ]
    for k, lbl, col, lw, ls in specs:
        r = rec[k]
        pl = np.array(r["plddt_B_model0"])
        csl0 = r["csl_pos_B"] - 1
        xs = np.arange(len(pl)) - csl0
        ax.plot(xs, pl, lw=lw, ls=ls, color=col, label=lbl, zorder=3)

    ax.text(-16, 104, "P1 64aa\nSD+ATG", ha="right", va="top", fontsize=5.4, fontweight="bold", color=PAL["ink"])
    ax.text(10, 104, "CSLSDILE..KYFLS", ha="center", va="top", fontsize=5.4, fontweight="bold", color=PAL["plus1"])
    ax.text(-112, 71.5, "pLDDT 70", ha="left", va="bottom", fontsize=5.2, color=PAL["muted"])
    ax.text(-94, 86, "Low-pLDDT N-term\n(ΔN41–46 in M.EcoHK31I)", ha="left", va="center",
            fontsize=5.2, style="italic", color=PAL["muted"],
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.8))

    ax.set_xlim(-115, 55)
    ax.set_ylim(18, 108)
    ax.set_xlabel("Chain B position relative to CSLSDILE start (aa)", fontsize=6.4)
    ax.set_ylabel("Per-residue pLDDT", fontsize=6.4)
    ax.legend(frameon=True, facecolor="white", framealpha=0.90, edgecolor="none",
              fontsize=5.3, loc="lower right", handlelength=1.4, labelspacing=0.22)
