"""Figure 4, Panel b: All-5-model AF3 ipTM and min inter-chain PAE benchmark across 6 complexes."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "structures/trans_vs_cis_3d_audit.json").read_text())
    rec = {r["name"]: r for r in d["trans_complexes"]}

    order = [
        ("posctrl_arch2_eaei_trans", "Pos. ctrl: M.EaeI (Arch 2)", PAL["arch2"]),
        ("phagep1_arch1_dmt_trans_64aa_core", "Phage P1 dmt + 64aa core", PAL["plus1"]),
        ("phagep1_arch1_dmt_trans_157aa", "Phage P1 dmt + 157aa", PAL["plus1"]),
        ("bact_arch1_asoli_trans", "Bacteria: A. soli (Arch 1)", PAL["plus1"]),
        ("arch_arch1_nanob_trans", "Archaea: Nanobdellota", PAL["plus1"]),
        ("negctrl_hhai_asoli_trans", "Neg. ctrl: M.HhaI (+1-less)", PAL["ctrl"]),
    ]

    y = np.arange(len(order))[::-1]
    rng = np.random.default_rng(42)

    ax.axvline(0.50, ls=":", lw=0.75, color=PAL["muted"], zorder=0)
    for yi, (k, lbl, col) in zip(y, order):
        r = rec[k]
        vals = np.array(r["iptm_5"])
        mean = r["iptm_mean"]
        pae = r["pae_min_mean"]
        ax.barh(yi, mean, height=0.58, color=col, alpha=0.82, zorder=2)
        ax.scatter(vals, yi + rng.uniform(-0.12, 0.12, len(vals)),
                   s=12, color=PAL["ink"], alpha=0.75, zorder=4)
        pae_txt = f"{mean:.2f}  (PAE {pae:.1f} Å)"
        ax.text(max(mean, vals.max()) + 0.025, yi, pae_txt, va="center", ha="left",
                fontsize=5.6, fontweight="bold" if mean >= 0.6 else "normal",
                color=PAL["ink"] if mean >= 0.6 else PAL["muted"],
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.90), zorder=5)

    ax.set_yticks(y)
    ax.set_yticklabels([x[1] for x in order], fontsize=6.0)
    ax.set_xlim(0, 1.24)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Interface pTM (5 predicted models each)", fontsize=6.4)
    ax.set_title("Trans docking confidence (0/30 clashes)", fontsize=6.7, fontweight="bold", pad=5)
