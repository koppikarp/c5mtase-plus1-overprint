"""Figure 4, Panel d: Sub-Angstrom conservation of Motif IX R1 <-> KYFLS F/Y across all 3 architectures."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "structures/trans_vs_cis_3d_audit.json").read_text())
    rec = {r["name"]: r for r in d["trans_complexes"]}

    # Highlight the conserved 3.5-4.2 A lock zone
    ax.axhspan(3.4, 4.3, color=PAL["motif9_lt"], alpha=0.75, zorder=0)
    ax.axhline(4.02, ls="--", lw=0.75, color=PAL["arch3"], zorder=1)

    entries = [
        ("P1 157\n(A1)", rec["phagep1_arch1_dmt_trans_157aa"]["r1_kyf_dists_5"], PAL["plus1"], "trans"),
        ("P1 64\n(A1)", rec["phagep1_arch1_dmt_trans_64aa_core"]["r1_kyf_dists_5"], PAL["plus1"], "trans"),
        ("Archaea\n(A1)", rec["arch_arch1_nanob_trans"]["r1_kyf_dists_5"], PAL["plus1"], "trans"),
        ("A. soli\n(A1)", rec["bact_arch1_asoli_trans"]["r1_kyf_dists_5"], PAL["plus1"], "trans"),
        ("M.EaeI\n(A2)", rec["posctrl_arch2_eaei_trans"]["r1_kyf_dists_5"], PAL["arch2"], "intra-β"),
        ("A0A6N9\n(A3 cis)", [4.02], PAL["arch3"], "cis"),
        ("M.HhaI\n(neg)", rec["negctrl_hhai_asoli_trans"]["r1_kyf_dists_5"], PAL["ctrl"], "none"),
    ]

    rng = np.random.default_rng(19)
    x = np.arange(len(entries))
    for xi, (lbl, vals, col, mode) in zip(x, entries):
        arr = np.array(vals)
        m = arr.mean()
        ax.bar(xi, m, width=0.52, color=col, alpha=0.80, zorder=2)
        if len(arr) > 1:
            ax.scatter(xi + rng.uniform(-0.11, 0.11, len(arr)), arr, s=11, color=PAL["ink"], alpha=0.75, zorder=4)
        else:
            ax.scatter([xi], arr, s=22, marker="D", color=PAL["ink"], zorder=4)
        ax.text(xi, max(arr.max(), m) + 0.9, f"{m:.1f}" if m > 10 else f"{m:.2f}",
                ha="center", va="bottom", fontsize=5.6, fontweight="bold",
                color=PAL["ink"] if m < 5 else PAL["muted"])

    ax.text(2.5, 9.2, "Conserved 3.7–4.0 Å lock\n(trans → intra-β → cis)", ha="center", va="center",
            fontsize=5.8, fontweight="bold", color=PAL["motif9"],
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=PAL["motif9"], lw=0.65))
    ax.annotate("", xy=(2.5, 4.4), xytext=(2.5, 7.2),
                arrowprops=dict(arrowstyle="-|>", color=PAL["motif9"], lw=0.75, mutation_scale=7))

    ax.set_xticks(x)
    ax.set_xticklabels([e[0] for e in entries], fontsize=5.6)
    ax.set_ylim(0, 46)
    ax.set_ylabel("Motif IX R1 ↔ KYFLS F/Y Cα dist. (Å)", fontsize=6.3)
    ax.set_title("Conserved 3D lock across architectures", fontsize=6.6, fontweight="bold", pad=5)
