"""Figure 3, Panel c: Quantitative congruence between Frame-0 (Tree 1) and CSLSDILE/KYFLS (Tree 2) ML trees."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    d = json.loads((DATA_DIR / "phylogeny/ml_bootstrap_summary.json").read_text())
    pairs = d["patristic_pairs"]
    a3m = d["mantel_subchecks"]["arch3_only_zero_overlap"]
    nolb = d["mantel_subchecks"]["all_shared_no_longbranch_outlier"]

    order = [
        ("Arch1-Arch1", PAL["plus1"], "Arch 1 – Arch 1", 0.30, 9),
        ("Arch1-Arch3", PAL["frame0"], "Arch 1 – Arch 3", 0.38, 10),
        ("Arch2-any", PAL["arch2"], "Arch 2 pairs", 0.65, 13),
        ("Arch3-Arch3", PAL["arch3"], "Arch 3 – Arch 3 (0 overlap)", 0.92, 18),
    ]
    all_x, all_y = [], []
    # Plot non-outlier points cleanly + keep full scale readable
    for pt, col, lbl, alp, sz in order:
        xs = np.array([p["d_f0"] for p in pairs if p["pair_type"] == pt])
        ys = np.array([p["d_p1"] for p in pairs if p["pair_type"] == pt])
        all_x.extend(xs); all_y.extend(ys)
        ax.scatter(xs, ys, s=sz, color=col, alpha=alp, edgecolors="none", label=lbl, rasterized=False)

    all_x, all_y = np.array(all_x), np.array(all_y)
    # Fit trend line on main distribution (excluding extreme bl=3.77 outlier pairs for clean visual slope)
    mask_main = all_x < 4.5
    m, b = np.polyfit(all_x[mask_main], all_y[mask_main], 1)
    xg = np.linspace(0, min(max(all_x), 4.6), 50)
    ax.plot(xg, m * xg + b, ls="--", lw=1.0, color=PAL["ink"], zorder=5)

    p_all = d["mantel_p_5000perms"]
    p_a3 = a3m["mantel_p_pearson"]
    box_txt = (
        f"All {d['n_shared']} shared ({d['n_pairs']:,} pairs):\n"
        f"  r = {d['pearson_r']:.2f} (excl. 1 outlier: {nolb['pearson_r']:.2f}), P = {p_all:.4f}\n"
        f"Arch 3 only (0 DNA overlap, n = {a3m['n_taxa']}):\n"
        f"  r = {a3m['pearson_r']:.2f}, ρ = {a3m['spearman_rho']:.2f}, P = {p_a3:.4f}"
    )
    ax.text(0.04, 0.96, box_txt,
            transform=ax.transAxes, ha="left", va="top", fontsize=5.15, color=PAL["ink"],
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=PAL["ink"], lw=0.7))

    ax.set_xlabel("Frame-0 patristic distance (Tree 1)", fontsize=6.4)
    ax.set_ylabel("CSLSDILE/KYFLS distance (Tree 2)", fontsize=6.4)
    ax.set_xlim(-0.1, 5.2)
    ax.set_ylim(-0.1, max(all_y) * 1.08)
    ax.legend(frameon=False, fontsize=5.1, loc="lower right", markerscale=1.05, handletextpad=0.3, labelspacing=0.20)
