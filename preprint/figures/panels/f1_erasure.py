"""Figure 1, Panel d: Systematic database erasure of the +1 CSLSDILE/KYFLS family."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw(fig):
    ax = fig.add_subplot(111)
    u = json.loads((DATA_DIR / "uniref100_deletion_audit.json").read_text())
    n_aud = u["n_audited"]
    n_del = u["n_deleted"]
    n_act = u["n_active"]
    pct_del = u["pct_deleted"]
    pct_act = 100.0 - pct_del

    cats = [
        "PDB / CATH\nstructures",
        "Pfam / InterPro\ndomains",
        f"UniProtKB hits\n(n = {n_aud})",
    ]
    y = np.arange(len(cats))
    missing = np.array([100.0, 100.0, pct_del])
    present = np.array([0.0, 0.0, pct_act])

    ax.barh(y, present, height=0.52, color=PAL["frame0"], label="Retained")
    ax.barh(y, missing, left=present, height=0.52, color=PAL["plus1_lt"],
            edgecolor=PAL["plus1"], linewidth=0.75, label="Deleted / missing")

    ax.text(56, 2, f"{n_del}/{n_aud} deleted ({pct_del:.0f}%)", ha="center", va="center",
            fontsize=5.9, fontweight="bold", color=PAL["plus1"])
    ax.text(pct_act / 2.0, 2, str(n_act), ha="center", va="center", fontsize=5.9, fontweight="bold", color="white")
    ax.text(50, 1, f"{u['pfam_hits']} Pfam / {u['interpro_hits']} InterPro", ha="center", va="center",
            fontsize=5.9, fontweight="bold", color=PAL["plus1"])
    ax.text(50, 0, f"{u['pdb_cath_foldseek_hits_e_lt_1']} PDB / CATH (E < 1)", ha="center", va="center",
            fontsize=5.9, fontweight="bold", color=PAL["plus1"])

    ax.set_yticks(y)
    ax.set_yticklabels(cats, fontsize=5.6)
    ax.set_xlim(0, 103)
    ax.set_xlabel("Fraction (%)", fontsize=6.4)
    ax.set_title("Database blind spot", fontsize=6.7, fontweight="bold", pad=5)
    ax.legend(frameon=False, fontsize=5.5, loc="lower center", bbox_to_anchor=(0.45, -0.35), ncol=2, handlelength=1.0)
    fig.tight_layout(pad=0.35)
