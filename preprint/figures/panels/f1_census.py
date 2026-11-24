"""Figure 1, Panels b & c: Prokaryotic phyla (b) and INPHARED phage prevalence (c)."""
import json
import numpy as np
from style import DATA_DIR, PAL


def draw_gtdb(fig):
    ax = fig.add_subplot(111)
    gtdb = json.loads((DATA_DIR / "gtdb2339_c5mtase_plus1_census.json").read_text())
    pba = gtdb["summary"]["phylum_by_arch"]

    # Order phyla by total loci descending; show top 7 individual phyla (incl. both DPANN archaeal phyla!) + Other (6 phyla)
    ordered = sorted(pba.items(), key=lambda kv: (-kv[1]["total"], kv[0]))
    # Ensure both archaeal phyla (Nanobdellota, Aenigmatarchaeota) are shown explicitly!
    top_keys = [k for k, _ in ordered[:6]]
    if "p__Aenigmatarchaeota" not in top_keys:
        top_keys.append("p__Aenigmatarchaeota")
    other_keys = [k for k, _ in ordered if k not in top_keys]

    def fmt_ph(k):
        nm = k.replace("p__", "")
        if nm in ("Nanobdellota", "Aenigmatarchaeota"):
            return f"{nm} (Arc.)"
        return nm

    labels = [fmt_ph(k) for k in top_keys] + [f"Other ({len(other_keys)} phyla)"]

    def get_counts(keys):
        a1 = sum(pba[k]["Arch1_intact"] + pba[k]["Arch1_interrupted"] for k in keys)
        a2 = sum(pba[k]["Arch2_split_IX"] + pba[k]["Arch2_split_TRD_IX_X"] for k in keys)
        a3 = sum(pba[k]["Arch3_inframe_cterm"] for k in keys)
        return a1, a2, a3

    rows = [get_counts([k]) for k in top_keys] + [get_counts(other_keys)]
    v1 = np.array([r[0] for r in rows])
    v2 = np.array([r[1] for r in rows])
    v3 = np.array([r[2] for r in rows])

    y = np.arange(len(labels))[::-1]
    ax.barh(y, v1, height=0.62, color=PAL["plus1"], label="Arch 1 (+1 IX-indep.)")
    ax.barh(y, v2, left=v1, height=0.62, color=PAL["arch2"], label="Arch 2 / split")
    ax.barh(y, v3, left=v1 + v2, height=0.62, color=PAL["arch3"], label="Arch 3 (fusion)")
    for yi, t in zip(y, v1 + v2 + v3):
        ax.text(t + 0.35, yi, str(int(t)), va="center", ha="left", fontsize=5.9, color=PAL["ink"])
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.4)
    ax.set_xlabel(f"GTDB R232 loci (n = {gtdb['summary']['n_loci']}, {gtdb['summary']['n_unique_phyla']} phyla)", fontsize=6.1)
    ax.set_xlim(0, max(v1 + v2 + v3) + 3.4)
    ax.set_title("Prokaryotes (2,339 genomes)", fontsize=6.7, fontweight="bold", pad=5)
    ax.legend(frameon=False, fontsize=5.2, loc="lower right", handlelength=1.0, labelspacing=0.25)
    fig.tight_layout(pad=0.35)


def draw_phages(fig):
    ax = fig.add_subplot(111)
    inph = json.loads((DATA_DIR / "inphared_phage_c5mtase_census.json").read_text())
    m4_d = inph["motif4_denominators"]
    m4_o = inph["motif4_overprinted"]
    n_tot = inph["n_mtases"]
    n_ov = len(inph["overprinted"])
    n_other = sum(v for k, v in m4_d.items() if k not in ("TPCQ", "FPCQ", "SPCQ"))
    n_intact = inph["refined_arch_counts"].get("Arch1_intact", 0)
    n_interr = inph["refined_arch_counts"].get("Arch1_interrupted", 0)

    specs = [
        ("TPCQ", m4_o.get("TPCQ", 0), m4_d.get("TPCQ", 0), PAL["plus1"]),
        ("FPCQ", m4_o.get("FPCQ", 0), m4_d.get("FPCQ", 0), PAL["plus1"]),
        ("All", n_ov, n_tot, PAL["frame0"]),
        ("SPCQ", m4_o.get("SPCQ", 0), m4_d.get("SPCQ", 0), PAL["ctrl"]),
        ("Other", 0, n_other, PAL["ctrl"]),
    ]
    cats = [f"{nm}\n({o}/{t})" for nm, o, t, _ in specs]
    pcts = [100.0 * o / max(1, t) for _, o, t, _ in specs]
    cols = [c for _, _, _, c in specs]

    x = np.arange(len(cats))
    ax.bar(x, pcts, width=0.58, color=cols, edgecolor="none")
    for xi, pv in zip(x, pcts):
        ax.text(xi, pv + 2.5, f"{pv:.1f}%" if pv > 0 else "0%",
                ha="center", va="bottom", fontsize=5.8, fontweight="bold", color=PAL["ink"])
    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=5.4)
    ax.set_ylabel("+1 Arch 1 prevalence (%)", fontsize=6.4)
    ax.set_ylim(0, 106)
    ax.set_title(f"INPHARED phages ({inph['n_phages']:,})", fontsize=6.7, fontweight="bold", pad=5)
    ax.text(0.97, 0.90, f"{n_intact}/{n_ov} intact IX–X\n({n_interr} split; 0 Arch 2/3)", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.3, color=PAL["plus1"],
            bbox=dict(boxstyle="round,pad=0.25", fc=PAL["plus1_lt"], ec="none"))
