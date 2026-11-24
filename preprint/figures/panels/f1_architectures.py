"""Figure 1, Panel a: Conceptual schematic of the three C5-MTase architectures vs conventional enzymes."""
from matplotlib.patches import FancyBboxPatch
from style import PAL


def _bar(ax, x, y, w, h, fc, ec=None, lw=0.55, r=0.08, zorder=3):
    ec = ec if ec is not None else PAL["ink"]
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
    ax.add_patch(p)


def draw(fig):
    ax = fig.add_subplot(111)
    ax.set_xlim(0.0, 17.6)
    ax.set_ylim(0.15, 6.25)
    ax.axis("off")

    # Header coordinate guides
    ax.text(5.70, 6.02, "Catalytic core", ha="center", va="center", fontsize=6.4, fontweight="bold", color=PAL["muted"])
    ax.text(9.85, 6.02, "TRD", ha="center", va="center", fontsize=6.4, fontweight="bold", color=PAL["muted"])
    ax.text(12.25, 6.02, "C-terminal collar", ha="center", va="center", fontsize=6.4, fontweight="bold", color=PAL["muted"])

    for gx in (4.92, 8.15, 11.45, 13.00):
        ax.plot([gx, gx], [0.25, 5.72], ls=":", lw=0.55, color="#D5D9E0", zorder=0)

    def _motifs(y, core_fc, include_ix=True):
        specs = [(3.78, 0.48, "I", core_fc), (4.65, 0.56, "IV", PAL["motif4"]),
                 (6.18, 0.52, "VI", core_fc), (7.86, 0.62, "VIII", core_fc),
                 (12.72, 0.56, "X", core_fc)]
        if include_ix:
            specs.append((11.12, 0.66, "IX", PAL["motif9"]))
        for mx, mw, ml, mc in specs:
            _bar(ax, mx, y, mw, 0.44, fc=mc, ec="none", r=0.04, zorder=4)
            ax.text(mx + mw / 2, y + 0.22, ml, ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=5)

    # Row 1: Conventional
    y = 5.02
    ax.text(0.05, y + 0.26, "Conventional", ha="left", va="center", fontsize=7.0, fontweight="bold", color=PAL["ink"])
    ax.text(0.05, y - 0.08, "e.g., M.HhaI, Dcm", ha="left", va="center", fontsize=5.8, style="italic", color=PAL["muted"])
    _bar(ax, 3.45, y, 10.20, 0.44, fc=PAL["ctrl_lt"], ec=PAL["ctrl"], lw=0.7)
    _motifs(y, PAL["ctrl"], include_ix=True)
    ax.text(14.00, y + 0.22, "No +1 gene (stop-rich)", ha="left", va="center", fontsize=6.0, color=PAL["muted"])

    # Row 2: Arch 1
    y = 3.60
    ax.text(0.05, y + 0.26, "Intact dual-coding", ha="left", va="center", fontsize=7.0, fontweight="bold", color=PAL["plus1"])
    ax.text(0.05, y - 0.14, "Arch 1 (37 GTDB +\n39 phage, incl. P1 dmt)", ha="left", va="center", fontsize=5.7, color=PAL["muted"])
    _bar(ax, 3.45, y, 10.20, 0.44, fc=PAL["frame0"], ec=PAL["ink"], lw=0.55)
    _motifs(y, "#365685", include_ix=True)
    _bar(ax, 4.65, y - 0.52, 3.92, 0.42, fc=PAL["plus1"], ec=PAL["ink"], lw=0.55, zorder=4)
    ax.text(6.61, y - 0.31, "+1  CSLSDILE / KYFLS", ha="center", va="center", fontsize=6.0, fontweight="bold", color="white", zorder=5)
    ax.text(14.00, y + 0.02, "Motif IX intact in F0;\n+1 stops before TRD", ha="left", va="center", fontsize=6.0, color=PAL["ink"])

    # Row 3: Arch 2
    y = 1.98
    ax.text(0.05, y + 0.26, "Split-subunit", ha="left", va="center", fontsize=7.0, fontweight="bold", color=PAL["arch2"])
    ax.text(0.05, y - 0.14, "Arch 2 (M.EcoHK31I,\nM.EaeI; 2 + 2 split)", ha="left", va="center", fontsize=5.7, color=PAL["muted"])
    _bar(ax, 3.45, y, 10.20, 0.44, fc=PAL["frame0"], ec=PAL["ink"], lw=0.55)
    _motifs(y, "#365685", include_ix=False)
    _bar(ax, 11.12, y + 0.02, 0.66, 0.40, fc="white", ec=PAL["arch2"], lw=0.8, r=0.04, zorder=4)
    ax.text(11.45, y + 0.22, "ΔIX", ha="center", va="center", fontsize=5.6, fontweight="bold", color=PAL["arch2"], zorder=5)
    _bar(ax, 5.05, y - 0.52, 6.73, 0.42, fc=PAL["arch2"], ec=PAL["ink"], lw=0.55, zorder=4)
    _bar(ax, 5.05, y - 0.52, 3.52, 0.42, fc=PAL["plus1"], ec="none", r=0.06, zorder=5)
    _bar(ax, 11.12, y - 0.52, 0.66, 0.42, fc=PAL["motif9"], ec="none", r=0.05, zorder=5)
    ax.text(6.81, y - 0.31, "+1 β  CSLSDILE / ESFLN", ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=6)
    ax.text(11.45, y - 0.31, "IX", ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=6)
    ax.text(14.00, y + 0.02, "Motif IX transferred\nonto +1 β subunit", ha="left", va="center", fontsize=6.0, color=PAL["arch2"])

    # Row 4: Arch 3
    y = 0.44
    ax.text(0.05, y + 0.28, "C-terminal fusion", ha="left", va="center", fontsize=7.0, fontweight="bold", color=PAL["arch3"])
    ax.text(0.05, y - 0.12, "Arch 3 (A0A6N9NLZ4;\n4 GTDB + UniProt)", ha="left", va="center", fontsize=5.7, color=PAL["muted"])
    _bar(ax, 3.45, y, 10.20, 0.44, fc=PAL["frame0"], ec=PAL["ink"], lw=0.55)
    _motifs(y, "#365685", include_ix=True)
    _bar(ax, 13.70, y, 2.75, 0.44, fc=PAL["arch3"], ec=PAL["ink"], lw=0.55, zorder=4)
    ax.text(15.08, y + 0.22, "CSLS / KYFLS (cis)", ha="center", va="center", fontsize=5.8, fontweight="bold", color="white", zorder=5)
    ax.text(16.60, y + 0.22, "0 overlap", ha="left", va="center", fontsize=5.7, color=PAL["arch3"])
