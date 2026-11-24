"""Shared Nature journal styling and matte semantic palette for all figures.

Edit colors, fonts, or sizes here and every figure in the project updates together.
Compatible with Fignimble (pachterlab/Fignimble) CLI and macOS visual editor.
"""
from pathlib import Path
from fignimble import Labels, Style

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR.parent.parent / "data"

# Harmonious matte semantic palette (consistent across Graphical Abstract + Figures 1-5)
PAL = {
    "frame0": "#4A6FA5",       # Matte slate blue: Frame-0 C5-MTase catalytic core (Motifs I-VIII, X)
    "frame0_lt": "#DCE6F2",
    "plus1": "#C85A3F",        # Matte terracotta: +1 overprinted CSLSDILE/KYFLS gene/domain (Arch 1)
    "plus1_lt": "#F6DFD9",
    "motif9": "#4E8A6E",       # Matte sage green: Conserved Motif IX (and R1 lock residue)
    "motif9_lt": "#D9ECE3",
    "motif4": "#C9933B",       # Matte ochre: Conserved Motif IV (PCQ catalytic Cys)
    "motif4_lt": "#F7EBD5",
    "arch2": "#8E4B78",        # Matte dusty plum: Derived Arch 2 split-Motif-IX (M.EaeI / M.EcoHK31I)
    "arch2_lt": "#EFE2EB",
    "arch3": "#22859A",        # Matte petrol teal: Derived Arch 3 in-frame C-terminal fusion
    "arch3_lt": "#D7EFF4",
    "ctrl": "#7E848C",         # Matte warm gray: +1-less conventional C5-MTase controls (M.HhaI/Dcm)
    "ctrl_lt": "#EAECEF",
    "dna": "#5A524C",          # Matte charcoal-sand: DNA duplex / backbone
    "dna_lt": "#E6E1DC",
    "ink": "#111111",          # Primary black ink
    "muted": "#5C6168",        # Secondary gray text
    "grid": "#E3E6EB",         # Subtle axis gridline
    "bg_box": "#FFFFFF",       # Pure white panel inset fill (black rounded border)
}

STYLE = Style.journal(
    "nature",
    font="Arial",
    size=7.0,
    linewidth=0.65,
    labels=Labels.journal("nature"),
    palette=[
        PAL["frame0"],
        PAL["plus1"],
        PAL["motif9"],
        PAL["arch3"],
        PAL["arch2"],
        PAL["motif4"],
        PAL["ctrl"],
    ],
    despine=True,
    rc={
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    },
)
