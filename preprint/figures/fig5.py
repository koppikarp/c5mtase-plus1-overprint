"""Figure 5: Bacteriophage P1 dmt locus, synonymous validation switches, and compact microprotein core."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=106, style=STYLE)
fig.panel("panels/f5_p1_locus.py", x=4, y=2, w=175, h=52)
fig.panel("panels/f5_plddt.py", x=4, y=56, w=175, h=48)
