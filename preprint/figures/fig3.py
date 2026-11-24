"""Figure 3: DNMT-rooted LG+Gamma ML phylogeny and Frame-0 vs CSLSDILE/KYFLS congruence."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=128, style=STYLE)
fig.panel("panels/f3_concept.py", x=3, y=2, w=76, h=55)
fig.panel("panels/f3_congruence.py", x=3, y=62, w=76, h=63)
fig.panel("panels/f3_tree.py", x=82, y=2, w=98, h=123)
