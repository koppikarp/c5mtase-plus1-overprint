"""Graphical Abstract (Nature double-column width, 183 x 78 mm)."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=78, style=STYLE)
fig.panel("panels/ga_panel.py", x=0, y=0, w=183, h=78, label=None)
