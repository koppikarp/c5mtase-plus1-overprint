"""Figure 1: Three genomic architectures, pan-prokaryotic/viral census, and database erasure."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=124, style=STYLE)
fig.panel("panels/f1_architectures.py", x=4, y=2, w=175, h=60)
fig.panel("panels/f1_census.py", function="draw_gtdb", x=4, y=67, w=70, h=53)
fig.panel("panels/f1_census.py", function="draw_phages", x=77, y=67, w=46, h=53)
fig.panel("panels/f1_erasure.py", x=126, y=67, w=53, h=53)
