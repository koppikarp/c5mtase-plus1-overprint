"""Figure 2: Orthogonal purifying selection signatures in the +1 reading frame."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=122, style=STYLE)
fig.panel("panels/f2_concept.py", x=4, y=2, w=76, h=56)
fig.panel("panels/f2_sliding.py", x=84, y=2, w=95, h=56)
fig.panel("panels/f2_paired.py", x=4, y=63, w=84, h=56)
fig.panel("panels/f2_stops.py", x=92, y=63, w=87, h=56)
