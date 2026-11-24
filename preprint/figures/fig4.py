"""Figure 4: AlphaFold 3 trans co-folding reproduces the cis Motif VIII/IX clamp across life."""
from fignimble import Figure
from style import STYLE

fig = Figure(width=183, height=130, style=STYLE)
fig.panel("panels/f4_concept.py", x=4, y=2, w=82, h=55)
fig.panel("panels/f4_benchmark.py", x=90, y=2, w=89, h=55)
fig.panel("panels/f4_3d_superposition.py", x=4, y=61, w=85, h=66)
fig.panel("panels/f4_lock_dist.py", x=93, y=61, w=86, h=66)
