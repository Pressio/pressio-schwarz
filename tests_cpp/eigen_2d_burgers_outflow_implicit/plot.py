
from pschwarz.vis_utils import plot_contours

# ----- START USER INPUTS -----

# 0: x-velocity
# 1: y-velocity
varplot = 1

# ----- END USER INPUTS -----

if varplot == 0:
    varlabel = r"X-velocity"
    nlevels = 21
    skiplevels = 2
    contourbounds = [0.0, 0.5]
elif varplot == 1:
    varlabel = r"Y-velocity"
    nlevels = 21
    skiplevels = 2
    contourbounds = [0.0, 0.5]

plot_contours(
    varplot,
    meshdirs="./",
    datadirs="./",
    nvars=2,
    dataroot="burgers_outflow2d_solution",
    nlevels=nlevels,
    skiplevels=skiplevels,
    contourbounds=contourbounds,
    plotskip=2,
    varlabel=varlabel,
)
