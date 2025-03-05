import os

from pschwarz.vis_utils import plot_contours

# ----- START USER INPUTS -----

# 0: x-velocity
# 1: y-velocity
varplot = 0

# ----- END USER INPUTS -----

exe_dir = os.path.dirname(os.path.realpath(__file__))
order = os.path.basename(os.path.normpath(exe_dir))

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

# TODO: modify monolithic directory to correct stencil order
plot_contours(
    varplot,
    meshdirs=[f"../../eigen_2d_burgers_outflow_implicit/{order}/", "./mesh",],
    datadirs=[f"../../eigen_2d_burgers_outflow_implicit/{order}/", "./"],
    nvars=2,
    dataroot="burgers_outflow2d_solution",
    plotlabels=["Monolithic", "Schwarz 2x2"],
    nlevels=nlevels,
    skiplevels=skiplevels,
    contourbounds=contourbounds,
    plotskip=2,
    varlabel=varlabel,
    plotbounds=True,
    bound_colors=["b", "r", "m", "c"],
    figdim_base=[8, 9],
    vertical=False,
)
