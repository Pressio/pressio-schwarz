import os

import numpy as np

from pschwarz.prom_utils import gen_pod_bases


exe_dir = os.path.dirname(os.path.realpath(__file__))
order = os.path.basename(os.path.normpath(exe_dir))

data = np.loadtxt(f"../../../eigen_2d_swe_slip_wall_implicit/{order}/solution_full_gold.txt")
data = np.reshape(data, (30, 30, 3, -1), order="C")
data = np.transpose(data, (1, 0, 3, 2))

gen_pod_bases(
    outdir="./trial_space",
    meshdir="./full_mesh",
    datalist=[data],
    nvars=3,
    dataroot="swe_slipWall2d_solution",
    center_method="init_cond",
    norm_method="one",
)
