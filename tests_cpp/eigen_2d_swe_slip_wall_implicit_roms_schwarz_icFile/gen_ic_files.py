import os

import numpy as np

from pschwarz.data_utils import load_meshes, decompose_domain_data, write_to_binary

ic_idx = 10

exe_dir = os.path.dirname(os.path.realpath(__file__))
order = os.path.basename(os.path.normpath(exe_dir))

# load data and mesh
data = np.loadtxt(f"../../../eigen_2d_swe_slip_wall_implicit/{order}/solution_full_gold.txt")
data = np.reshape(data, (30, 30, 3, -1), order="C")
data = np.transpose(data, (1, 0, 3, 2))
data_snap = data[:, :, ic_idx, :]
_, meshlist_decomp = load_meshes("./mesh_decomp", merge_decomp=False)

# decompose initial conditions
data_decomp = decompose_domain_data(
    data_snap,
    meshlist_decomp,
    6,
    is_ts=False,
    is_ts_decomp=False,
)

# write to binary
for j in range(2):
    for i in range(2):
        dom_idx = i + 2 * j
        data_out = np.transpose(data_decomp[i][j][0], (2, 0, 1))
        data_out = data_out.flatten(order="F")
        outfile = f"./ic_file_{dom_idx}.bin"
        write_to_binary(data_out, outfile)

print("Finished")


