
from pschwarz.samp_utils import gen_sample_mesh

gen_sample_mesh(
    "random",
    "./full_mesh_decomp",
    0.2,
    "./sample_mesh_decomp",
    seed_qdeim=False,
    seed_phys_bounds=False,
    seed_dom_bounds=False,
    samp_phys_bounds=False,
    samp_dom_bounds=False,
    randseed=2,
)

