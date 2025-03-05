
from pschwarz.samp_utils import gen_sample_mesh

gen_sample_mesh(
    "random",
    "./full_mesh",
    0.2,
    "./sample_mesh",
    seed_qdeim=False,
    seed_phys_bounds=False,
    samp_phys_bounds=False,
    randseed=2,
)
