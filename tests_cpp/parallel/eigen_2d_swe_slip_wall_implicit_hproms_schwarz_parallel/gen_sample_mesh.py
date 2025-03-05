import os
from argparse import ArgumentParser

from pschwarz.samp_utils import gen_sample_mesh

parser = ArgumentParser()
parser.add_argument("--outdir", dest="outdir")
args = parser.parse_args()

gen_sample_mesh(
    "random",
    os.path.join(args.outdir, "full_mesh_decomp"),
    0.2,
    os.path.join(args.outdir, "sample_mesh_decomp"),
    seed_qdeim=False,
    seed_phys_bounds=True,
    seed_dom_bounds=True,
    samp_phys_bounds=False,
    samp_dom_bounds=False,
    randseed=2,
)

