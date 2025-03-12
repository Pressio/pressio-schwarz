import os
from argparse import ArgumentParser, RawTextHelpFormatter

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors

from pschwarz.data_utils import load_info_domain, load_mesh_single


def make_subdomain_list(ndom_list, default=None):
    return [[[default for _ in range(ndom_list[2])] for _ in range(ndom_list[1])] for _ in range(ndom_list[0])]


def main(decompdir, sampdir, outdir=None, outstr=""):

    colors = ["r", "b", "c", "m"]

    # read decomposition info
    ndom_list, overlap = load_info_domain(decompdir)
    ndomains = np.prod(ndom_list)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    xmin_full = np.inf
    xmax_full = -np.inf
    ymin_full = np.inf
    ymax_full = -np.inf

    coords_full_sub = make_subdomain_list(ndom_list)
    connect_full_sub = make_subdomain_list(ndom_list)
    dcell_sub = make_subdomain_list(ndom_list)
    connect_stencil_sub = make_subdomain_list(ndom_list)
    gids_samp_sub = make_subdomain_list(ndom_list)
    gids_stencil_sub = make_subdomain_list(ndom_list)
    for dom_idx in range(ndomains):
        i = dom_idx % ndom_list[0]
        j = int(dom_idx / ndom_list[0])
        k = int(dom_idx / (ndom_list[0] * ndom_list[1]))

        # load full mesh info
        subdom_full_dir = os.path.join(decompdir, f"domain_{dom_idx}")
        coords_full_sub[i][j][k] = load_mesh_single(subdom_full_dir)
        nx, ny = coords_full_sub[i][j][k].shape[:-1]
        connect_full_sub[i][j][k] = np.loadtxt(os.path.join(subdom_full_dir, "connectivity.dat"), dtype=np.int32)[:, 1:]
        ndim = coords_full_sub[i][j][k].shape[-1]
        dcell_sub[i][j][k] = [None for _ in range(ndim)]
        with open(os.path.join(subdom_full_dir, "info.dat"), "r") as f:
            for line in f.readlines():
                if "dx" in line:
                    val = float(line.strip().split(" ")[1])
                    dcell_sub[i][j][k][0] = val
                if "dy" in line:
                    val = float(line.strip().split(" ")[1])
                    dcell_sub[i][j][k][1] = val
                if "dz" in line:
                    val = float(line.strip().split(" ")[1])
                    dcell_sub[i][j][k][2] = val

        # load sampled stencil mesh info
        sampdir_sub = os.path.join(sampdir, f"domain_{dom_idx}")
        connect_stencil_sub[i][j][k] = np.loadtxt(os.path.join(sampdir_sub, "connectivity.dat"), dtype=np.int32)[:, 1:]
        gids_samp_sub[i][j][k] = np.loadtxt(os.path.join(sampdir_sub, "sample_mesh_gids.dat"), dtype=np.int32)
        gids_stencil_sub[i][j][k] = np.loadtxt(os.path.join(sampdir_sub, "stencil_mesh_gids.dat"), dtype=np.int32)

        xmin = np.min(coords_full_sub[i][j][k][:, :, 0]) - dcell_sub[i][j][k][0] / 2
        xmax = np.max(coords_full_sub[i][j][k][:, :, 0]) + dcell_sub[i][j][k][0] / 2
        ymin = np.min(coords_full_sub[i][j][k][:, :, 1]) - dcell_sub[i][j][k][1] / 2
        ymax = np.max(coords_full_sub[i][j][k][:, :, 1]) + dcell_sub[i][j][k][1] / 2

        # get full domain boundaries
        xmin_full = min(xmin_full, xmin)
        xmax_full = max(xmax_full, xmax)
        ymin_full = min(ymin_full, ymin)
        ymax_full = max(ymax_full, ymax)

        for gid_idx, gid_stencil in enumerate(gids_stencil_sub[i][j][k]):
            if gid_stencil in gids_samp_sub[i][j][k]:
                is_sample = True
            else:
                is_sample = False

            xidx = gid_stencil % nx
            yidx = int(gid_stencil / nx)

            coords_samp = coords_full_sub[i][j][k][xidx, yidx, :]

            xmin_cell = coords_samp[0] - dcell_sub[i][j][k][0] / 2
            ymin_cell = coords_samp[1] - dcell_sub[i][j][k][1] / 2
            edgetuple = mcolors.to_rgb("k") + (1.0,)
            facetuple = mcolors.to_rgb(colors[dom_idx]) + (0.5,)
            rect = Rectangle(
                (xmin_cell, ymin_cell),
                dcell_sub[i][j][k][0], dcell_sub[i][j][k][1],
                edgecolor=edgetuple,
                facecolor=facetuple, fill=is_sample,
                linewidth=1.5,
            )
            ax.add_patch(rect)

    for dom_idx in range(ndomains):
        i = dom_idx % ndom_list[0]
        j = int(dom_idx / ndom_list[0])
        k = int(dom_idx / (ndom_list[0] * ndom_list[1]))
        xmin = np.min(coords_full_sub[i][j][k][:, :, 0]) - dcell_sub[i][j][k][0] / 2
        xmax = np.max(coords_full_sub[i][j][k][:, :, 0]) + dcell_sub[i][j][k][0] / 2
        ymin = np.min(coords_full_sub[i][j][k][:, :, 1]) - dcell_sub[i][j][k][1] / 2
        ymax = np.max(coords_full_sub[i][j][k][:, :, 1]) + dcell_sub[i][j][k][1] / 2
        rect = Rectangle(
            (xmin, ymin),
            xmax-xmin, ymax-ymin,
            linestyle="--", edgecolor=colors[dom_idx],
            fill=False, linewidth=2)
        ax.add_patch(rect)

    ax.set_xlim([xmin_full, xmax_full])
    ax.set_ylim([ymin_full, ymax_full])
    if outdir is None:
        plt.show()
    else:
        outfile = os.path.join(outdir, f"samp_decomp_mesh{outstr}.png")
        print(f"Saving image to {outfile}")
        plt.savefig(outfile)

if __name__ == "__main__":

    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)

    # location of full decomposed mesh
    parser.add_argument(
        "--decompDir", "--decompdir", "--decomp_dir",
        dest="decompdir",
        help="Full path to base directory of FULL decomposed mesh.",
    )
    parser.add_argument(
        "--sampDir", "--sampdir", "--samp_dir",
        dest="sampdir",
        help="Full path to base directory of SAMPLED decomposed mesh.",
    )
    parser.add_argument(
        "--outDir", "--outdir", "--out_dir",
        dest="outdir",
        help="Full path to directory where image will be saved.",
    )

    argobj = parser.parse_args()

    main(argobj.decompdir, argobj.sampdir, argobj.outdir)
