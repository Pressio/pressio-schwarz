import os
from math import floor, ceil

import numpy as np
from scipy.linalg import qr

from pschwarz.data_utils import load_meshes, load_info_domain
from pschwarz.prom_utils import load_pod_basis


def get_bound_indices(
    meshdims,
    dom_idx=None,
    ndom_list=None,
    phys_bounds=True,
    phys_rate=1,
    dom_bounds=True,
    dom_rate=1,
):

    # for decomposed mesh
    if dom_idx is not None:
        assert ndom_list is not None
        i = dom_idx % ndom_list[0]
        j = int(dom_idx / ndom_list[0])
        k = int(dom_idx / (ndom_list[0] * ndom_list[1]))
        assert phys_bounds or dom_bounds

    ndim = len(meshdims)

    phys_samples = []
    dom_samples = []
    if ndim == 1:
        samples_left = [0]
        samples_right = [meshdims[0]]

        if dom_idx is not None:
            if phys_bounds:
                if i == 0:
                    phys_samples += samples_left
                if i == (ndom_list[0] - 1):
                    phys_samples += samples_right
            if dom_bounds:
                if i != 0:
                    dom_samples += samples_left
                if i != (ndom_list[0] - 1):
                    dom_samples += samples_right
        else:
            if phys_bounds:
                phys_samples += samples_left + samples_right

    elif ndim == 2:
        samples_left  = [i * meshdims[0] for i in range(meshdims[1])]
        samples_top   = [i for i in range((meshdims[1] - 1) * meshdims[0], meshdims[0] * meshdims[1])]
        samples_right = [(meshdims[0] - 1) + i * meshdims[0] for i in range(meshdims[1])]
        samples_bot   = [i for i in range(meshdims[0])]

        if dom_idx is not None:
            if phys_bounds:
                if i == 0:
                    phys_samples += samples_left[::phys_rate]
                if i == (ndom_list[0] - 1):
                    phys_samples += samples_right[::phys_rate]
                if j == 0:
                    phys_samples += samples_bot[::phys_rate]
                if j == (ndom_list[1] - 1):
                    phys_samples += samples_top[::phys_rate]

            if dom_bounds:
                if i != 0:
                    dom_samples += samples_left[::dom_rate]
                if i != (ndom_list[0] - 1):
                    dom_samples += samples_right[::dom_rate]
                if j!= 0:
                    dom_samples += samples_bot[::dom_rate]
                if j != (ndom_list[1] - 1):
                    dom_samples += samples_top[::dom_rate]
        else:
            if phys_bounds:
                phys_samples += samples_left[::phys_rate] + samples_top[::phys_rate] + samples_right[::phys_rate] + samples_bot[::phys_rate]

    elif ndim == 3:
        raise ValueError("3D seed samples not implemented yet")

    return list(np.unique(phys_samples)), list(np.unique(dom_samples))


def sample_domain(
    samptype,
    coords_full,
    coords_sub,
    percpoints,
    outdir,
    basis=None,
    seed_qdeim=False,
    seed_phys_bounds=False,
    seed_phys_rate=1,
    seed_dom_bounds=False,
    seed_dom_rate=1,
    samp_phys_bounds=True,
    samp_dom_bounds=True,
    randseed=0,
    dom_idx=None,
    ndom_list=None
):

    assert not (samp_phys_bounds and seed_phys_bounds)
    assert not (samp_dom_bounds and seed_dom_bounds)
    if seed_phys_bounds:
        assert seed_phys_rate > 0
    if seed_dom_bounds:
        assert seed_dom_rate > 0

    # monolithic
    if dom_idx is None:
        decomp = False
        assert (not seed_dom_bounds) and (not samp_dom_bounds), \
            "Passed samp_dom_bounds=True or seed_dom_bounds=True for a monolithic domain, check your inputs"
        meshdims = coords_full.shape[:-1]

    # decomposed
    else:
        decomp = True
        assert ndom_list is not None
        ndomains = np.prod(ndom_list)

        i = dom_idx % ndom_list[0]
        j = int(dom_idx / ndom_list[0])
        k = int(dom_idx / (ndom_list[0] * ndom_list[1]))
        coords_local = coords_sub[i][j][k]
        meshdims = coords_local.shape[:-1]

    ncells = np.prod(meshdims)

    # do bounds without subsampling for now
    all_ids = np.arange(ncells)
    phys_bound_ids, dom_bound_ids = get_bound_indices(
        meshdims,
        dom_idx=dom_idx,
        ndom_list=ndom_list
    )
    interior_ids = np.setdiff1d(np.setdiff1d(all_ids, phys_bound_ids), dom_bound_ids)
    ncells_phys = len(phys_bound_ids)
    ncells_dom = len(dom_bound_ids)
    ncells_int  = len(interior_ids)

    if seed_phys_bounds or seed_dom_bounds:
        phys_bound_ids_samp, dom_bound_ids_samp = get_bound_indices(
            meshdims,
            dom_idx=dom_idx,
            ndom_list=ndom_list,
            phys_bounds=seed_phys_bounds,
            phys_rate=seed_phys_rate,
            dom_bounds=seed_dom_bounds,
            dom_rate=seed_dom_rate,
        )

    # "seed" sample indices
    points_seed = []
    if seed_qdeim:
        # don't bother with separating boundaries for QDEIM
        points_seed += calc_qdeim_samples(basis, ncells)
    if seed_phys_bounds:
        points_seed += phys_bound_ids_samp
    if seed_dom_bounds:
        points_seed += dom_bound_ids_samp
    points_seed = list(np.unique(points_seed))
    npoints_seed = len(points_seed)

    # divvy up remaining cells
    npoints_tot = floor(ncells * percpoints)
    assert npoints_tot >= npoints_seed, f"{npoints_tot} < {npoints_seed}"
    npoints_remain = npoints_tot - npoints_seed

    samp_ratio_phys = ncells_phys / ncells
    samp_ratio_dom = ncells_dom / ncells
    samp_ratio_int = 1.0

    # have to do these first to avoid messing math up later
    if samp_phys_bounds:
        samp_ratio_int -= samp_ratio_phys
    if samp_dom_bounds:
        samp_ratio_int -= samp_ratio_dom

    if not samp_phys_bounds:
        if samp_dom_bounds:
            samp_ratio_dom += samp_ratio_phys / 2
            samp_ratio_int -= samp_ratio_phys / 2
        samp_ratio_phys = 0.0
    if not samp_dom_bounds:
        if samp_phys_bounds:
            samp_ratio_phys += samp_ratio_dom / 2
            samp_ratio_int -= samp_ratio_dom / 2
        samp_ratio_dom = 0.0

    # override if boundaries already sampled
    if seed_phys_bounds:
        samp_ratio_phys = 0.0
    if seed_dom_bounds:
        samp_ratio_dom = 0.0
    npoints_phys = round(npoints_remain * samp_ratio_phys)
    npoints_dom  = round(npoints_remain * samp_ratio_dom)
    npoints_int  = round(npoints_remain * samp_ratio_int)

    if samptype == "random":
        # don't bother separating boundaries for random sampling
        samples = gen_random_samples(0, ncells-1, npoints_tot, randseed=randseed, points_seed=points_seed)
    elif samptype == "eigenvec":
        samples = calc_eigenvec_samples(
            basis, ncells,
            npoints_seed+npoints_int,
            points_seed=points_seed, search_cell_ids=interior_ids
        )
        if npoints_phys > 0:
            samples = calc_eigenvec_samples(
                basis, ncells,
                npoints_seed+npoints_int+npoints_phys,
                points_seed=samples, search_cell_ids=phys_bound_ids
            )
        if npoints_dom > 0:
            samples = calc_eigenvec_samples(
                basis, ncells,
                npoints_seed+npoints_int+npoints_phys+npoints_dom,
                points_seed=samples, search_cell_ids=dom_bound_ids
            )
        assert samples.shape[0] == (npoints_seed + npoints_int + npoints_phys + npoints_dom)

    samples = np.unique(samples)
    outfile = os.path.join(outdir, "sample_mesh_gids.dat")
    print(f"Saving sample mesh global indices to {outfile}")
    np.savetxt(outfile, samples, fmt='%8i')


def gen_sample_mesh(
    samptype,
    meshdir,
    percpoints,
    outdir,
    basis_dir=None,
    nmodes=None,
    seed_qdeim=False,
    seed_phys_bounds=False,
    seed_phys_rate=1,
    seed_dom_bounds=False,
    seed_dom_rate=1,
    samp_phys_bounds=True,
    samp_dom_bounds=True,
    randseed=0,
    recalc=True,
):

    # expand as necessary
    assert samptype in ["random", "eigenvec", "gnat"]
    assert (percpoints > 0.0) and (percpoints <= 1.0)

    load_basis = False
    if seed_qdeim or (samptype in ["eigenvec", "gnat"]):
        load_basis = True
        assert basis_dir is not None
        assert os.path.isdir(basis_dir), f"No basis directory at {basis_dir}"
        assert nmodes is not None

    if os.path.isdir(outdir) and (not recalc):
        print(f"Not recalculating sample mesh at {outdir}")
        return
    else:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

    # get full mesh dimensions
    coords_full, coords_sub = load_meshes(meshdir)

    # monolithic sample mesh
    if coords_sub is None:

        if load_basis:
            assert isinstance(nmodes, int)
            basis, = load_pod_basis(basis_dir, nmodes=nmodes)
        else:
            basis = None

        sample_domain(
            samptype,
            coords_full,
            coords_sub,
            percpoints,
            outdir,
            basis=basis,
            seed_qdeim=seed_qdeim,
            seed_phys_bounds=seed_phys_bounds,
            seed_phys_rate=seed_phys_rate,
            seed_dom_bounds=False,
            samp_phys_bounds=samp_phys_bounds,
            samp_dom_bounds=False,
            randseed=randseed,
        )

    # decomposed sample mesh
    else:

        ndom_list, _ = load_info_domain(meshdir)
        ndomains = np.prod(ndom_list)

        if load_basis:
            if isinstance(nmodes, int):
                nmodes = [nmodes] * ndomains
            else:
                assert isinstance(nmodes, list)
            assert len(nmodes) == ndomains
            basis_list, = load_pod_basis(basis_dir, nmodes=nmodes)
        else:
            basis_list = [None for _ in range(ndomains)]

        for dom_idx in range(ndomains):

            outdir_samps = os.path.join(outdir, f"domain_{dom_idx}")
            if not os.path.isdir(outdir_samps):
                os.mkdir(outdir_samps)

            sample_domain(
                samptype,
                coords_full,
                coords_sub,
                percpoints,
                outdir_samps,
                basis=basis_list[dom_idx],
                seed_qdeim=seed_qdeim,
                seed_phys_bounds=seed_phys_bounds,
                seed_phys_rate=seed_phys_rate,
                seed_dom_bounds=seed_dom_bounds,
                seed_dom_rate=seed_dom_rate,
                samp_phys_bounds=samp_phys_bounds,
                samp_dom_bounds=samp_dom_bounds,
                randseed=randseed+dom_idx,
                dom_idx=dom_idx,
                ndom_list=ndom_list,
            )


def gen_random_samples(
    low,
    high,
    numsamps,
    randseed=0,
    points_seed=[],
):

    rng = np.random.default_rng(seed=randseed)

    all_samples = np.arange(low, high+1)
    rng.shuffle(all_samples)

    numsamps_seed = len(points_seed)
    if numsamps_seed == 0:
        samples = np.sort(all_samples[:numsamps])
    else:
        samples = points_seed
        idx = 0
        while len(samples) < numsamps:
            if all_samples[idx] not in samples:
                samples.append(all_samples[idx])
            idx += 1
        samples = np.sort(samples)

    return samples


def calc_qdeim_samples(basis, ncells, search_cell_ids=None):
    nmodes = basis.shape[1]

    _, _, dof_ids = qr(basis.T, pivoting=True, mode="economic")
    samp_ids = dof_ids % ncells

    if search_cell_ids is None:
        samp_ids_out = samp_ids[:nmodes]
    else:
        assert len(search_cell_ids) >= nmodes
        samp_ids_out = []
        samp_idx = 0
        while len(samp_ids_out) < nmodes:
            if samp_ids[samp_idx] in search_cell_ids:
                samp_ids_out.append(samp_ids[samp_idx])
            samp_idx += 1

    return list(np.unique(samp_ids_out))


def calc_eigenvec_samples(
    basis,
    ncells,
    numsamps,
    points_seed=[],
    search_cell_ids=None,
    randseed=0,
    printrate=0.01,
):

    print("Starting eigenvector sampling")

    ndof_per_cell = round(basis.shape[0] / ncells)
    nmodes = basis.shape[-1]

    # if no seed provided, just sample nmodes/ndof_per_cell random cells to start things
    # need that many cells to make sure right svecs have correct shape
    if len(points_seed) == 0:
        rng = np.random.default_rng(seed=randseed)
        all_samples = np.arange(ncells)
        rng.shuffle(all_samples)
        points_seed += list(all_samples[:ceil(nmodes / ndof_per_cell)])

    samp_ids = np.array(points_seed, dtype=np.int32)

    if search_cell_ids is None:
        search_cell_ids = [i for i in range(ncells)]
    search_dof_ids  = np.concatenate([search_cell_ids * (i + 1) for i in range(ndof_per_cell)])
    ncells_search = len(search_cell_ids)

    rateperc = max(int(numsamps * printrate), 1)
    while samp_ids.shape[0] < numsamps:

        if (samp_ids.shape[0] % rateperc) == 0:
            print(f"Sampled {(samp_ids.shape[0] / numsamps) * 100:#.2f}% samples")

        dof_ids = np.concatenate([samp_ids * (i + 1) for i in range(ndof_per_cell)])
        basis_samp = basis[dof_ids, :]
        _, _, svecs_right = np.linalg.svd(basis_samp, full_matrices=False)
        resvec = np.squeeze(np.square(svecs_right[:, [-1]].T @ basis[search_dof_ids, :].T))

        # compute cell average
        resvec = np.reshape(resvec, (ncells_search, -1), order="F")
        resvec = np.mean(resvec, axis=-1)
        sort_idxs = np.flip(np.argsort(resvec))

        for sort_idx in sort_idxs:
            cell_id_samp = search_cell_ids[sort_idx]
            if cell_id_samp not in samp_ids:
                samp_ids = np.append(samp_ids, [cell_id_samp])
                samp_ids = np.sort(samp_ids)
                break

    return samp_ids

