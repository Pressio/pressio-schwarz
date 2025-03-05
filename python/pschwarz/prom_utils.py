import os

import numpy as np
from scipy.linalg import svd

from pschwarz.data_utils import load_meshes, load_info_domain, decompose_domain_data, merge_domain_data
from pschwarz.data_utils import load_unified_helper, write_to_binary, read_from_binary, make_empty_domain_list


def center(data_in, centervec=None, method=None):
    # data assumed to have shape (spatial, time, var)

    if centervec is None:
        assert method is not None
        if method == "zero":
            centervec = np.zeros((data_in.shape[0], 1, data_in.shape[-1]), dtype=np.float64)
        elif method == "init_cond":
            centervec = data_in[:, [0], :]
        elif (method == "mean"):
            centervec = np.mean(data_in, axis=1, keepdims=True)
        else:
            raise ValueError(f"Invalid centering method: {method}")
    else:
        # just get in shape to broadcast
        assert centervec.ndim == 2
        centervec = centervec[:, None, :]

    data_out = data_in - centervec

    return data_out, np.squeeze(centervec)


def normalize(data_in, normvec=None, method=None):
    # data assumed to have shape (spatial, time, var)

    if normvec is None:
        if method == "one":
            normvec = np.ones((data_in.shape[0], 1, data_in.shape[-1]), dtype=np.float64)
        elif method == "l2":
            normvec = np.mean(
                np.square(np.linalg.norm(data_in, axis=0, ord=2, keepdims=True)),
                axis=1,
                keepdims=True,
            ) / data_in.shape[0]
            normvec = np.repeat(normvec, data_in.shape[0], axis=0)
        else:
            raise ValueError(f"Invalid normalization method: {method}")
    else:
        # just get in shape to broadcast
        assert normvec.ndim == 2
        normvec = normvec[:, None, :]

    data_out = data_in / normvec

    return data_out, np.squeeze(normvec)


def calc_pod_single(
    data_in,
    centervec_in=None,
    center_method=None,
    normvec_in=None,
    norm_method=None,
    nmodes=None,
):
    assert (centervec_in is not None) or (center_method is not None)
    assert (normvec_in is not None) or (norm_method is not None)

    # flatten spatial dimension (I/O is column-major)
    dim = data_in.ndim - 2
    if dim == 2:
        data = np.reshape(data_in, (-1,) + data_in.shape[-2:], order="F")
    else:
        raise ValueError(f"Unsupported dimension: {dim}")

    # center, normalize data
    data_proc, centervec = center(data, centervec=centervec_in, method=center_method)
    data_proc, normvec = normalize(data_proc, normvec=normvec_in, method=norm_method)

    # TODO: could do scalar POD here
    # flatten variable dimensions
    nsamps = data_proc.shape[1]
    data_proc = np.transpose(data_proc, (0, 2, 1)) # spatial, variable, time
    data_proc = np.reshape(data_proc, (-1, nsamps), order="C")

    # compute POD basis, truncate
    U, S, _ = svd(data_proc, full_matrices=False)
    U = U[:, :nmodes]

    # bake normalization into basis
    U = normvec.flatten(order="F")[:, None] * U

    # return centering/normalization vectors and basis
    return U, S, centervec, normvec


def gen_pod_bases(
    outdir,
    meshlist=None,
    datalist=None,
    meshdir=None,
    datadir=None,
    nvars=None,
    dataroot=None,
    concat=False,
    pod_decomp=False,
    meshdir_decomp=None,
    idx_start=0,
    idx_end=None,
    idx_skip=1,
    centervec_in=None,
    center_method=None,
    normvec_in=None,
    norm_method=None,
    nmodes=None,
):

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    meshlist, datalist = load_unified_helper(
        meshlist,
        datalist,
        meshdir,
        datadir,
        nvars,
        dataroot,
        merge_decomp=False,
    )
    ndim = datalist[0].ndim - 2

    # downsample in time
    for idx, _ in enumerate(datalist):
        if ndim == 2:
            datalist[idx] = datalist[idx][:, :, idx_start:idx_end:idx_skip, :]
        else:
            raise ValueError(f"Unsupported ndim = {ndim}")

    # concatenate datasets along time dimension
    if concat:
        datalist = [np.concatenate(datalist, axis=ndim)]

    # if doing decomposed POD, either must be mono solution or share mesh
    ndata_in = len(datalist)
    assert ndata_in >= 1
    if pod_decomp:

        if ndata_in == 1:
            assert meshdir_decomp is not None
            _, overlap = load_info_domain(meshdir_decomp)
            _, meshlist_decomp = load_meshes(meshdir_decomp, merge_decomp=False)

            datalist_decomp = decompose_domain_data(datalist[0], meshlist_decomp, overlap, is_ts=True, is_ts_decomp=False)

            # unroll from 3D list
            datalist = []
            for k in range(len(datalist_decomp[0][0])):
                for j in range(len(datalist_decomp[0])):
                    for i in range(len(datalist_decomp)):
                        datalist.append(datalist_decomp[i][j][k])

        else:
            # don't need do do anything, already have decomposed solution
            assert (meshlist_decomp is None) or (meshdir_decomp is None)

    ndata_out = len(datalist)
    if centervec_in is None:
        centervec_in = [None for _ in range(ndata_out)]
    if normvec_in is None:
        normvec_in = [None for _ in range(ndata_out)]

    print(f"Writing basis to {outdir}")
    for data_idx, data in enumerate(datalist):

        # compute basis and feature scaling vectors
        basis, svals, centervec, normvec = calc_pod_single(
            data,
            centervec_in=centervec_in[data_idx],
            center_method=center_method,
            normvec_in=normvec_in[data_idx],
            norm_method=norm_method,
            nmodes=nmodes,
        )

        # write to disk
        if (ndata_out == 1) and (not pod_decomp):
            numstr = ""
        else:
            numstr = f"_{data_idx}"
        basis_file = os.path.join(outdir, f"basis{numstr}.bin")
        # FIXME: this transpose and "reverse" is bad practice
        write_to_binary(basis.T, basis_file, reverse=True)
        svec_file = os.path.join(outdir, f"svals{numstr}.npy")
        np.save(svec_file, svals)
        center_file = os.path.join(outdir, f"center{numstr}.bin")
        write_to_binary(centervec.flatten(order="C"), center_file)
        norm_file = os.path.join(outdir, f"norm{numstr}.bin")
        write_to_binary(normvec.flatten(order="C"), norm_file)


        # determine 99%, 99.9%, 99.99%
        sumsq = np.sum(np.square(svals))
        res = 1.0 - np.cumsum(np.square(svals)) / sumsq
        energy_file = os.path.join(outdir, f"pod_power{numstr}.dat")
        with open(energy_file, "w") as f:
            f.write(f"99.00%: {np.argwhere(res < 0.01)[0][0]}\n")
            f.write(f"99.90%: {np.argwhere(res < 0.001)[0][0]}\n")
            f.write(f"99.99%: {np.argwhere(res < 0.0001)[0][0]}\n")


def load_pod_basis(
    basisdir,
    return_basis=True,
    return_center=False,
    return_norm=False,
    return_svals=False,
    nmodes=None,
):

    # REMINDER: merging decomposed bases is NOT VALID

    assert os.path.isdir(basisdir), f"No basis directory at {basisdir}"
    print(f"Loading basis from {basisdir}")

    out_tuple = ()
    if os.path.isfile(os.path.join(basisdir, "basis.bin")):
        print("Monolithic basis detected")

        if return_basis:
            basis = read_from_binary(os.path.join(basisdir, "basis.bin"))[:, :nmodes]
            out_tuple += (basis,)
        if return_center:
            center = read_from_binary(os.path.join(basisdir, "center.bin"))
            out_tuple += (center,)
        if return_norm:
            norm = read_from_binary(os.path.join(basisdir, "norm.bin"))
            out_tuple += (norm,)
        if return_svals:
            svals = np.load(os.path.join(basisdir, "svals.npy"))
            out_tuple += (svals,)

    elif os.path.isfile(os.path.join(basisdir, "basis_0.bin")):
        print("Decomposed basis detected")

        if nmodes is not None:
            if isinstance(nmodes, int):
                nmodes = [nmodes] * 100
            assert isinstance(nmodes, list)

        dom_idx = 0
        basis = []
        center = []
        norm = []
        svals = []
        while os.path.isfile(os.path.join(basisdir, f"basis_{dom_idx}.bin")):
            if return_basis:
                basis_in = read_from_binary(os.path.join(basisdir, f"basis_{dom_idx}.bin"))[:, :nmodes[dom_idx]]
                basis.append(basis_in.copy())
            if return_center:
                center_in = read_from_binary(os.path.join(basisdir, f"center_{dom_idx}.bin"))
                center.append(center_in.copy())
            if return_norm:
                norm_in = read_from_binary(os.path.join(basisdir, f"norm_{dom_idx}.bin"))
                norm.append(norm_in.copy())
            if return_svals:
                svals_in = np.load(os.path.join(basisdir, f"svals_{dom_idx}.npy"))
                svals.append(svals_in.copy())
            dom_idx += 1
        print(f"{dom_idx} domain bases found")

        if return_basis:
            out_tuple += (basis,)
        if return_center:
            out_tuple += (center,)
        if return_norm:
            out_tuple += (norm,)
        if return_svals:
            out_tuple += (svals,)

    else:
        raise ValueError(f"No basis file(s) found in {basisdir}")

    return out_tuple


def load_reduced_data(
    datadir,
    fileroot,
    nvars,
    meshdir,
    trialdir,
    basisroot,
    centerroot,
    normroot,
    nmodes,
    basis_in=None,
    center_in=None,
    norm_in=None,
    merge_decomp=False,
):

    assert os.path.isdir(datadir), f"No data directory at {datadir}"
    print(f"Loading data from {datadir}")

    # detect monolithic vs decomposed
    if os.path.isfile(os.path.join(datadir, fileroot + ".bin")):

        print("Monolithic solution detected")
        coords, coords_sub = load_meshes(meshdir)
        assert coords_sub is None
        ndim = coords.shape[-1]
        meshdims = coords.shape[:-1]

        if basis_in is None:
            basis_file = os.path.join(trialdir, basisroot + ".bin")
            basis = read_from_binary(basis_file)[:, :nmodes]
        else:
            basis = basis_in[:, :nmodes]
        if center_in is None:
            center_file = os.path.join(trialdir, centerroot + ".bin")
            center = read_from_binary(center_file)
        else:
            center = center_in.copy()
        if norm_in is None:
            norm_file = os.path.join(trialdir, normroot + ".bin")
            norm = read_from_binary(norm_file)
        else:
            norm = norm_in.copy()

        data_red = np.fromfile(os.path.join(datadir, fileroot + ".bin"))
        data_red = np.reshape(data_red, (nmodes, -1), order="F")

        data_full = center + norm * (basis @ data_red)
        nsnaps = data_full.shape[-1]
        data_full = np.reshape(data_full, ((nvars,) + meshdims + (nsnaps,)), order="F")
        data_full = np.transpose(data_full, tuple(np.arange(1,ndim+1)) + (ndim+1, 0,))

    elif os.path.isfile(os.path.join(datadir, fileroot + "_0.bin")):

        print("Decomposed solution detected")

        ndom_list, overlap = load_info_domain(meshdir)
        coords, coords_sub = load_meshes(meshdir)
        ndim = coords_sub[0][0][0].shape[-1]

        assert isinstance(nmodes, list)
        ndom = len(nmodes)

        if any([val is None for val in [basis_in, center_in, norm_in]]):
            basis, center, norm = load_pod_basis(
                trialdir,
                return_basis=True,
                return_center=True,
                return_norm=True,
                nmodes=nmodes,
            )
        else:
            basis  = [val.copy() for val in basis_in]
            center = [val.copy() for val in center_in]
            norm   = [val.copy() for val in norm_in]

        assert len(basis) == ndom
        assert len(center) == ndom
        assert len(norm) == ndom

        # truncate
        for dom_idx in range(ndom):
            basis[dom_idx] = basis[dom_idx][:, :nmodes[dom_idx]]

        dom_idx = 0
        data_full = make_empty_domain_list(ndom_list)
        while os.path.isfile(os.path.join(datadir, fileroot + f"_{dom_idx}.bin")):

            i = dom_idx % ndom_list[0]
            j = int(dom_idx / ndom_list[0])
            k = int(dom_idx / (ndom_list[0] * ndom_list[1]))
            meshdims = coords_sub[i][j][k].shape[:-1]

            data_red = np.fromfile(os.path.join(datadir, fileroot + f"_{dom_idx}.bin"))
            data_red = np.reshape(data_red, (nmodes[dom_idx], -1), order="F")

            data_full_in = center[dom_idx] + norm[dom_idx] * (basis[dom_idx] @ data_red)
            nsnaps = data_full_in.shape[-1]
            data_full_in = np.reshape(data_full_in, ((nvars,) + meshdims + (nsnaps,)), order="F")
            data_full_in = np.transpose(data_full_in, tuple(np.arange(1,ndim+1)) + (ndim+1, 0,))

            data_full[i][j][k] = data_full_in.copy()

            dom_idx += 1

        assert dom_idx == ndom

        if merge_decomp:
            data_full = merge_domain_data(data_full, overlap, is_ts=True)

    else:
        raise ValueError(f"Could not find reduced data at {datadir}")

    return data_full


def project_single(
    data,
    basis,
    center,
    norm,
):

    assert basis.ndim == 2
    ndof = basis.shape[0]
    if center.ndim == 1:
        center = center[:, None]
    assert center.shape == (ndof, 1)
    if norm.ndim == 1:
        norm = norm[:, None]
    assert norm.shape == (ndof, 1)

    ndim = data.ndim - 2
    spatial_dims = data.shape[:ndim]
    nsnaps, nvars = data.shape[-2:]

    # flatten data
    data_flat = np.transpose(data, (ndim+1,) + tuple(np.arange(ndim+1)))
    data_flat = np.reshape(data_flat, (-1, nsnaps), order="F")

    # compute projection
    # TODO: modify for scalar POD
    data_red = basis.T @ ((data_flat - center) / norm)
    data_flat = center + norm * (basis @ data_red)

    # reshape data
    data_flat = np.reshape(data_flat, (nvars,) + spatial_dims + (-1,), order="F")
    data_out = np.transpose(data_flat, tuple(np.arange(1,ndim+2)) + (0,))

    return data_out, data_red


def calc_projection(
    nmodes,
    meshlist=None,
    datalist=None,
    meshdirs=None,
    datadirs=None,
    nvars=None,
    dataroot=None,
    basisdir=None,
    basis_in=None,
    center=None,
    norm=None,
    meshdir_decomp=None,
    merge_decomp=False,
):
    # NOTE: this is ONE basis applied to MULTIPLE datasets
    # Data to be projected will always be monolithic/merged
    # Basis may be decomposed, solution will be decomposed before projection

    assert nmodes >= 1

    # load data (if not provided)
    # always merge, will be decomposed later
    _, datalist = load_unified_helper(
        meshlist=meshlist,
        datalist=datalist,
        meshdirs=meshdirs,
        datadirs=datadirs,
        nvars=nvars,
        dataroot=dataroot,
        merge_decomp=True,
    )
    assert all([isinstance(data, np.ndarray) for data in datalist])

    # check dimensions, should all share same spatial dimensions
    for data_idx, data in enumerate(datalist):
        if data_idx == 0:
            ndim = data.ndim - 2
            spatial_dims = data.shape[:ndim]
            ndof = np.prod(spatial_dims) * nvars
        else:
            assert (data.ndim - 2) == ndim
            assert data.shape[:ndim] == spatial_dims

    # load basis, center, norm (if not provided)
    if any([inp is not None for inp in [basis_in, center, norm]]):
        assert all([inp is not None for inp in [basis_in, center, norm]])
        if isinstance(basis_in, list):
            basis = [basis_arr[:, :nmodes] for basis_arr in basis_in]
        elif isinstance(basis_in, np.ndarray):
            basis = basis_in[:, :nmodes]
        else:
            raise ValueError("Unexpected basis type")
    else:
        assert basisdir is not None
        basis, center, norm = load_pod_basis(
            basisdir,
            return_basis=True,
            return_center=True,
            return_norm=True,
            nmodes=nmodes,
        )

    datalist_out = []

    # with decomposed basis
    if isinstance(basis, list):
        assert meshdir_decomp is not None
        ndom_list, overlap = load_info_domain(meshdir_decomp)
        _, meshlist_decomp = load_meshes(meshdir_decomp, merge_decomp=False)
        ndomains = np.prod(ndom_list)
        assert len(basis) == ndomains
        assert len(center) == ndomains
        assert len(norm) == ndomains

        for data_idx, data in enumerate(datalist):

            # decompose data
            data_sub = decompose_domain_data(data, meshlist_decomp, overlap, is_ts=True, is_ts_decomp=False)

            # project each domain
            for dom_idx in range(ndomains):
                i = dom_idx % ndom_list[0]
                j = int(dom_idx / ndom_list[0])
                k = int(dom_idx / (ndom_list[0] * ndom_list[1]))

                data_sub[i][j][k], _ = project_single(
                    data_sub[i][j][k],
                    basis[dom_idx],
                    center[dom_idx],
                    norm[dom_idx],
                )

            if merge_decomp:
                # merge data
                datalist_out.append(merge_domain_data(data_sub, overlap, is_ts=True))
            else:
                datalist_out.append(data_sub)

    # with monolithic basis
    elif isinstance(basis, np.ndarray):
        assert basis.shape[0] == ndof
        assert center.shape[0] == ndof
        assert norm.shape[0] == ndof
        assert basis.shape[1] >= nmodes
        basis = basis[:, :nmodes]

        # project data
        for data_idx, data in enumerate(datalist):
            data_out, _ = project_single(
                datalist[data_idx],
                basis,
                center,
                norm,
            )
            datalist_out.append(data_out.copy())

    else:
        raise ValueError("Unexpected basis type")

    return datalist_out
