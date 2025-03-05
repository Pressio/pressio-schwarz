from copy import copy, deepcopy
import numpy as np

from pschwarz.data_utils import get_nested_decomp_dims, load_info_domain, make_empty_domain_list
from pschwarz.data_utils import load_unified_helper, load_meshes, decompose_domain_data

def calc_shared_samp_interval(
    dtlist=None,
    samplist=None,
    startlist=None,
    ndata=1,
):
    if (dtlist is not None) and (samplist is not None):
        if isinstance(dtlist, list):
            assert len(dtlist) == ndata
            dtlist = np.array(dtlist, dtype=np.float64)
        else:
            dtlist = dtlist * np.ones(ndata, dtype=np.float64)
        if isinstance(samplist, list):
            assert len(samplist) == ndata
            samplist = np.array(samplist, dtype=np.float64)
        else:
            samplist = samplist * np.ones(ndata, dtype=np.float64)

        # make sure there's a universally shared interval
        samplengths = dtlist * samplist
        sampintervals = np.amax(samplengths) / samplengths
        assert all([samp.is_integer() for samp in sampintervals])
        sampintervals = sampintervals.astype(np.int32)

    else:
        sampintervals = np.ones(ndata, dtype=np.int32)

    # make sure the starting index falls on an interval
    if startlist is not None:
        assert len(startlist) == ndata
        assert all([(start % interval) == 0 for start, interval in zip(startlist, sampintervals)])

    return samplengths, sampintervals

def calc_error_fields(
    meshlist=None,
    datalist=None,
    meshdirs=None,
    datadirs=None,
    nvars=None,
    dataroot=None,
    dtlist=None,
    samplist=None,
    startlist_in=None,
    mismatch_nan=False,
    merge_decomp=False,
):
    # NOTE: the FIRST data datalist/datadirs element is treated as "truth"

    _, datalist_in = load_unified_helper(
        meshlist=meshlist,
        datalist=datalist,
        meshdirs=meshdirs,
        datadirs=datadirs,
        nvars=nvars,
        dataroot=dataroot,
        merge_decomp=merge_decomp,
    )
    ndata = len(datalist_in)
    assert ndata > 1
    ndim = datalist_in[0].ndim - 2

    # check if any datasets are decomposed, check for meshdirs required for decomposition
    for data_idx, data in enumerate(datalist_in):
        if isinstance(data, list):
            assert meshdirs is not None
            assert isinstance(meshdirs, list)
            assert len(meshdirs) == ndata
            assert isinstance(meshdirs[data_idx], str)

    # If samplist and dtlist provided, comparison interval is explicit
    # Otherwise, same dt and sampling interval assumed
    samplengths, sampintervals = calc_shared_samp_interval(
        dtlist=dtlist,
        samplist=samplist,
        startlist=startlist_in,
        ndata=ndata,
    )

    if startlist_in is None:
        startlist = [0 for _ in range(ndata)]
    else:
        startlist = startlist_in


    # compute SIGNED errors (comparison - truth, not absolute)
    nsamps = None
    errorlist = []
    for data_idx in range(1, ndata):

        if isinstance(datalist_in[data_idx], list):
            ndom_list, overlap = load_info_domain(meshdirs[data_idx])
            _, meshlist_decomp = load_meshes(meshdirs[data_idx], merge_decomp=False)

            datalist_fom = decompose_domain_data(datalist_in[0], meshlist_decomp, overlap, is_ts=True, is_ts_decomp=False)
            datalist_comp = datalist_in[data_idx]

        else:
            ndom_list = [1, 1, 1]
            datalist_fom = [[[datalist_in[0]]]]
            datalist_comp = [[[datalist_in[data_idx]]]]

        ndomains = np.prod(ndom_list)

        errorlist_doms = make_empty_domain_list(ndom_list)
        for dom_idx in range(ndomains):
            i = dom_idx % ndom_list[0]
            j = int(dom_idx / ndom_list[0])
            k = int(dom_idx / (ndom_list[0] * ndom_list[1]))

            # short circuit any simulation that exploded
            if np.isnan(datalist_comp[i][j][k]).any():
                if ndim == 2:
                    fomshape = datalist_fom[i][j][k][:, :, startlist[0]::sampintervals[0], :].shape
                    error = np.empty(fomshape)
                    error[:] = np.nan
                else:
                    raise ValueError(f"Invalid ndim: {ndim}")
            else:
                if ndim == 2:
                    try:
                        error = datalist_comp[i][j][k][:, :, startlist[data_idx]::sampintervals[data_idx], :] - \
                            datalist_fom[i][j][k][:, :, startlist[0]::sampintervals[0], :]
                    except ValueError as e:
                        if mismatch_nan:
                            fomshape = datalist_fom[i][j][k][:, :, startlist[0]::sampintervals[0], :].shape
                            error = np.empty(fomshape)
                            error[:] = np.nan
                        else:
                            print("===============")
                            print("There was a dimension mismatch, check startlist and samplist")
                            print("===============")
                            raise(e)
                else:
                    raise ValueError(f"Invalid ndim: {ndim}")

            # double check that everything matches up
            if nsamps is None:
                nsamps = error.shape[-2]
            else:
                assert error.shape[-2] == nsamps

            errorlist_doms[i][j][k] = error.copy()

        if ndomains == 1:
            errorlist.append(errorlist_doms[0][0][0].copy())
        else:
            errorlist.append(deepcopy(errorlist_doms))

    # get sample times (ignoring any possible offset), for later plotting
    if (dtlist is not None) and (samplist is not None):
        samptimes = np.arange(nsamps) * np.amax(samplengths)
    else:
        samptimes = np.arange(nsamps)

    return errorlist, samptimes

def calc_error_norms(
    errorlist=None,
    samptimes=None,
    meshlist=None,
    datalist=None,
    meshdirs=None,
    datadirs=None,
    nvars=None,
    dataroot=None,
    dtlist=None,
    samplist=None,
    startlist=None,
    timenorm=False,
    spacenorm=False,
    relative=False,
    mismatch_nan=False,
    merge_decomp=False,
):
    assert timenorm or spacenorm

    # if computing relative norm, need data for denominator
    if relative:
        # this assert is a little weird, but don't want to potentially mix up error fields from different sources
        # can see an instance of passing an error list that doesn't correspond to the same datalist
        assert errorlist is None

        _, datalist = load_unified_helper(
            meshlist=meshlist,
            datalist=datalist,
            meshdirs=meshdirs,
            datadirs=datadirs,
            nvars=nvars,
            dataroot=dataroot,
            merge_decomp=merge_decomp,
        )
        # only need the truth value
        datatruth = datalist[0]

    if errorlist is None:
        errorlist, samptimes = calc_error_fields(
            meshlist=meshlist,
            datalist=datalist,
            meshdirs=meshdirs,
            datadirs=datadirs,
            nvars=nvars,
            dataroot=dataroot,
            dtlist=dtlist,
            samplist=samplist,
            startlist_in=startlist,
            mismatch_nan=mismatch_nan,
            merge_decomp=merge_decomp,
        )
        nsamps = samptimes.shape[0]
    else:
        if not isinstance(errorlist, list):
            errorlist = [errorlist]
        if isinstance(errorlist[0], list):
            nsamps = errorlist[0][0][0][0].shape[-2]
        else:
            nsamps = errorlist[0].shape[-2]
        assert all([error.shape[-2] == nsamps for error in errorlist])
        if samptimes is None:
            samptimes = np.arange(nsamps)
        else:
            assert samptimes.shape[0] == nsamps

    if isinstance(errorlist[0], list):
        ndim = errorlist[0][0][0][0].ndim - 2
        nvars_in = errorlist[0][0][0][0].shape[-1]
    else:
        ndim = errorlist[0].ndim - 2
        nvars_in = errorlist[0].shape[-1]
    space_axes = tuple(range(ndim))
    time_axis = ndim

    # need same sampling rate for time norm
    if relative:
        _, sampintervals = calc_shared_samp_interval(
            dtlist=dtlist,
            samplist=samplist,
            startlist=startlist,
            ndata=len(datalist),
        )

    # relative error scaling factors
    if relative:
        if startlist is None:
            startlist = [0]
        datanorm = datatruth[:, :, startlist[0]::sampintervals[0], :]
        if timenorm:
            relfacs = np.linalg.norm(datanorm, ord=2, axis=time_axis, keepdims=True)
        else:
            relfacs = datanorm.copy()
        if spacenorm:
            relfacs = np.linalg.norm(relfacs, ord=2, axis=space_axes, keepdims=True)
    else:
        relfacs = 1.0

    # compute norms
    for error_idx, error in enumerate(errorlist):

        if isinstance(error, list):
            ndom_list = get_nested_decomp_dims(error)
            error_doms = error
        else:
            ndom_list = [1, 1, 1]
            error_doms = [[[error]]]

        ndomains = np.prod(ndom_list)

        errorlist_doms = make_empty_domain_list(ndom_list)
        for dom_idx in range(ndomains):
            i = dom_idx % ndom_list[0]
            j = int(dom_idx / ndom_list[0])
            k = int(dom_idx / (ndom_list[0] * ndom_list[1]))

            error_dom = error_doms[i][j][k]

            if timenorm:
                error_dom = np.sqrt(np.sum(np.square(error_dom), axis=time_axis, keepdims=True))
            if spacenorm:
                error_dom = np.sqrt(np.sum(np.square(error_dom), axis=space_axes, keepdims=True))

            if np.isnan(error_dom).any():
                error_dom[:] = 1.0
            else:
                error_dom = error_dom / relfacs

            error_dom = np.squeeze(error_dom)
            if nvars_in == 1:
                error_dom = np.expand_dims(error_dom, -1)

            errorlist_doms[i][j][k] = copy(error_dom)

        if ndomains == 1:
            errorlist[error_idx] = copy(errorlist_doms[0][0][0])
        else:
            errorlist[error_idx] = deepcopy(errorlist_doms)

    return errorlist, samptimes