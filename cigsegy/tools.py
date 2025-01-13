# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import Dict, List, Tuple
import numpy as np
from cigsegy.cpp._CXX_SEGY import Pysegy
from cigsegy import utils
import warnings


def read_header(fname: str, type, n=0, printstr=True):
    """
    Read binary or trace header

    Parameters
    ----------
    segy : str
        input segy file
    type: str
        can be one of ['bh', 'th'],
            'bt' means binary header, 'th' means trace header,
    n : int
        trace number when type is 'th'
    printstr : bool
        print header information, if False, return a dict of header's infomation

    Returns
    -------
    Dict or None
    """
    if type == 'bh':
        try:
            segy = Pysegy(fname)
            arr = segy.get_binary_header()
        except:
            warnings.warn("The SEG-Y file is broken, try to read binary header by numpy")
            arr = np.fromfile(fname, dtype=np.uint8, count=400, offset=3200)

        out, hstring = utils.parse_bheader(arr)
    elif type == 'th':
        segy = Pysegy(fname)
        arr = segy.get_trace_header(n)
        out, hstring = utils.parse_theader(arr)
    else:
        raise ValueError("type must be one of ['bh', 'th']")

    if printstr:
        print('\n'.join(hstring))
        return
    return out


def get_metaInfo(
    segyname: str,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
    xloc: int = None,
    yloc: int = None,
    *,
    is4d: bool = None,
    apply_scalar: bool = False,
) -> Dict:
    """
    get metainfo dict of `segyname` file

    Parameters
    ----------
    segyname : str
        input segy file
    iline : int
        iline location in trace header
    iline : int
        iline location in trace header
    istep : int
        iline step, 2 means iline is like 100, 102, 104, ...
    xstep : int
        xline step
    xloc : int
        x (real world) value location in trace header
    yloc : int
        y (real world) value location in trace header
    apply_scalar : bool
        apply scalar to 'i-interval' and 'x-interval'

    Returns
    -------
    Dict
        Dict of meta information 
    """
    if isinstance(segyname, Pysegy):
        segy = segyname
    else:
        segy = Pysegy(str(segyname))

    [iline, xline, offset, istep, xstep, ostep, xloc, yloc, _is4d] = utils.guess(segy, iline, xline, offset, istep, xstep, ostep, xloc, yloc) # yapf: disable
    if is4d is None:
        is4d = _is4d
    # print(iline, xline, offset, istep, xstep, ostep, xloc, yloc, is4d)
    segy.setLocations(iline, xline, offset)
    segy.setSteps(istep, xstep, ostep)
    segy.setXYLocations(xloc, yloc)
    ndim = 4 if is4d else 3
    segy.set_segy_type(ndim)
    segy.scan()
    keys = segy.get_keylocs()
    meta = segy.get_metainfo()
    meta = {**keys, **meta}

    meta = utils.post_process_meta(segy, meta, apply_scalar)

    if not isinstance(segyname, Pysegy):
        segy.close()

    return meta


def get_lineInfo(
    fname: str,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
    xloc: int = None,
    yloc: int = None,
    mode: str = 'raw',
):
    """
    mode can be one of ['raw', 'geom']
    """
    if isinstance(fname, Pysegy):
        segy = fname
    else:
        segy = Pysegy(str(fname))

    [iline, xline, offset, istep, xstep, ostep, xloc, yloc, is4d] = utils.guess(segy, iline, xline, offset, istep, xstep, ostep, xloc, yloc) # yapf: disable
    segy.setLocations(iline, xline, offset)
    segy.setSteps(istep, xstep, ostep)
    segy.setXYLocations(xloc, yloc)
    ndim = 4 if is4d else 3
    segy.set_segy_type(ndim)
    segy.scan()
    lineinfo = segy.get_lineInfo()
    ndim = segy.ndim

    out = None
    if mode == 'raw':
        out = lineinfo
    elif mode == 'geom':
        if ndim == 4:
            raise NotImplementedError("4D geometry is not supported yet")
        lineinfo = lineinfo[~np.any(lineinfo == -1, axis=1)]
        N = lineinfo.shape[0]
        out = np.zeros((N * 2 + 1, 4), dtype=np.int32)
        out[:N, :2] = lineinfo[:, :2]
        for i in range(N):
            idx = lineinfo[i, 3]
            out[i, 2:] = [segy.coordx(idx), segy.coordy(idx)]

        lineinfo = lineinfo[::-1, ...]
        out[N:2 * N, :2] = lineinfo[:, [0, 2]]
        for i in range(N):
            idx = lineinfo[i, 4]
            out[N + i, 2:] = [segy.coordx(idx), segy.coordy(idx)]

        out[-1] = out[0]
    else:
        raise ValueError("`mode` only support 'raw' and 'geom'")

    if not isinstance(fname, Pysegy):
        segy.close()

    # HACK: add more mode
    return out


def trace_count(fname: str) -> int:
    """
    Count the total numbers of a segy file

    Parameters
    ----------
    segy: str
        input segy file

    Returns
    -------
    int
        The total numbers of a segy file
    """
    if isinstance(fname, Pysegy):
        return fname.ntrace
    segy = Pysegy(str(fname))
    count = segy.ntrace
    segy.close()
    return count


def full_scan(fname: str,
              iline: int = None,
              xline: int = None,
              offset: int = 37,
              keys: np.ndarray = None,
              is4d: bool = None) -> dict:
    """
    Scan all keys of the SEG-Y file. This is useful for unsorted SEG-Y file.
    """
    if isinstance(fname, Pysegy):
        segy = fname
    else:
        segy = Pysegy(str(fname))

    if keys is None:
        if iline is None or xline is None:
            raise ValueError("keys is None, so iline and xline must be inputed, but got None")
        keys = segy.get_trace_keys([iline, xline, offset], [4] * 3, 0, segy.ntrace)
    else:
        assert keys.ndim == 2
        if keys.shape[1] == 2:
            is4d = False 
        elif keys.shape[1] == 3:
            is4d = True
        else:
            raise ValueError("keys' shape must be (N, 2) or (N, 3)")

    ib = keys[:, 0].min()
    ie = keys[:, 0].max()
    diff = np.diff(np.sort(keys[:, 0]))
    diff = diff[diff != 0]
    istep = diff.min()
    if (ie - ib) % istep != 0:
        raise RuntimeError("can not create geomtry (error when determine iline/istep)") # yapf: disable

    xb = keys[:, 1].min()
    xe = keys[:, 1].max()
    diff = np.diff(np.sort(keys[:, 1]))
    diff = diff[diff != 0]
    xstep = diff.min()
    if (xe - xb) % xstep != 0:
        raise RuntimeError("can not create geomtry (error when determine xline/xstep)") # yapf: disable

    if is4d is None or is4d == True:
        is4d = True
        ob = keys[:, 2].min()
        oe = keys[:, 2].max()
        if oe == ob:
            is4d = False
        elif np.unique(keys[:, 2]).size > 500:
            is4d = False
        else:
            try:
                diff = np.diff(np.sort(keys[:, 2]))
                diff = diff[diff != 0]
                ostep = diff.min()
                if (oe - ob) % ostep != 0:
                    raise RuntimeError("can not create geomtry (error when determine xline/xstep)") # yapf: disable

                no = int((oe - ob) // ostep + 1)
            except:
                is4d = False

    ni = int((ie - ib) // istep + 1)
    nx = int((xe - xb) // xstep + 1)
    nt = segy.nt
    location = [iline, xline, offset] if is4d else [iline, xline]
    shape = [ni, nx, no, nt] if is4d else [ni, nx, nt]
    ir = dict(min_iline=ib, max_iline=ie, istep=istep)
    xr = dict(min_xline=xb, max_xline=xe, xstep=xstep)

    keys[:, 0] = (keys[:, 0] - ib) / istep
    keys[:, 1] = (keys[:, 1] - xb) / xstep
    if is4d:
        keys[:, 2] = (keys[:, 2] - ob) / ostep
    keys = np.round(keys).astype(np.int32)

    geom = np.full(shape[:-1], -1, np.int32)
    if is4d:
        geom[keys[:, 0], keys[:, 1], keys[:, 2]] = np.arange(segy.ntrace)
    else:
        geom[keys[:, 0], keys[:, 1]] = np.arange(segy.ntrace)

    geominfo = {
        'location': location,
        'shape': shape,
        'iline': ir,
        'xline': xr,
    }
    if is4d:
        geominfo['offset'] = dict(min_offset=ob, max_offset=oe, ostep=ostep)

    geominfo['geom'] = geom

    return geominfo


def load_by_geom(
    fname,
    geominfo: dict,
    *,
    ib: int = 0,
    ie: int = -1,
    xb: int = 0,
    xe: int = -1,
    ob: int = 0,
    oe: int = -1,
    tb: int = 0,
    te: int = -1,
) -> np.ndarray:
    """
    Using the given geom, extract data without scanning. This is useful for unsorted SEG-Y file.
    
    `geominfo` must contain the key locations, shape, and ranges. For example:
    
    geominfo = {
        'location': [189, 193],
        'shape': [650, 781, 951],
        'iline': dict(min_iline=2201, max_iline=2850, istep=1),
        'xline': dict(min_xline=5650, max_xline=7210, xstep=2),
        'geom': geom, # geom is a (ni, nx) array for 3D or (ni, nx, no) for 4D}
    """
    is4d = False
    shape = geominfo['shape']
    if 'offset' in geominfo:
        is4d = True
        assert len(shape) == 4
    else:
        assert len(shape) == 3

    ie = shape[0] if ie == -1 else ie
    xe = shape[1] if xe == -1 else xe
    te = shape[3] if te == -1 else te
    if is4d:
        oe = shape[2] if oe == -1 else oe

    # check bound
    assert ib >= 0 and ib < ie and ie <= shape[0]
    assert xb >= 0 and xb < xe and xe <= shape[1]
    assert tb >= 0 and tb < te and te <= shape[3]
    if is4d:
        assert ob >= 0 and ob < oe and oe <= shape[2]

    if isinstance(fname, Pysegy):
        segy = fname
    else:
        segy = Pysegy(str(fname))

    if not is4d:
        x, y = np.meshgrid(np.arange(ib, ie), np.arange(xb, xe), indexing='ij')
        shape = x.shape
        index = geominfo['geom'][x.flatten(), y.flatten()]
        d = segy.collect(index, tb, te).reshape(*shape, -1)
    else:
        x, y, z = np.meshgrid(np.arange(ib, ie),
                              np.arange(xb, xe),
                              np.arange(ob, oe),
                              indexing='ij')
        shape = x.shape
        index = geominfo['geom'][x.flatten(), y.flatten(), z.flatten()]
        d = segy.collect(index, tb, te).reshape(*shape, -1)

    return d


################# atomic operation ####################

def set_keyi2(header: np.ndarray, loc: int, value: int, endian: str = '>'):
    loc = loc - 1
    header[loc:loc + 2] = np.array([value], dtype=f'{endian}i2').view(np.uint8)
    return header

def set_keyi4(header: np.ndarray, loc: int, value: int, endian: str = '>'):
    loc = loc - 1
    header[loc:loc + 4] = np.array([value], dtype=f'{endian}i4').view(np.uint8)
    return header

def get_keyi2(header: np.array, loc: int, endian: str = '>'):
    loc = loc - 1
    return np.frombuffer(header[loc:loc + 2], dtype=f'{endian}i2')[0]

def get_keyi4(header: np.array, loc: int, endian: str = '>'):
    loc = loc - 1
    return np.frombuffer(header[loc:loc + 4], dtype=f'{endian}i4')[0]



############### Deprecated functions ####################


def fromfile_by_guess(*args, **kwargs):
    raise RuntimeError(
        "`fromfile_by_guess` is deprecated and removed. Please use `fromfile` instead."
    )


def tofile_by_guess(*args, **kwargs):
    raise RuntimeError(
        "`tofile_by_guess` is deprecated and removed. Please use `tofile` instead."
    )


def scan_prestack(*args, **kwargs):
    raise RuntimeError(
        "`scan_prestack` is deprecated and removed. Please use `metaInfo` instead."
    )


def load_prestack3D(*args, **kwargs):
    raise RuntimeError(
        "`load_prestack3D` is deprecated and removed. Please use `fromfile` instead."
    )


def create_by_sharing_header_guess(*args, **kwargs):
    raise RuntimeError(
        "`create_by_sharing_header_guess` is deprecated and removed. Please use `create_by_sharing_header` instead."
    )


def scan_unsorted3D(*args, **kwargs):
    """
    Scan an unsored 3D SEG-Y file and get the geometry
    """
    raise RuntimeError(
        "`scan_unsorted3D` is deprecated and removed. Please use `full_scan` instead."
    )


def load_unsorted3D(*args, **kwargs):
    """
    using fromfile_without_scan to load unsorted 3D segy file
    """
    raise RuntimeError(
        "`load_unsorted3D` is deprecated and removed. Please use `load_by_geom` instead."
    )
