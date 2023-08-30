# Copyright (c) 2023 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
#
# github: https://github.com/JintaoLee-Roger

import warnings
import numpy as np
from typing import List, Tuple, Dict
from .cigsegy import (Pysegy, fromfile, tofile, create_by_sharing_header,
                      _load_prestack3D, kBinaryHeaderHelp, kTraceHeaderHelp)
from . import utils


def create(segy_out: str,
           binary_in: str or np.ndarray,
           shape: Tuple = None,
           format: int = 5,
           dt: int = 2000,
           start_time: int = 0,
           iline_interval: float = 25,
           xline_interval: float = 25,
           min_iline: int = 1,
           min_xline: int = 1,
           custom_info: List[str] = []) -> None:
    """
    Create a segy format file from a binary file or np.ndarray
    
    Parameters
    ----------
    segy_out : str
        out segy format file path
    binary_in : str or np.array
        the input binary file or array
    shape : Tuple
        len == 3
    format : int
        the data format code, 1 for 4 bytes IBM float, 5 for 4 bytes IEEE float
    dt : int
        data sample interval, 2000 means 2ms
    start_time : int
        start time for each trace
    iline_interval : int
        inline interval, will affect cdp x and cdp y
    xline_interval : int
        crossline interval, will affect cdp x and cdp y
    min_iline : int
        the start inline number
    min_xline : int 
        the start crossline number
    custom_info : List[str]
        textual header info by user custom, max: 12 rows each row is less than 76 chars
    """
    if isinstance(binary_in, str):
        assert shape is not None
        assert len(shape) == 3
        segy_create = Pysegy(binary_in, shape[2], shape[1], shape[0])
    elif isinstance(binary_in, np.ndarray):
        assert len(binary_in.shape) == 3
        sizeZ, sizeY, sizeX = binary_in.shape
        segy_create = Pysegy(sizeX, sizeY, sizeZ)
    else:
        raise ValueError(
            f'the input argument: binary_in must be a string or np array')
    segy_create.setDataFormatCode(format)
    segy_create.setSampleInterval(dt)
    segy_create.setStartTime(start_time)
    segy_create.setInlineInterval(iline_interval)
    segy_create.setCrosslineInterval(xline_interval)
    segy_create.setMinInline(min_iline)
    segy_create.setMinCrossline(min_xline)
    if isinstance(binary_in, str):
        segy_create.create(segy_out, custom_info)
    else:
        segy_create.create(segy_out, binary_in, custom_info)


def textual_header(segy_name: str, coding: str = None) -> None:
    """
    Print segy file's 3200 bytes textual header.

    Parameters
    ----------
    segy_name : str
        input segy file
    coding : str
        force the coding as 'a': ascii or 'e': ebcdic. If coding is None, cigsegy will guess the coding
    """
    if coding is not None:
        if coding != 'a' and coding != 'e' and coding != 'u':
            raise ValueError(
                f"`coding` can be one of " +
                f"{{ 'a': 'ascii', 'e': 'EBCDIC', 'u': 'Unkown' }}, " +
                f"but your input is `coding='{coding}'`")

    coding = 'u' if coding is None else coding
    segy = Pysegy(segy_name)
    print(segy.textual_header(coding))
    segy.close_file()


def metaInfo(segy_name: str,
             iline: int = 189,
             xline: int = 193,
             istep: int = 1,
             xstep: int = 1,
             xloc: int = 73,
             yloc: int = 77,
             use_guess: bool = False) -> None:
    """
    print meta info of `segy_name` file

    Parameters
    ----------
    iline : int
        iline location in trace header
    iline : int
        iline location in trace header
    istep : int
        iline step, 2 means iline is like 100, 102, 104, ...
    xstep : int
        xline step
    xloc : int
        cdp x (real world) value location in trace header
    yloc : int
        cdp y (real world) value location in trace header
    use_guess : bool
        if iline/xline/istep/xstep are unknow, you can set use_guess = True to guess them
    """
    if use_guess:
        [iline, xline, istep, xstep] = utils.guess(segy_name)[0]
    segy = Pysegy(segy_name)
    segy.setInlineLocation(iline)
    segy.setCrosslineLocation(xline)
    segy.setSteps(istep, xstep)
    segy.setXLocation(xloc)
    segy.setYLocation(yloc)
    print(segy.metaInfo())
    segy.close_file()


def fromfile_by_guess(segy_name: str) -> np.ndarray:
    """
    reading from a segy file.

    Parameters
    ----------
    segy_name : str
        the input segy file name

    Returns
    -------
    np.ndarray
        3D array data
    """

    loc = utils.guess(segy_name)

    for l in loc:
        try:
            metaInfo(segy_name, l[0], l[1], l[2], l[3])
            d = fromfile(segy_name, l[0], l[1], l[2], l[3])
            return d
        except:
            continue

    raise RuntimeError(
        "Cannot read by guess location, please specify the location")


def tofile_by_guess(segy_name: str, out_name: str) -> None:
    """
    convert a segy file to a binary file

    Parameters
    -----------
    segy_name : str
        the input segy file name
    out_name : str
        the output binary file name
    """
    loc = utils.guess(segy_name)
    finish = False

    for l in loc:
        try:
            metaInfo(segy_name, l[0], l[1], l[2], l[3])
            tofile(out_name, l[0], l[1], l[2], l[3])
            finish = True
            break
        except:
            continue

    if not finish:
        raise RuntimeError(
            "Cannot read by guess location, please specify the location")


def create_by_sharing_header_guess(segy_name: str,
                                   header_segy: str,
                                   src: np.ndarray or str,
                                   shape: list or tuple = None,
                                   offset: list or tuple or dict = None,
                                   custom_info: List[str] = []) -> None:
    """
    create a segy and its header is from an existed segy.

    Parameters
    ----------
    segy_name : str
        the out segy name
    header_segy : str 
        the header segy file
    src : np.ndarray
        source data
    shape : Tuple or List
        if src is str, shape must be specify
    offset : Tuple or List
        the offset of a sub data from a original data, e.g., dsub = d[256:400, 500: 100:], offset = [256, 500, 100]
    custom_info : List[str]
        textual header info by user custom, max: 12 rows each row is less than 76 chars, use it when offset is not None
    """
    if isinstance(src, str) and shape is None:
        raise ValueError("Shape is None!")

    loc = utils.guess(header_segy)
    finish = False

    for l in loc:
        try:
            if isinstance(src, str):
                create_by_sharing_header(segy_name,
                                         header_segy,
                                         src,
                                         shape,
                                         l[0],
                                         l[1],
                                         l[2],
                                         l[3],
                                         offset=offset,
                                         custom_info=custom_info)
            else:
                create_by_sharing_header(segy_name,
                                         header_segy,
                                         src,
                                         l[0],
                                         l[1],
                                         l[2],
                                         l[3],
                                         offset=offset,
                                         custom_info=custom_info)
            finish = True
            break
        except:
            continue

    if not finish:
        raise RuntimeError(
            "Cannot read by guess location, please specify the location")


def read_header(segy: str, type, n=0, printstr=True):
    """
    Read binary or trace header

    Parameters
    ----------
    segy : str
        input segy file
    type: str
        can be one of ['bh', 'th', 't'],
            'bt' means binary header, 'th' means trace header,
            't' means trace (include trace header and trace)
    n : int
        trace number when type is 'th' or t
    printstr : bool
        print header information, if False, return a dict of header's infomation

    Returns
    -------
    Dict or None
    """
    segy = Pysegy(segy)

    if type == 'bh':
        arr = segy.get_binary_header()
        out = utils.convert_header(arr, kBinaryHeaderHelp, printstr, 3200)
        if printstr:
            for l in out:
                print(l)
            return
        return out
    elif type == 'th':
        arr = segy.get_trace_header(n)
        out = utils.convert_header(arr, kTraceHeaderHelp, printstr)
        if printstr:
            for l in out:
                print(l)
            return
        return out
    elif type == 't':
        arr = segy.get_trace(n)
        out = utils.convert_header(arr[:240], kTraceHeaderHelp, printstr)
        d = utils.convert_trace(arr[240:])
        if printstr:
            for l in out:
                print(l)
            return d
        return out, d


def get_metaInfo(segy_name: str,
                 iline: int = 189,
                 xline: int = 193,
                 istep: int = 1,
                 xstep: int = 1,
                 xloc: int = 73,
                 yloc: int = 77,
                 use_guess: bool = False,
                 apply_scalar: bool = False) -> Dict:
    """
    get metainfo dict of `segy_name` file

    Parameters
    ----------
    segy_name : str
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
    use_guess : bool
        if iline/xline/istep/xstep are unknow, 
            you can set use_guess = True to guess them
    apply_scalar : bool
        apply scalar to 'i-interval' and 'x-interval'

    Returns
    -------
    Dict
        Dict of meta information 
    """
    if use_guess:
        [iline, xline, istep, xstep] = utils.guess(segy_name)[0]
    segy = Pysegy(segy_name)
    segy.setInlineLocation(iline)
    segy.setCrosslineLocation(xline)
    segy.setSteps(istep, xstep)
    segy.setXLocation(xloc)
    segy.setYLocation(yloc)
    segy.scan()
    m = segy.get_metaInfo()
    segy.close_file()

    return utils.metainfo_to_dict(m, apply_scalar)


def trace_count(segy: str or Pysegy) -> int:
    """
    Count the total numbers of a segy file

    Parameters
    ----------
    segy: str or Pysegy
        input segy file

    Returns
    -------
    int
        The total numbers of a segy file
    """
    if isinstance(segy, str):
        segy = Pysegy(segy)
        count = segy.trace_count
        segy.close_file()
        return count

    return segy.trace_count


def scan_prestack(segy: str or Pysegy,
                  iline: int,
                  xline: int,
                  offset: int = 37) -> Dict:
    """
    scan the 3D prestack gather and get the geometry of the SEG-Y file

    Parameters
    ----------
    segy : str or Pysegy
        input segy file
    iline : int
    xline : int
    offset : int

    Returns
    --------
    geom : Dict
    """
    if isinstance(segy, str):
        segyc = Pysegy(segy)

    if segyc.get_metaInfo().trace_sorting_code == 4:
        warnings.warn("trace sorting code is 4, this means the " +
                      "segy is a horizontally stacked (post-stack) data." +
                      "`scan_prestack` may not be correct.")

    tracecount = segyc.trace_count
    i0 = utils.get_trace_keys(segyc, iline, 0)
    ie = utils.get_trace_keys(segyc, iline, tracecount - 1)
    ic = utils.get_trace_keys(segyc, iline, tracecount // 2)
    # tracecount = int(segyc.trace_count * ratio)

    if not (i0 < ie and ic < ie):
        raise RuntimeError(
            f"cannot analyse this segy, " +
            f"as inline numbers in 0-th, center-th, end-th traces are " +
            f"{i0}, {ic}, {ie}. Inline numbers should be in ascending order.\n"
            + "the issue can arise due to the following reasons:\n" +
            "1. `iline` is incorrect, set correct `iline` to avoid this error\n"
            + "2. The SEG-Y file is unsorted\n" +
            "3. This is a 2D SEG-Y file (only contains one line)")

    pos = [iline, xline, offset]

    def _scan(il, xl, of, beg=0):
        out = []
        # scan inline
        diff_il = np.diff(il)
        error_idx = np.where(diff_il < 0)[0]  # index of decrease elements
        if len(error_idx) > 0:
            raise RuntimeError(
                "inline numbers must be incremental, but in trace " +
                f"{error_idx[0]+beg}->{error_idx[0]+1+beg}, " +
                f"inlines are " +
                f"{il[error_idx[0]]}->{il[error_idx[0]+1]}." +
                "This SEG-Y file may be unsorted.")

        sep_il = np.where(diff_il > 0)[0] + 1  # index of increase elements
        unique_il = np.unique(il)  # unique elements of inline
        diff_uil = np.diff(unique_il)
        istep = diff_uil[0]
        assert np.all(diff_uil == istep), f"inline steps must be same"
        out += [istep]

        # scan crossline
        diff_xl = np.diff(xl)
        error_xl = np.where(diff_xl < 0)[0]  # decrease elements within lines
        # exclude sep of inline
        error_xl = error_xl[~np.isin(error_xl, sep_il - 1)]
        if len(error_xl) > 0:
            raise RuntimeError(
                "crossline numbers must be incremental within an inline, " +
                f"but in trace {error_xl[0]+beg}->{error_xl[0]+1+beg}" +
                f", inline: {il[error_xl[0]]}->{il[error_xl[0]+1]}, "
                f"crossline: {xl[error_xl[0]]}->{xl[error_xl[0]+1]}")

        sep_xl = np.where(diff_xl != 0)[0] + 1  # index of sep
        unique_xl = np.unique(xl)
        xstep = np.diff(unique_xl).min()
        x0 = unique_xl[0]
        xe = unique_xl[-1]
        # nx = (xe - x0) / xstep + 1
        out += [x0, xe, xstep]

        # scan offset
        diff_of = np.diff(of)
        # decrease elements within crossline
        # offset must be strongly incremental within a crossline
        error_of = np.where(diff_of < 1)[0]  # < 1 for stongly
        error_of = error_of[~np.isin(error_of, sep_xl - 1)]
        if len(error_of) > 0:
            raise RuntimeWarning(
                "offset must be strongly incremental within a crossline, " +
                f"but in trace {error_of[0]+beg}->{error_of[0]+1+beg}, " +
                f"crossline: {xl[error_of[0]]}->{xl[error_of[0]+1]}, " +
                f"offset: {of[error_of[0]]}->{of[error_of[0]+1]}")

        # sep_of = np.where(diff_of < 0)[0] + 1  # sep in decreased pos
        unique_of = np.unique(of)
        ostep = np.diff(unique_of).min()
        o0 = unique_of[0]
        oe = unique_of[-1]
        # no = (oe - o0) / ostep + 1
        out += [o0, oe, ostep]

        return out

    # out = [istep, x0, xe, xstep, o0, oe, ostep]
    per_read = 100000

    range_func = range
    try:
        from tqdm import trange
        range_func = trange
    except:
        pass

    if tracecount < per_read * 2:
        keys = utils.get_trace_keys(segyc, pos, 0, tracecount)
        out = _scan(keys[:, 0], keys[:, 1], keys[:, 2])
        istep, x0, xe, xstep, o0, oe, ostep = out

    else:
        N = tracecount // per_read
        istep, x0, xe, xstep, o0, oe, ostep = [0] * 7

        for i in range_func(N):
            beg = i * per_read
            end = (i + 1) * per_read if i < N - 1 else tracecount
            ext = end + 100 if i < N - 1 else end
            keys = utils.get_trace_keys(segyc, pos, beg, ext)
            out = _scan(keys[:, 0], keys[:, 1], keys[:, 2], beg)

            if i == 0:
                istep, x0, xe, xstep, o0, oe, ostep = out
            else:
                istep = min(istep, out[0])
                x0 = min(x0, out[1])
                xe = max(xe, out[2])
                xstep = min(xstep, out[3])
                o0 = min(o0, out[4])
                oe = max(oe, out[5])
                ostep = min(ostep, out[6])

    ni = (ie - i0) // istep + 1
    nx = (xe - x0) // xstep + 1
    no = (oe - o0) // ostep + 1
    nt = segyc.get_metaInfo().sizeX

    geom = {
        'shape': [ni, nx, no, nt],
        'iline': dict(min_iline=i0, max_iline=ie, istep=istep),
        'xline': dict(min_xline=x0, max_xline=xe, xstep=xstep),
        'offset': dict(min_offset=o0, max_offset=oe, ostep=ostep),
    }

    if isinstance(segy, str):
        segyc.close_file()

    return geom


def load_prestack3D(segy: str,
                    iline: int,
                    xline: int,
                    offset: int = 37,
                    geom: Dict = None,
                    fill: float = 0) -> np.ndarray:
    """
    Load 3D prestack gather (4D array)

    Parameters
    ----------
    segy : str
        input segy file
    iline : int
        iline location
    xline : int
        crossline location
    offset : int
        offset location
    geom : Dict
        geometry of the segy file, it can be obtained by `scan_prestack` function,
        if None, will call `scan_prestack` to get geom
    fill : float
        fill value for missing traces

    Returns
    -------
    np.ndarray
        shape as (ni, nx, no, nt), no means number of offset
    """
    if geom is None:
        print('Scan 3D prestack:')
        geom = scan_prestack(segy, iline, xline, offset)

    shape = geom['shape']
    iline_info = geom['iline']
    i0, ie, istep = iline_info['min_iline'], iline_info[
        'max_iline'], iline_info['istep']
    xline_info = geom['xline']
    x0, xe, xstep = xline_info['min_xline'], xline_info[
        'max_xline'], xline_info['xstep']
    offset_info = geom['offset']
    o0, oe, ostep = offset_info['min_offset'], offset_info[
        'max_offset'], offset_info['ostep']

    print('Load 3D prestack:')
    out = _load_prestack3D(segy, shape, i0, ie, x0, xe, o0, oe, istep, xstep,
                           ostep, iline, xline, offset, fill)

    return out
