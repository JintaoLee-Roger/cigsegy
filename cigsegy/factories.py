# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
"""
Core functions for easy and direct access to common SEG-Y file operations.

This module provides a set of factory functions that offer a simplified interface for working with SEG-Y files. 
These functions are designed for quick and direct usage, allowing you to perform essential operations 
like header reading, file scanning, data retrieval, and file creation using a convenient `cigsegy.xxx` syntax.
"""

import gc
import os
import warnings
import numpy as np
from cigsegy.cpp import _CXX_SEGY
from cigsegy.constinfo import kBinaryHeaderHelp, kTraceHeaderHelp
from cigsegy.tools import get_metaInfo
from cigsegy import utils, createtool


def _normalize_textual_header(meta: dict, textual):
    if textual in ("", None):
        return ""
    if isinstance(textual, (bytes, bytearray)):
        textual = bytes(textual)
        if len(textual) != 3200:
            raise ValueError(
                f"textual header length must be exactly 3200 bytes, got {len(textual)}"
            )
        return textual
    if isinstance(textual, str) and len(textual) == 3200:
        return textual
    return createtool.generate_textual(meta, textual)


def textual_header(segy_name: str,
                   coding: str = None,
                   *,
                   printtext: bool = True):
    """
    Print segy file's 3200 bytes textual header by default.
    Set `printtext=False` to return the textual header string instead.

    Parameters
    ----------
    segy_name : str
        input segy file
    coding : str
        force the coding as 'a': ascii or 'e': ebcdic. If coding is None, cigsegy will guess the coding
    printtext : bool
        if True, print the textual header and return None. If False, return the textual header string.
    """
    if coding is not None:
        if coding != 'a' and coding != 'e' and coding != 'u':
            raise ValueError(f"`coding` can be one of {{ 'a': 'ascii', 'e': 'EBCDIC', 'u': 'Unkown' }}, but your input is `coding='{coding}'`") # yapf: disable

    coding = 'u' if coding is None else coding
    if isinstance(segy_name, _CXX_SEGY.Pysegy):
        segy = segy_name
    else:
        segy = _CXX_SEGY.Pysegy(str(segy_name))

    text = segy.textual_header(coding)

    if not isinstance(segy_name, _CXX_SEGY.Pysegy):
        segy.close()

    if printtext:
        print(text)
        return
    return text


def metaInfo(
    segy_name: str,
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
    apply_scalar: bool = True,
    printtext: bool = True,
):
    """
    Print meta info of `segy_name` by default.
    Set `printtext=False` to return the formatted meta-info string instead.

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
    printtext : bool
        if True, print meta info and return None. If False, return the formatted meta-info string.
    """

    meta = get_metaInfo(segy_name, iline, xline, offset, istep, xstep, ostep, xloc, yloc, is4d=is4d, apply_scalar=apply_scalar) # yapf: disable
    out = utils.parse_metainfo(meta)
    if printtext:
        print(out)
        return
    return out


def fromfile(
    segy_name: str,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
    *,
    is4d: bool = None,
) -> np.ndarray:
    """
    reading from a segy file.

    Parameters
    ----------
    segy_name : str
        the input segy file name
    iline : int
        the inline number field in each trace header
    xline : int
        the crossline number field in each trace header
    istep : int
        the step of inline numbers
    xstep : int
        the step of crossline numbers
    
    Returns
    -------
    numpy.ndarray :
        shape as (n-inline, n-crossline, n-time)
    """
    [iline, xline, offset, istep, xstep, ostep, xloc, yloc, _is4d] = utils.guess(segy_name, iline, xline, offset, istep, xstep, ostep, 181, 185) # yapf: disable

    if is4d is None:
        is4d = _is4d

    if isinstance(segy_name, _CXX_SEGY.Pysegy):
        segy = segy_name
    else:
        segy = _CXX_SEGY.Pysegy(str(segy_name))

    segy.setLocations(iline, xline, offset)
    segy.setSteps(istep, xstep, ostep)
    ndim = 4 if is4d else 3
    segy.set_segy_type(ndim)
    segy.scan()
    d = segy.read()

    if not isinstance(segy_name, _CXX_SEGY.Pysegy):
        segy.close()
    return d


def collect(
    segy_in: str,
    beg: int = -1,
    end: int = 0,
    tbeg: int = -1,
    tend: int = 0,
    indices: np.ndarray = None,
) -> np.ndarray:
    """
    collect traces as a 2D data from the `segy_in` file in
    range of [beg, end), beg < 0 means collect all traces,
    end < 0 means collect traces from beg to the last trace,
    end == 0 means read the beg-th trace (one trace).

    Parameters
    ----------
    segy_in : str
        the input segy file
    beg : int
        the begin index of traces, < 0 means collect all traces
    end : int
        the end index of traces (not include), < 0 means collect 
            traces from beg to the last trace, == 0 means read 
            the beg-th trace (one trace).
    tbeg : int
        the start time index
    tend : int
        the end time index (not include)
    indices : np.ndarray
        the indices of traces, if `indices` is not None, `beg`, `end` will be ignored.

    Returns
    -------
    numpy.ndarray :
        its shape = (trace_count, n-time)
    """
    if isinstance(segy_in, _CXX_SEGY.Pysegy):
        segy = segy_in
    else:
        segy = _CXX_SEGY.Pysegy(str(segy_in))

    if indices is not None:
        if beg != -1 or end != 0:
            print("indices is not None, beg and end will be ignored.")
        indices = np.array(indices)
        assert np.issubdtype(indices.dtype, np.integer), "indices must be an integer type"
        indices = indices.astype(np.int32)
        assert indices.ndim == 1, "indices must be a 1D array"
    else:
        if beg < 0:
            beg = 0
            end = segy.ntrace
        if end == 0:
            end = beg + 1
        if end < 0:
            end = segy.ntrace

    if beg > end or end > segy.ntrace:
        raise ValueError(f"beg = {beg}, end = {end} is out of range.")

    if tbeg < 0:
        tbeg = 0
        tend = segy.nt
    if tend == 0:
        tend = tbeg + 1
    if tend < 0:
        tend = segy.nt
    if tbeg > tend or tend > segy.nt:
        raise ValueError(f"tbeg = {tbeg}, tend = {tend} is out of range.")

    if indices is None:
        data = segy.collect(beg, end, tbeg, tend).squeeze()
    else:
        data = segy.collect(indices, tbeg, tend).squeeze()

    if not isinstance(segy_in, _CXX_SEGY.Pysegy):
        segy.close()

    return data


def _prepare_segy_for_stream_export(
    segy_name,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
    *,
    is4d: bool = None,
    as2d: bool = False,
):
    if isinstance(segy_name, _CXX_SEGY.Pysegy):
        segy = segy_name
        need_close = False
    else:
        segy = _CXX_SEGY.Pysegy(str(segy_name))
        need_close = True

    try:
        if as2d:
            return segy, need_close, (int(segy.ntrace), int(segy.nt))

        [iline, xline, offset, istep, xstep, ostep, xloc, yloc, _is4d] = utils.guess(segy, iline, xline, offset, istep, xstep, ostep, 181, 185) # yapf: disable
        if is4d is None:
            is4d = _is4d
        segy.setLocations(iline, xline, offset)
        segy.setSteps(istep, xstep, ostep)
        ndim = 4 if is4d else 3
        segy.set_segy_type(ndim)
        segy.scan()
        return segy, need_close, tuple(int(v) for v in segy.shape)
    except Exception:
        if need_close:
            segy.close()
        raise


def tofile(
    segy_name: str,
    out_name: str,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
    *,
    is4d: bool = None,
    as2d: bool = False,
) -> None:
    """
    Stream SEG-Y samples to a raw binary file.

    This function writes little-endian IEEE float32 samples without textual,
    binary, or trace headers.  It is useful when the SEG-Y is too large to load
    with ``fromfile`` but a downstream tool can consume raw float32 data.

    Parameters
    ----------
    segy_name : str
        the input segy file name
    out_name : str
        the output binary file name
    iline : int
        the inline number field in each trace header
    xline : int
        the crossline number field in each trace header
    istep : int
        the step of inline numbers
    xstep : int
        the step of crossline numbers
    as2d : bool
        if True, just remove the header and convert data to IEEE 32 in litte endian
    """
    out_name = os.fspath(out_name)
    segy, need_close, _ = _prepare_segy_for_stream_export(
        segy_name,
        iline,
        xline,
        offset,
        istep,
        xstep,
        ostep,
        is4d=is4d,
        as2d=as2d,
    )
    try:
        segy.tofile(out_name, as2d)
    finally:
        if need_close:
            segy.close()


def to_npy(
    segy_name: str,
    out_name: str,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
    *,
    is4d: bool = None,
    as2d: bool = False,
) -> tuple:
    """
    Stream SEG-Y samples directly to a NumPy ``.npy`` file.

    Unlike ``fromfile(...); np.save(...)``, this function does not materialize
    the whole SEG-Y volume in memory.  It writes a NumPy header first and then
    streams little-endian float32 samples into the array payload.

    Returns
    -------
    tuple
        The output array shape.
    """
    out_name = os.fspath(out_name)
    segy, need_close, shape = _prepare_segy_for_stream_export(
        segy_name,
        iline,
        xline,
        offset,
        istep,
        xstep,
        ostep,
        is4d=is4d,
        as2d=as2d,
    )
    try:
        mmap = np.lib.format.open_memmap(
            out_name,
            mode="w+",
            dtype=np.dtype("<f4"),
            shape=shape,
        )
        payload_offset = int(mmap.offset)
        mmap.flush()
        del mmap
        gc.collect()
        segy.tofile(out_name, as2d, payload_offset)
        return shape
    finally:
        if need_close:
            segy.close()


def create_by_sharing_header(
    out_segy: str,
    header_segy: str,
    src,
    shape: tuple = None,
    *,
    start: list = None,
    keylocs: list = None,
    is4d: bool = None,
    as2d: bool = False,
    textual: str = "",

    # the following params for freeing the t dimension
    strict: bool = True,
    start_time_from_zero: bool = False,
    dt_new: int = 0,
) -> None:
    """
    create a segy and its header is from an existed segy.

    Parameters
    ----------
    segy_name : str
        the out segy name
    header_segy : str
        the header segy file
    src : numpy.ndarray
        source data
    shape : tuple
        shape of the data, only used when src is a path of a binary file
    start : list
        start of the geometry
    keylocs : list
        key locations, [iline, xline, istep, xstep], or [iline, xline, offset, istep, xstep, ostep]
    is4d : bool
        if True, the data is 4D, otherwise 3D
    as2d : bool
        if True, treat the data as 2D data
    textual : str
        textual header
    strict : bool
        if True, the src data must be a subset of the header_segy file
    start_time_from_zero : bool
        if True, the start time (i.e., start[-1]) is from zero, otherwise it is from the header_segy file.
        This parameter is only used when strict is False.
    dt_new : int
        the new sample interval, only used when strict is False
    """
    segy = _CXX_SEGY.Pysegy(str(header_segy))
    try:
        if as2d:
            if start is None:
                start = [0, 0]
            meta = utils.make_2d_meta(segy)
            textual = _normalize_textual_header(meta, textual)
            if isinstance(src, np.ndarray):
                src = np.ascontiguousarray(src, dtype=np.float32)
                shape = list(src.shape)
                if src.ndim != 2:
                    raise ValueError("src.ndim must be 2 when as2d is True")
                segy.create_by_sharing_header(out_segy, src, start, as2d,
                                              textual, strict,
                                              start_time_from_zero, dt_new)
            else:
                if shape is None:
                    raise ValueError("`shape` must be provided when src is a file")
                if len(shape) != 2:
                    raise ValueError(
                        f"`shape` must be 2D when as2d is True, got {shape}"
                    )
                segy.create_by_sharing_header(out_segy, src, shape, start,
                                              as2d, textual, strict,
                                              start_time_from_zero, dt_new)
            return

        if keylocs is not None and len(keylocs) not in (4, 6):
            raise ValueError("`keylocs` must contain 4 or 6 items")

        _is4d = False
        if keylocs is None:
            iline, xline, offset, istep, xstep, ostep = [None] * 6
        elif len(keylocs) == 4:
            iline, xline, istep, xstep = keylocs
            offset, ostep = None, None
        else:
            iline, xline, offset, istep, xstep, ostep = keylocs

        [iline, xline, offset, istep, xstep, ostep, xloc, yloc, _is4d] = utils.guess(segy, iline, xline, offset, istep, xstep, ostep, 181, 185) # yapf: disable
        if is4d is None:
            is4d = _is4d
        segy.setLocations(iline, xline, offset)
        segy.setSteps(istep, xstep, ostep)
        ndim = 4 if is4d else 3
        segy.set_segy_type(ndim)
        segy.scan()

        if isinstance(src, np.ndarray):
            src = np.ascontiguousarray(src, dtype=np.float32)
            shape = list(src.shape)
        elif shape is None:
            raise ValueError("`shape` must be provided when src is a file")

        if segy.ndim != len(shape):
            raise ValueError(
                f"header ndim is {segy.ndim}, but input shape has {len(shape)} dims"
            )
        if start is None:
            start = [0] * len(shape)
        if len(start) != len(shape):
            raise ValueError("`start` length must match `shape` length")

        meta = {**segy.get_keylocs(), **segy.get_metainfo()}
        meta = utils.post_process_meta(segy, meta, False)
        textual = _normalize_textual_header(meta, textual)

        if isinstance(src, np.ndarray):
            segy.create_by_sharing_header(out_segy, src, start, as2d, textual,
                                          strict, start_time_from_zero,
                                          dt_new)
        else:
            segy.create_by_sharing_header(out_segy, src, shape, start, as2d,
                                          textual, strict,
                                          start_time_from_zero, dt_new)
    finally:
        segy.close()


def get_trace_keys(segyname,
                   keyloc,
                   beg: int = -1,
                   end: int = 0,
                   force: int = None,
                   indices: np.ndarray = None) -> np.ndarray:
    """
    get values at key location of trace headers

    Parameters
    ----------
    segyname : str or Pysegy
        input segy file
    keyloc : int or List or 1D np.ndarray
        key locations, can input multi-values
    beg : int
        begin trace, if beg < 0, means read all values from all traces.
    end : int
        end trace, if end < 0, means read values from beg to the last trace,
            if end == 0, means read one value from beg trace.
    force : int
        force to read a value in keyloc (size = force), even if the keyloc is 
        not in the location of the standard SEG-Y Trace header's keys.
        When force is not None, can only input one key location, i.e., keyloc 
        must be a int or one element List
    indices : np.ndarray
        the indices of traces, if `indices` is not None, `beg`, `end` will be ignored.

    Returns
    -------
    np.ndarray
        shape as (end-beg, )
    """
    if isinstance(segyname, _CXX_SEGY.Pysegy):
        segy = segyname
    else:
        segy = _CXX_SEGY.Pysegy(str(segyname))

    if indices is not None:
        if beg != -1 or end != 0:
            print("indices is not None, beg and end will be ignored.")
        indices = np.array(indices)
        assert np.issubdtype(indices.dtype, np.integer), "indices must be an integer type"
        indices = indices.astype(np.int32)
        assert indices.ndim == 1, "indices must be a 1D array"
    else:
        if beg < 0:
            beg = 0
            end = segy.ntrace
        if end == 0:
            end = beg + 1
        if end < 0:
            end = segy.ntrace

    if isinstance(keyloc, (int, np.integer)):
        keyloc = [keyloc]

    if force is not None:
        if isinstance(force, (int, np.integer)):
            force = [force] * len(keyloc)
        assert len(keyloc) == len(force), "force must have the same length as keyloc" # yapf: disable
        if indices is None:
            d = segy.get_trace_keys(keyloc, force, beg, end).squeeze()
        else:
            d = segy.get_trace_keys(keyloc, force, indices).squeeze()
    else:
        try:
            length = [kTraceHeaderHelp[k][1] for k in keyloc]
        except KeyError:
            raise ValueError(
                f"keyloc = {keyloc} is not a valid trace header key location. "
                "If you want to force to read a key, please input force, e.g., `get_trace_keys('t.sgy', 221, 10, force=4)`, `force` is the length of the key."
            )
        if indices is None:
            d = segy.get_trace_keys(keyloc, length, beg, end).squeeze()
        else:
            d = segy.get_trace_keys(keyloc, length, indices).squeeze()

    if d.size == 1 and d.ndim == 0:
        d = int(d)
    elif d.size == 1 and d.ndim == 1:
        d = int(d[0])

    if not isinstance(segyname, _CXX_SEGY.Pysegy):
        segy.close()
    return d


def modify_bin_key(segyname: str, loc: int, value, force: int = None) -> None:
    """
    modify the value of the binary header key.

    Parameters
    -----------
    segyname : str
        segy file name
    loc : int
        location of the binary key
    value : float or int
        assigned value to the key
    force : int
        one of 2, 4
    """
    if isinstance(segyname, _CXX_SEGY.Pysegy):
        segy = segyname
    else:
        segy = _CXX_SEGY.Pysegy(str(segyname), True)

    if force is not None:
        assert isinstance(force, (int, np.integer))
        l = force
    else:
        try:
            l = kBinaryHeaderHelp[loc][1]
        except KeyError:
            raise ValueError(
                f"loc = {loc} is not a valid binary header key location. "
                "If you want to force to modify a key, please input force, e.g., `modify_bin_key('t.sgy', 221, 10, force=4)`, `force` is the length of the key."
            )

    if l == 2:
        segy.set_bkeyi2(loc, value)
    elif l == 4:
        segy.set_bkeyi4(loc, value)
    else:
        raise ValueError("only support int32 (l==4) and int16 (l==2)")

    if not isinstance(segyname, _CXX_SEGY.Pysegy):
        segy.close()


def modify_trace_key(segyname: str,
                     loc: int,
                     value,
                     idx: int,
                     force: int = None) -> None:
    """
    modify the value of the trace header key.

    Parameters
    -----------
    segyname : str
        segy file name
    loc : int
        location of the binary key
    value : float or int
        assigned value to the key
    idx : int
        trace index, if idx < 0, assign the value for all traces.
    force : int
        one of 2/4
    """
    if isinstance(segyname, _CXX_SEGY.Pysegy):
        segy = segyname
    else:
        segy = _CXX_SEGY.Pysegy(str(segyname))

    if force is not None:
        assert isinstance(force, (int, np.integer))
        l = force
    else:
        try:
            l = kTraceHeaderHelp[loc][1]
        except KeyError:
            raise ValueError(
                f"loc = {loc} is not a valid trace header key location. "
                "If you want to force to modify a key, please input force, e.g., `modify_trace_key('t.sgy', 221, 10, 100000, force=4)`, `force` is the length of the key."
            )

    if l == 2:
        segy.set_keyi2(idx, loc, value)
    elif l == 4:
        segy.set_keyi4(idx, loc, value)
    else:
        raise ValueError("only support int32 (l==4) and int16 (l==2)")

    if not isinstance(segyname, _CXX_SEGY.Pysegy):
        segy.close()


def ibm_to_ieee(value, is_big_endian: bool):
    """
    convert IBM floating point to IEEE floating point
    """
    if isinstance(value, (float, int)):
        return _CXX_SEGY.ibm_to_ieee(value, is_big_endian)
    elif isinstance(value, np.ndarray):
        if value.dtype != np.float32:
            warnings.warn(f"value.dtype = {value.dtype} is not np.float32, it will be converted to np.float32.")  # yapf: disable
            value = value.astype(np.float32)
        return _CXX_SEGY.ibms_to_ieees(value, is_big_endian)
    else:
        raise ValueError(f"value = {value} is not a valid type, it should be one of float, int, np.ndarray") # yapf: disable


def ieee_to_ibm(value, is_little_endian_input: bool, is_big_endian_output: bool):
    """
    convert IEEE floating point to IBM floating point
    """
    if isinstance(value, (float, int)):
        return _CXX_SEGY.ieee_to_ibm(value, is_little_endian_input, is_big_endian_output)
    elif isinstance(value, np.ndarray):
        if value.dtype != np.float32:
            warnings.warn(f"value.dtype = {value.dtype} is not np.float32, it will be converted to np.float32.") # yapf: disable
            value = value.astype(np.float32)
        return _CXX_SEGY.ieees_to_ibms(value, is_little_endian_input, is_big_endian_output)
    else:
        raise ValueError(f"value = {value} is not a valid type, it should be one of float, int, np.ndarray") # yapf: disable


def create(
    segy_out: str,
    src: np.ndarray,
    format: int = 5,
    dt: int = 2000,
    start_time: int = 0,
    iline_interval: float = 25,
    xline_interval: float = 25,
    min_iline: int = 1,
    min_xline: int = 1,
    custom_info="",
) -> None:
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
    warnings.warn(
        "create function is deprecated and will remove in the future version, please consider use `cigsegy.SegyCreate` class instead.",
        DeprecationWarning,
    )
    if not isinstance(src, np.ndarray):
        raise TypeError("src only support np.ndarray")
    src = np.ascontiguousarray(src, dtype=np.float32)

    if src.ndim != 3:
        raise ValueError("this function only supports 3D volumes")
    start = [min_iline, min_xline, start_time]
    interval = [iline_interval, xline_interval, dt]

    meta = createtool.assemble_metainfo(src.shape, dformat=format, start=start, interval=interval)

    segyc = createtool.SegyCreate(segy_out, src, meta)
    segyc.set_textual_header(custom_info)
    segyc.create()
