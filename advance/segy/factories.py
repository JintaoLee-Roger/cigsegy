# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

import warnings
import numpy as np
from segy.cpp import _CXX_SEGY


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
    segy = _CXX_SEGY.SegyRWpy(str(segy_name))
    print(segy.textual_header(coding))
    segy.close()


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
) -> None:
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
    """


def fromfile(
    segy_name: str,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
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


def collect(
    segy_in: str,
    beg: int = -1,
    end: int = 0,
    tbeg: int = -1,
    tend: int = 0,
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

    Returns
    -------
    numpy.ndarray :
        its shape = (trace_count, n-time)
    """
    segy = _CXX_SEGY.SegyRWpy(str(segy_in))
    if beg < 0:
        beg = 0
        end = segy.ntrace
    if end == 0:
        end = beg + 1
    if end < 0:
        end = segy.ntrace

    if beg > end or end >= segy.ntrace:
        raise ValueError(f"beg = {beg}, end = {end} is out of range.")

    if tbeg < 0:
        tbeg = 0
        tend = segy.nt
    if tend == 0:
        tend = tbeg + 1
    if tend < 0:
        tend = segy.nt
    if tbeg > tend or tend >= segy.nt:
        raise ValueError(f"tbeg = {tbeg}, tend = {tend} is out of range.")

    data = segy.collect(beg, end, tbeg, tend).squeeze()
    return data


def tofile(
    segy_name: str,
    out_name: str,
    **kwargs,
) -> None:
    """
    convert a segy file to a binary file

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
    as_2d : bool
        if True, just remove the header and convert data to IEEE 32 in litte endian
    """


def create_by_sharing_header(
    segy_name: str,
    header_segy: str,
    src,
    **kwargs,
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
    """


def get_trace_keys(segy,
                   keyloc,
                   beg: int = -1,
                   end: int = 0,
                   force: int = None) -> np.ndarray:
    """
    get values at key location of trace headers

    Parameters
    ----------
    segy : str or Pysegy
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

    Returns
    -------
    np.ndarray
        shape as (end-beg, )
    """


def modify_bin_key(segy_name: str,
                   loc: int,
                   value,
                   force: bool = False,
                   type: str = None) -> None:
    """
    modify the value of the binary header key.

    Parameters
    -----------
    segy_name : str
        segy file name
    loc : int
        location of the binary key
    value : float or int
        assigned value to the key
    force : bool
        force to write
    type : str
        value type of the assigned value when force is True, can be
        one of {'int8', 'int16', 'int32', 'int64', 'float32', 'float64'}
    """


def modify_trace_key(segy_name: str,
                     loc: int,
                     value,
                     idx: int,
                     force: bool = False,
                     type: str = None) -> None:
    """
    modify the value of the trace header key.

    Parameters
    -----------
    segy_name : str
        segy file name
    loc : int
        location of the binary key
    value : float or int
        assigned value to the key
    idx : int
        trace index, if idx < 0, assign the value for all traces.
    force : bool
        force to write
    type : str
        value type of the assigned value when force is True, can be
        one of {'int8', 'int16', 'int32', 'int64', 'float32', 'float64'}
    """


def ibm_to_ieee(value, is_big_endian):
    """
    convert IBM floating point to IEEE floating point
    """
    return _CXX_SEGY.ibm_to_ieee(value, is_big_endian)


def ieee_to_ibm(value, is_little_endian):
    """
    convert IEEE floating point to IBM floating point
    """
    return _CXX_SEGY.ieee_to_ibm(value, is_little_endian)


def ibms_to_ieees(ibms: np.ndarray, is_big_endian: bool) -> np.ndarray:
    """
    convert IBM floating array to IEEE floating points
    """
    if ibms.dtype != np.float32:
        warnings.warn(
            f"ibms.dtype = {ibms.dtype} is not np.float32, it will be converted to np.float32."
        )
        ibms = ibms.astype(np.float32)
    return _CXX_SEGY.ieees_to_ibms(ibms, is_big_endian)


def ieees_to_ibms(ieees: np.ndarray, is_little_endian: bool) -> np.ndarray:
    """
    convert IEEE floating array to IBM floating points
    """
    if ieees.dtype != np.float32:
        warnings.warn(
            f"ieees.dtype = {ieees.dtype} is not np.float32, it will be converted to np.float32."
        )
        ieees = ieees.astype(np.float32)
    return _CXX_SEGY.ieees_to_ibms(ieees, is_little_endian)
