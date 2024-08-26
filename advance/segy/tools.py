# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import Dict, List, Tuple
import numpy as np
from segy.cpp._CXX_SEGY import Pysegy
from segy import utils


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
    segy = Pysegy(fname)
    if type == 'bh':
        arr = segy.get_binary_header()
        out, hstring = utils.parse_bheader(arr)
    elif type == 'th':
        arr = segy.get_trace_header(n)
        out = utils.parse_theader(arr)
    else:
        raise ValueError("type must be one of ['bh', 'th']")

    if printstr:
        print('\n'.join(hstring))
        return
    return out


def get_metaInfo(
    segy_name: str,
    iline: int = None,
    xline: int = None,
    offset: int = None,
    istep: int = None,
    xstep: int = None,
    ostep: int = None,
    xloc: int = None,
    yloc: int = None,
    apply_scalar: bool = False,
) -> Dict:
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
    apply_scalar : bool
        apply scalar to 'i-interval' and 'x-interval'

    Returns
    -------
    Dict
        Dict of meta information 
    """


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
    count = segy.trace_count
    segy.close_file()
    return count


def scan_unsorted3D(
    fname,
    iline: int,
    xline: int,
):
    """
    Scan an unsored 3D SEG-Y file and get the geometry
    """


def load_unsorted3D(fname: str, geom: Dict = None):
    """
    using fromfile_without_scan to load unsorted 3D segy file
    """


def create(segy_out: str,
           binary_in,
           shape: Tuple = None,
           format: int = 5,
           dt: int = 2000,
           start_time: int = 0,
           iline_interval: float = 25,
           xline_interval: float = 25,
           min_iline: int = 1,
           min_xline: int = 1,
           custom_info: List = []) -> None:
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
