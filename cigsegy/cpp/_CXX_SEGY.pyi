# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
"""
_CXX_SEGY: A C++ core library for SEG-Y file operations, exposed to Python using pybind11.

This module provides high-performance, low-level functionality for handling SEG-Y files, a common format for geophysical data storage. The library is written in C++ for efficiency and is wrapped with pybind11 to be accessible from Python.

Key Features:
--------------
- **Scan SEG-Y Files**: Quickly scan SEG-Y files to extract metadata and other essential information without loading the entire file into memory.
  
- **Read SEG-Y Files**: Efficiently read seismic data and headers from SEG-Y files, allowing for fast data access.
  
- **Write SEG-Y Files**: Modify SEG-Y files with full control over file headers and format specifications (including header and data).
  
- **Create SEG-Y Files**: Initialize new SEG-Y files, specifying headers, formats, and data structures according to industry standards.

- **Efficient Functions**: data converting between IBM and IEEE floating-point formats.
"""

from typing import List, Tuple, overload

import numpy as np


class Pysegy:
    """
    Pysegy: A Python class for scanning, accessing, reading, and modifying SEG-Y file headers and data.

    The `Pysegy` class provides a high-level interface for working with SEG-Y files, which are widely used in the geophysical and seismic industries. 
    This class is implemented in C++ for performance reasons and is exposed to Python using pybind11, making it both powerful and easy to use.
    
    Example:
    --------
    ```python
    from cigsegy import Pysegy
    # from cigsegy.cpp._CXX_SEGY import Pysegy

    segy_file = Pysegy("example.segy")
    segy_file.setLocations(189, 193)
    segy_file.scan()
    # read the entire SEG-Y file
    d = segy_file.read()
    d.close()
    ```
    """

    def __init__(self, segy_name: str, write: bool = False) -> None:
        ...

    @property
    def ntrace():
        """
        number of traces
        """

    @property
    def nt():
        """
        nt
        """

    @property
    def ndim():
        """
        ndim
        """

    @property
    def shape():
        """
        ndim
        """

    @property
    def is_scanned():
        """
        is_scanned
        """

    def close() -> None:
        """
        close file
        """

    def show_progress(self, show: bool) -> None:
        """
        show progress
        """

    def setLocations(self, iline: int, xline: int, offset: int = 37) -> None:
        """ 
        set the inline, crossline and offset field of trace headers
        """

    def setInlineLocation(self, iline: int) -> None:
        """ 
        set the crossline field of trace headers (for reading segy)
        recommend: 189, 5, 9, default is 189.
        """

    def setCrosslineLocation(self, xline: int) -> None:
        """ 
        set the crossline field of trace headers (for reading segy)
        recommend: 193, 17, 21, default is 193.
        """

    def setOffsetLocation(self, offset: int) -> None:
        """ 
        set the offset field of trace headers (for reading segy)
        recommend: 37
        """

    def setXLocation(self, xloc: int) -> None:
        """
        set the x field of trace headers (for reading segy)
        recommend: 181, 73
        """

    def setYLocation(self, yloc: int) -> None:
        """
        set the y field of trace headers (for reading segy)
        recommend: 185, 77
        """

    def setXYLocations(self, xloc: int, yloc: int) -> None:
        """
        set the x and y field of trace headers (for reading segy)
        """

    def setInlineStep(self, istep: int) -> None:
        """ 
        set the step of inline (for reading segy)
        """

    def setCrosslineStep(self, xstep: int) -> None:
        """ 
        set the step of crossline (for reading segy)
        """

    def setOffsetStep(self, ostep: int) -> None:
        """ 
        set the step of offset (for reading segy)
        """

    def setSteps(self, istep: int, xstep: int, ostep: int = 1) -> None:
        """
        set inline, crossline and offset setps
        """

    def setFill(self, fill: float) -> None:
        """ 
        set a value for filling the missing trace (for reading segy),
        can be any float number or np.nan
        """

    def textual_header(self, coding='u') -> str:
        """
        obtain the 3200 bytes textual header.

        Parameters
        ---------
        coding : {'u', 'e', 'a'}, optional
            force the 3200 bytes in coding format, 'u' means guess by default,
            'e' means EBCDIC, 'a' means ASCII

        Returns
        -------
        str
            3200 bytes textual header
        """

    def bkeyi2(self, loc: int) -> int:
        """
        return the key value of the binary header in loc, 2 bytes
        """

    def bkeyi4(self, loc: int) -> int:
        """
        return the key value of the binary header in loc, 4 bytes
        """

    def keyi2(self, n: int, loc: int) -> int:
        """
        return the key value of the n-th trace header in loc, 2 bytes
        """

    def keyi4(self, n: int, loc: int) -> int:
        """
        return the key value of the n-th trace header in loc, 4 bytes
        """

    def iline(self, n: int) -> int:
        """
        return the inline value of the n-th trace
        """

    def xline(self, n: int) -> int:
        """
        return the crossline value of the n-th trace
        """

    def offset(self, n: int) -> int:
        """
        return the offset value of the n-th trace
        """

    def coordx(self, n: int) -> int:
        """
        return the X coordinate value of the n-th trace
        """

    def coordy(self, n: int) -> int:
        """
        return the Y coordinate value of the n-th trace
        """

    def get_binary_header(self) -> np.ndarray:
        """
        Obtain the 400 bytes binary header, return np.ndarray[np.uint8].

        Note: the binary header is raw data, i.e., big endian
        """

    def get_trace_header(self, n) -> np.ndarray:
        """
        Obtain the 240 bytes trace header of the n-th trace, return np.ndarray[np.uint8].

        Note: the trace header is raw data, i.e., big endian
        """

    def get_trace_keys(
        self,
        keys: List,
        length: List,
        beg: int,
        end: int,
    ) -> np.ndarray:
        """
        get the trace keys from beg to end

        Parameters
        ------------
        keys : List
            location
        length : List
            keys' length
        beg : int
            the start trace index
        end : int 
            the stop trace index (not included)

        Returns
        -------
        numpy.ndarray[numpy.int32]
            shape is (end-beg, len(keys))
        """

    def itrace(self, n: int) -> np.ndarray:
        """
        exract the data in the n-th trace
        """

    @overload
    def collect(self, beg: int, end: int, tbeg: int, tend: int) -> np.ndarray:
        """
        collect the data from the beg-th trace to the end-th trace (not included),
        time window is from tbeg to tend

        Parameters
        ----------
        beg : int
            the start trace index
        end : int
            the stop trace index (not included)
        tbeg : int
            the start time index
        tend : int
            the stop time index (not included)
        
        Returns
        -------
        numpy.ndarray
            shape is (end-beg, tend-tbeg)
        """

    @overload
    def collect(self, index: np.ndarray, tbeg: int, tend: int) -> np.ndarray:
        """
        Collect the data from the index array. This funcion is useful 
        when collect unsorted traces. if a index is negative, the trace 
        will be filled by fills (0 is default).

        Parameters
        -----------
        index : numpy.ndarray
            1D array, np.int32, 
        tbeg : int
            the start time index
        tend : int
            the stop time index (not included)
        
        Returns
        -------
        numpy.ndarray
            shape is (index.shape[0], tend-tbeg)

        Examples
        ---------
        >>> index = np.array([10, -1, 200, 300])
        >>> data = Pysegy('out.segy').collect(index, 0, 1000)
        >>> data.shape # (4, 1000)
        >>> np.allclose(data[1], 0) # True
        """

    def get_keylocs(self) -> dict:
        """
        get key locations in dict format
        """

    def get_metainfo(self) -> dict:
        """
        get metainfo in dict format.
        """

    def get_lineInfo(self) -> np.ndarray:
        """
        get lineinfo in np.ndarray format.
        if ndim == 3, each raw: 
            [iline, xstart, xend, trace_start, trace_end]
        if ndim == 4, each raw:
            [iline, xline, ostart, oend, trace_start, trace_end]


        Returns
        -------
        numpy.ndarray
            if ndim == 3, each raw: 
                [iline, xstart, xend, trace_start, trace_end]
            if ndim == 4, each raw:
                [iline, xline, ostart, oend, trace_start, trace_end]
        """

    def set_segy_type(self, ndim) -> None:
        """
        set the segy type, 2 or 3 or 4
        2 for a collect of traces,
        3 for poststack 3D volume,
        4 for prestack 4D volume/gather
        """

    def scan(self) -> None:
        """
        fast scan the file to obtain some meta infos
        """

    def read4d(
        self,
        ib: int,
        ie: int,
        xb: int,
        xe: int,
        ob: int,
        oe: int,
        tb: int,
        te: int,
    ) -> np.ndarray:
        """
        """

    def read3d(
        self,
        ib: int,
        ie: int,
        xb: int,
        xe: int,
        tb: int,
        te: int,
    ) -> np.ndarray:
        """
        """

    def read(self) -> np.ndarray:
        """
        """

    def read_tslice(self,
                    t: int,
                    stepi: int = 1,
                    stepx: int = 1) -> np.ndarray:
        """
        read time slice with stepi and stepx, which may faster than read()
        return a 2D array with shape: ((ni+stepi-1)//stepi, (nx+stepx-1)//stepx)
        """

    def tofile(
        self,
        binary_out_name: str,
        as_2d=False,
    ) -> None:
        """
        """

    def cut(self,
            outname: str,
            ranges: List,
            is2d: bool = False,
            textual: str = "") -> None:
        """
        """

    @overload
    def create_by_sharing_header(
        self,
        segy_name: str,
        src: np.ndarray,
        start: List,
        is2d: bool = False,
        textual: str = "",
    ) -> None:
        """
        """

    @overload
    def create_by_sharing_header(
        self,
        segy_name: str,
        src_file: str,
        shape: List,
        start: List,
        is2d: bool = False,
        textual: str = "",
    ) -> None:
        """
        """

    def set_bkeyi2(self, loc: int, val: int) -> None:
        """
        """

    def set_bkeyi4(self, loc: int, val: int) -> None:
        """
        """

    def set_keyi2(self, n: int, loc: int, val: int) -> None:
        """
        """

    def set_keyi4(self, n: int, loc: int, val: int) -> None:
        """
        """

    def set_iline(self, n: int, val: int) -> None:
        """
        """

    def set_xline(self, n: int, val: int) -> None:
        """
        """

    def set_offset(self, n: int, val: int) -> None:
        """
        """

    def set_coordx(self, n: int, val: int) -> None:
        """
        """

    def set_coordy(self, n: int, val: int) -> None:
        """
        """

    def write_itrace(self, data: np.ndarray, n: int) -> None:
        """
        """

    @overload
    def write_traces(
        self,
        data: np.ndarray,
        beg: int,
        end: int,
        tbeg: int,
        tend: int,
    ) -> None:
        """
        """

    @overload
    def write_traces(
        self,
        data: np.ndarray,
        index: np.ndarray,
        tbeg: int,
        tend: int,
    ) -> None:
        """
        """

    def write(self, data: np.ndarray) -> None:
        """
        """

    def write3d(
        self,
        data: np.ndarray,
        ib: int,
        ie: int,
        xb: int,
        xe: int,
        tb: int,
        te: int,
    ) -> None:
        """
        """

    def write4d(
        self,
        data: np.ndarray,
        ib: int,
        ie: int,
        xb: int,
        xe: int,
        ob: int,
        oe: int,
        tb: int,
        te: int,
    ) -> None:
        """
        """


def create_segy(
    segyname: str,
    src: np.ndarray,
    keys: np.ndarray,
    textual: str,
    bheader: np.ndarray,
    theader: np,
    ndarray,
) -> None:
    """
    create a new segy file from src data

    Parameters
    ----------
    segyname : str
        the name of the segy file
    src : np.ndarray, np.float32
        the source data in np.float32
    keys : np.ndarray, np.int32
        the keys of the traces, shape is (N, keysize), where N is the number of traces,
        and keysize can be 4, 5. N = src.shape[0]*...*src.shape[-2].
        if kesize is 4, each row contains [iline, xline, X, Y]
        if keysize is 5, each row contains [iline, xline, offset, X, Y]
    textual : str
        the textual header, 3200 bytes
    bheader : np.ndarray, np.uint8
        the binary header, 400 bytes
    theader : np.ndarray, np.uint8
        the trace header template, 240 bytes
    """


def ibm_to_ieee(value: float, is_big_endian: bool) -> float:
    """
    convert IBM floating point to IEEE floating point
    """


def ieee_to_ibm(value: float,
                is_litte_endian_input: bool,
                is_big_endian_output: bool = True) -> float:
    """
    convert IEEE floating point to IBM floating point
    """


def ibms_to_ieees(ibms: np.ndarray, is_big_endian: bool) -> np.ndarray:
    """
    convert IBM floating array to IEEE floating points
    """


def ieees_to_ibms(ieees: np.ndarray,
                  is_litte_endian_input: bool) -> np.ndarray:
    """
    convert IEEE floating array to IBM floating points
    """


def set_progress_callback(func: callable[[int, int], None]) -> None:
    """
    set the progress callback function

    Parameters
    ----------
    func : callable
        A callback function that takes two integer arguments: current and total.
    """


def set_global_show_progress(show: bool) -> None:
    """
    set the global progress bar visiable or not
    """
