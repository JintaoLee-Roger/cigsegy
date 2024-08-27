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

    def __init__(self, segy_name: str) -> None:
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

    def close() -> None:
        """
        close file
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
        """

    def bkeyi4(self, loc: int) -> int:
        """
        """

    def keyi2(self, n: int, loc: int) -> int:
        """
        """

    def keyi4(self, n: int, loc: int) -> int:
        """
        """

    def iline(self, n: int) -> int:
        """
        """

    def xline(self, n: int) -> int:
        """
        """

    def offset(self, n: int) -> int:
        """
        """

    def coordx(self, n: int) -> int:
        """
        """

    def coordy(self, n: int) -> int:
        """
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
            the stop trace index

        Returns
        -------
        numpy.ndarray[numpy.int32]
            shape is (end-beg, len(keys))
        """

    def itrace(self, n: int) -> np.ndarray:
        """
        """

    @overload
    def collect(self, index: np.ndarray, tbeg: int, tend: int) -> np.ndarray:
        """
        """

    @overload
    def collect(self, beg: int, end: int, tbeg: int, tend: int) -> np.ndarray:
        """
        """

    def get_keylocs(self) -> dict:
        """
        """

    def get_metainfo(self) -> dict:
        """
        """

    def set_segy_type(self, ndim) -> None:
        """
        """

    def scan(self) -> None:
        """
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
        shape: List,
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


def ibm_to_ieee(value: float, is_big_endian: bool) -> float:
    """
    convert IBM floating point to IEEE floating point
    """


def ieee_to_ibm(value: float, is_little_endian: bool) -> float:
    """
    convert IEEE floating point to IBM floating point
    """


def ibms_to_ieees(ibms: np.ndarray, is_big_endian: bool) -> np.ndarray:
    """
    convert IBM floating array to IEEE floating points
    """


def ieees_to_ibms(ieees: np.ndarray, is_little_endian: bool) -> np.ndarray:
    """
    convert IEEE floating array to IBM floating points
    """
