from __future__ import annotations
import cigsegy
import typing
from typing import Dict, List, Tuple, TypedDict, Union
import numpy

_Shape = typing.Tuple[int, ...]

__all__ = ["Pysegy", "fromfile", "tofile", "create_by_sharing_header"]

kBinaryHeaderHelp: Dict[int, Tuple[str, int]]  # binary header help

# trace header help
kTraceHeaderHelp: Dict[int, Tuple[str, int]]


class MetaInfo(TypedDict):
    """
    Meta information of a segy file
    """

    sizeX: int  # nt
    sizeY: int  # nx
    sizeZ: int  # ni
    trace_count: int  # trace count
    sample_interval: int  # dt
    data_format: int  # 1 (>4f-ibm) or 5 (>4f-ieee)
    Y_interval: float  # interval of crossline
    Z_interval: float  # interval of inline
    start_time: int  # start time for each trace
    scalar: int  # scalar applied to Y_interval and Z_interval / (real world) x and y

    min_inline: int  # min inline number
    max_inline: int  # max inline number
    min_crossline: int  # min crossline number
    max_crossline: int  # max crossline number

    isNormalSegy: bool  # is a normal segy file? normal means a rectangle region

    fillNoValue: float  # the fill value for missing traces

    # field information
    inline_field: int  # inline number location, default is 189
    crossline_field: int  # crossline number location, defualt is 193
    X_field: int  # x (real world) value location, default is 73
    Y_field: int  # y (real world) value location, default is 77

    inline_step: int  # inline step, default is 1, step=3 means the inline is like 100, 103, 106, ...
    crossline_step: int  # crossline step, default is 1


class Pysegy():

    @typing.overload
    def __init__(self, segy_name: str) -> None:
        """
        Reading mode

        Parameters
        ----------
        segy_name : str
            the input segy format file
        """

    @typing.overload
    def __init__(self, sizeX: int, sizeY: int, sizeZ: int) -> None:
        """
        Creating mode

        Parameters
        ----------
        sizeX : int
            the number of samples per trace,
        sizeY : int
            the number of crossline,
        sizeZ : int
            the number of inline.
        """

    @typing.overload
    def __init__(self, binary_name: str, sizeX: int, sizeY: int,
                 sizeZ: int) -> None:
        """
        Creating mode

        Parameters
        ----------
        binary_name : str
            the input binary file name,
        sizeX : int
            the number of samples per trace,
        sizeY : int
            the number of crossline,
        sizeZ : int
            the number of inline.
        """

    @typing.overload
    def create(self, segy_out_name: str, custom_info: List[str] = []) -> None:
        """
        create a segy

        Parameters
        ----------
        segy_out_name : str
            the output file name to create segy format
        custom_info : List[str]
            textual header info by user custom, max: 12 rows
                each row is less than 76 chars
        """

    @typing.overload
    def create(self,
               segy_out_name: str,
               src: numpy.ndarray[numpy.float32],
               custom_info: List[str] = []) -> None:
        """
        create a segy from src

        Parameters
        ----------
        segy_out_name : str
            the output file name to create segy format
        src : numpy.ndarray[numpy.float32]
            data in numpy array format to create
        custom_info: List[str]
            textual header info by user custom, max: 12 rows
                each row is less than 76 chars
        """

    def metaInfo(self) -> str:
        """ 
        the meta info for the whole segy file in string format

        Returns
        -------
        str
            meta information string
        """

    def collect(self,
                beg: int = -1,
                end: int = 0) -> numpy.ndarray[numpy.float32]:
        """
        collect traces as a 2D data from the `segy_in` file in
        range of [beg, end), beg < 0 means collect all traces,
        end < 0 means collect traces from beg to the last trace,
        end == 0 means read the beg-th trace (one trace).

        Parameters
        ----------
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

    @typing.overload
    def read(self) -> numpy.ndarray[numpy.float32]:
        """
        read whole volume from segyfile

        Returns
        -------
        numpy.ndarray[numpy.float32]
            data, 3D array
        """

    @typing.overload
    def read(self, startZ: int, endZ: int, startY: int, endY: int, startX: int,
             endX: int) -> numpy.ndarray[numpy.float32]:
        """ 
        read a volume with index, the volume size is 
        [startZ:endZ, startY:endY, startX:endX]

        Returns
        -------
        numpy.ndarray
            data, 3D array
        """

    @typing.overload
    def read(self, sizeY: int, sizeZ: int, minY: int, minZ: int) -> numpy.ndarray[numpy.float32]:
        """ 
        read without scanning

        Returns
        -------
        numpy.ndarray
            data, 3D array
        """

    def read_cross_slice(self, iY: int) -> numpy.ndarray[numpy.float32]:
        """
        read a crossline slice with index

        Returns
        -------
        numpy.ndarray
            data, 2D array
        """

    def read_inline_slice(self, iZ: int) -> numpy.ndarray[numpy.float32]:
        """
        read an inline slice with index

        Returns
        -------
        numpy.ndarray
            data, 2D array
        """

    def read_time_slice(self, iX: int) -> numpy.ndarray[numpy.float32]:
        """
        read a time slice with index

        Returns
        -------
        numpy.ndarray
            data, 2D array
        """

    def read_trace(self, iZ: int, iY: int) -> numpy.ndarray[numpy.float32]:
        """
        read a trace with index

        Returns
        -------
        numpy.ndarray
            data, 1D array
        """

    def cut(self,
            outname: str,
            startZ: int,
            endZ: int,
            startY: int,
            endY: int,
            startX: int,
            endX: int,
            custom_info: List[str] = []) -> None:
        """
        cut a sub volume data in segy format,
        range is [startZ:endZ, startY:endY, startX:endX],
        endZ/Y/X are not included.

        Parameters
        ----------
        outname : str
            out segy file name for the cutted data
        startZ : int
            start index in inline dimension
        endZ : int
            end index in inline dimension
        startY : int
            start index in crossline dimension
        endY : int
            end index in crossline dimension
        startX : int
            start index in time dimension
        endX : int
            end index in time dimension
        custom_info : List[str]
            textual header info by user custom, max: 12 rows
                each row is less than 76 chars
        """

    def cut(self,
            outname: str,
            startZ: int,
            endZ: int,
            startY: int,
            endY: int,
            custom_info: List[str] = []) -> None:
        """
        cut a sub volume data in segy format,
        range is [startZ:endZ, startY:endY, :],
        endZ/Y are not included.

        Parameters
        ----------
        outname : str
            out segy file name for the cutted data
        startZ : int
            start index in inline dimension
        endZ : int
            end index in inline dimension
        startY : int
            start index in crossline dimension
        endY : int
            end index in crossline dimension
        custom_info : List[str]
            textual header info by user custom, max: 12 rows
                each row is less than 76 chars
        """

    def cut(self,
            outname: str,
            startX: int,
            endX: int,
            custom_info: List[str] = []) -> None:
        """
        cut a sub volume data in segy format,
        range is [:, :, startX:endX],
        endX is not included.

        Parameters
        ----------
        outname : str
            out segy file name for the cutted data
        startX : int
            start index in time dimension
        endX : int
            end index in time dimension
        custom_info : List[str]
            textual header info by user custom, max: 12 rows
                each row is less than 76 chars
        """

    def scan(self) -> None:
        """
        scan the whole segy file
        """

    def setCrosslineLocation(self, xline: int) -> None:
        """ 
        set the crossline field of trace headers (for reading segy)
        recommend: 193, 17, 21, default is 193.
        """

    def setDataFormatCode(self, format: int) -> None:
        """
        set data format code (for creating segy).
        Only support: 
        1 for 4 bytes IBM float
        5 for 4 bytes IEEE float
        default is 5.
        """

    def setFillNoValue(self, fills: float) -> None:
        """ 
        set a value for filling the missing trace (for reading segy),
        can be any float number or np.nan
        """

    def setInlineLocation(self, iline: int) -> None:
        """ 
        set the crossline field of trace headers (for reading segy)
        recommend: 189, 5, 9, default is 189.
        """

    def setXLocation(self, xfield: int) -> None:
        """
        set the x field of trace headers (for reading segy)
        recommend: 73, 181
        """

    def setYLocation(self, yfield: int) -> None:
        """
        set the y field of trace headers (for reading segy)
        recommend: 77, 185
        """

    def setMinCrossline(self, minXline: int) -> None:
        """
        set min crossline number (for creating segy)
        """

    def setMinInline(self, minInline: int) -> None:
        """
        set min inline number (for creating segy)
        """

    def setSampleInterval(self, dt: int) -> None:
        """ 
        set sample interval for creating segy,
        the defualt is 4000, i.e. 4 ms
        """

    def setStartTime(self, start_time: int) -> None:
        """
        set start time for creating segy
        the defualt is 0
        """

    def setInlineInterval(self, dx: float) -> None:
        """ 
        set X interval for creating segy
        """

    def setCrosslineInterval(self, dy: float) -> None:
        """ 
        set Y interval for creating segy
        """

    def set_size(self, sizeX: int, sizeY: int, sizeZ: int) -> None:
        """
        set data shape for creating segy or reading by ignoring headers.

        Parameters
        ----------
        sizeX : int
            number of samples per trace
        sizeY : int
            number of crossline
        sizeZ : int
            number of inline
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

    def tofile(self, binary_out_name: str) -> None:
        """ 
        read a segy file and convert it to a binary file.

        Parameters
        ----------
        binary_out_name : str
            the output binary file name
        """

    def close_file(self) -> None:
        """
        close mmap for reading mode
        """

    def setInlineStep(self, step: int) -> None:
        """
        set inline step
        """

    def setCrosslineStep(self, step: int) -> None:
        """
        set inline step
        """

    def setSteps(self, istep: int, xstep: int) -> None:
        """
        set inline and crossline step
        """

    def get_lineInfo(self) -> numpy.ndarray:
        """
        get lineinfo: [sizeZ, 6]
        each raw: [inline, crossline_start, crossline_end, 
        trace_start, trace_end, count]

        Returns
        -------
        numpy.ndarray
            each row is [inline, crossline_start, crossline_end, 
                trace_start, trace_end, count]
        """

    @property
    def trace_count(self) -> int:
        """
        get the total traces
        """

    @property
    def nt(self) -> int:
        """
        get dimension nt
        """

    @property
    def is_crossline_fast_order(self) -> bool:
        """
        The fast order is crossline order or not
        """

    def get_binary_header(self, raw=False) -> numpy.ndarray[numpy.uint8]:
        """
        return binary header chars, 400 chars, 

        Parameters
        ----------
        raw : bool
            if raw is True, return the raw stream (big endian),
            if raw is False, each element of the stream is swap to little endian

        Returns
        -------
        numpy.ndarray[numpy.uint8]
            400 elements
        """

    def get_trace_header(self,
                         n: int,
                         raw: bool = False) -> numpy.ndarray[numpy.uint8]:
        """
        return the n-th trace header, 240 chars

        Parameters
        ----------
        n : int
            the n-th trace
        raw : bool
            if raw is True, return the raw stream (big endian),
            if raw is False, each element of the stream is swap to little endian

        Returns
        -------
        numpy.ndarray[numpy.uint8]
            240 elements
        """

    def get_trace(self,
                  n: int,
                  raw: bool = False) -> numpy.ndarray[numpy.uint8]:
        """
        return the n-th trace, (240 + sizeX * 4) chars

        Parameters
        ----------
        n : int
            the n-th trace
        raw : bool
            if raw is True, return the raw stream (big endian),
            if raw is False, each element of the stream is swap to little endian

        Returns
        -------
        numpy.ndarray[numpy.uint8]
            240 + sizeX * 4 elements
        """

    def get_trace_keys(self, keys: List, length: List, beg: int,
                       end: int) -> numpy.ndarray[numpy.int32]:
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

    def get_metaInfo() -> MetaInfo:
        """
        get metainfo in class `MetaInfo` format
        """

    pass


def disable_progressbar() -> None:
    """
    disable progress bar
    """


def fromfile(segy_name: str,
             iline: int = 189,
             xline: int = 193,
             istep: int = 1,
             xstep: int = 1) -> numpy.ndarray[numpy.float32]:
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


def fromfile_without_scan(segy_name: str,
                          ni: int,
                          nx: int,
                          il_min: int,
                          xl_min: int,
                          iline: int = 189,
                          xline: int = 193,
                          istep: int = 1,
                          xstep: int = 1) -> numpy.ndarray[numpy.float32]:
    """
    reading from a segy file without scanning.

    Parameters
    ----------
    segy_name : str
        the input segy file name
    ni : int
        number of inline
    nx : int
        number of crossline
    il_min : int
        the min number of inline
    xl_min : int
        the min number of crossline
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


def tofile(segy_name: str,
           out_name: str,
           iline: int = 189,
           xline: int = 193,
           istep: int = 1,
           xstep: int = 1) -> None:
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
    """


@typing.overload
def create_by_sharing_header(segy_name: str,
                             header_segy: str,
                             src: numpy.ndarray[numpy.float32],
                             iline: int = 189,
                             xline: int = 193,
                             istep: int = 1,
                             xstep: int = 1,
                             offset: Union[list, dict] = None,
                             custom_info: List[str] = []) -> None:
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
    iline : int
        the inline number field of header segy
    xline : int
        the crossline number field of header segy
    istep : int
        the step of inline numbers
    xstep : int
        the step of crossline numbers
    offset : int
        the offset of a sub data from a original data, e.g., 
            dsub = d[256:400, 500: 100:], offset = [256, 500, 100]
    custom_info : List[str]
        textual header info by user custom, max: 12 rows
            each row is less than 76 chars, use it when offset is not None
    """


@typing.overload
def create_by_sharing_header(segy_name: str,
                             header_segy: str,
                             src_file: str,
                             shape: Union[tuple, list],
                             iline: int = 189,
                             xline: int = 193,
                             istep: int = 1,
                             xstep: int = 1,
                             offset: Union[list, dict] = None,
                             custom_info: List[str] = []) -> None:
    """
    create a segy and its header is from an existed segy.

    Parameters
    ----------
    segy_name : str
        the out segy name
    header_segy : str
        the header segy file
    src_file : str
        source data file
    shape : tuple or list
        like: (sizeZ, sizeY, sizeX) or (iline, xline, nt)
    iline : int
        the inline number field of header segy
    xline : int
        the crossline number field of header segy
    istep : int
        the step of inline numbers
    xstep : int
        the step of crossline numbers
    offset : tuple or list
        the offset of a sub data from a original data, e.g., 
            dsub = d[256:400, 500: 100:], offset = [256, 500, 100]
    custom_info : List[str]
        textual header info by user custom, max: 12 rows
            each row is less than 76 chars, use it when offset is not None
    """


def _load_prestack3D(segy_name: str,
                     shape: Union[List, Tuple],
                     min_iline: int,
                     max_iline: int,
                     min_xline: int,
                     max_xline: int,
                     min_offset: int,
                     max_offset: int,
                     istep: int,
                     xstep: int,
                     ostep: int,
                     iline: int,
                     xline: int,
                     offset: int = 37,
                     fill: float = 0):
    """
    Load 3D prestack data (4D array)
    """


def modify_bin_key(segy_name: str,
                   loc: int,
                   value: Union[float, int],
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
                     value: Union[float, int],
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
