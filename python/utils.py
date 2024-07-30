# Copyright (c) 2023 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
#
# github: https://github.com/JintaoLee-Roger

from pathlib import Path
import sys
import struct
import numpy as np
from typing import List, Tuple, Dict, Union
from .cigsegy import (Pysegy, disable_progressbar, MetaInfo, kTraceHeaderHelp)  # type: ignore


def get_trace_keys(segy: Union[str, Pysegy],
                   keyloc: Union[int, List, np.ndarray],
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

    if isinstance(keyloc, int):
        keyloc = [keyloc]
    elif isinstance(keyloc, np.ndarray):
        assert keyloc.ndm == 1
    assert isinstance(keyloc, List) or isinstance(keyloc, np.ndarray)

    if isinstance(segy, (str, Path)):
        segyc = Pysegy(str(segy))
    else:
        segyc = segy

    if beg < 0:
        beg = 0
        end = segyc.trace_count
    if end < 0:
        end = segyc.trace_count
    if end == 0:
        end = beg + 1

    assert beg < end, f"beg ({beg} > end ({end}))"
    assert beg >= 0, f"beg ({beg}) >= 0"
    assert end <= segyc.trace_count, f"end <= trace_count"
    if force is None:
        for key in keyloc:
            assert key in kTraceHeaderHelp.keys(
            ), f"keyloc ({key}) in tracehelper"

    if force is None:
        length = [kTraceHeaderHelp[key][1] for key in keyloc]
        raw = False
    else:
        assert len(keyloc) == 1
        length = [force]
        raw = True

    range_func = range
    if end - beg > 1000000:
        try:
            from tqdm import trange
            range_func = trange
        except:
            pass

    unpackstr = {2: 'h', 4: 'i', 5: 'f', 8: 'q', 9: 'd'}

    out = []

    for i in range_func(beg, end):  # TODO: 需要优化, end-beg特别大时，打印费时间
        th = segyc.get_trace_header(i, raw=raw)

        for key, l in zip(keyloc, length):
            if l > 9 or l == 6:
                out.append(0)
            elif l == 1:
                out.append(th[key - 1])
            else:
                upstr = unpackstr[l]
                if l % 2 != 0:
                    l -= 1
                sub = th[key - 1:key + l - 1]

                if force is not None:  # swap endian, th is raw header
                    sub = sub.copy()[::-1]

                v = struct.unpack(upstr, struct.pack(f'{l}B', *sub))[0]
                out.append(v)

    if isinstance(segy, (str, Path)):
        segyc.close_file()

    out = np.array(out)
    if len(keyloc) > 1:  # multi-locations
        out = out.reshape(-1, len(keyloc))

    if end - beg == 1 and len(keyloc) == 1:  # one value
        return out[0]

    return out


def get_trace_keys2(segy: Union[str, Pysegy],
                    keyloc: Union[int, List, np.ndarray],
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

    if isinstance(keyloc, int):
        keyloc = [keyloc]
    elif isinstance(keyloc, np.ndarray):
        assert keyloc.ndm == 1
    assert isinstance(keyloc, List) or isinstance(keyloc, np.ndarray)

    if isinstance(segy, (str, Path)):
        segyc = Pysegy(str(segy))
    else:
        segyc = segy

    if beg < 0:
        beg = 0
        end = segyc.trace_count
    if end < 0:
        end = segyc.trace_count
    if end == 0:
        end = beg + 1

    assert beg < end, f"beg ({beg} > end ({end}))"
    assert beg >= 0, f"beg ({beg}) >= 0"
    assert end <= segyc.trace_count, f"end <= trace_count"
    if force is None:
        for key in keyloc:
            assert key in kTraceHeaderHelp.keys(
            ), f"keyloc ({key}) in tracehelper"

    if force is None:
        length = [kTraceHeaderHelp[key][1] for key in keyloc]
        raw = False
    else:
        assert len(keyloc) == 1
        length = [force]
        raw = True

    print(keyloc, length, beg, end)
    out = segyc.get_trace_keys(keyloc, length, beg, end)

    if isinstance(segy, (str, Path)):
        segyc.close_file()

    if len(keyloc) == 1:
        out = out.flatten()

    if end - beg == 1 and len(keyloc) == 1:  # one value
        return out[0]

    return out


def eval_xline(segy: Pysegy) -> List:
    """
    To guess the crossline location of the segy

    Parameters
    ----------
    segy : Pysegy
        input Pysegy class

    Returns
    -------
    List
        possible result
    """
    tracecout = segy.trace_count
    options = [193, 17, 21, 13]
    select = [193, 17, 21, 13]
    for op in options:
        if get_trace_keys(segy, op, tracecout - 1) - get_trace_keys(
                segy, op, 0) > tracecout // 2:
            select.remove(op)
            continue
        d = get_trace_keys(segy, op, tracecout // 2, tracecout // 2 + 10)
        if sum(d >= 0) != 10:
            select.remove(op)
            continue
        if len(np.where(np.diff(d) == 0)[0]) > 1:
            select.remove(op)
            continue

    if select == []:
        raise RuntimeError("Cannot evaluate crossline location")

    return select


def eval_iline(segy: Pysegy) -> List:
    """
    To guess the inline location of the segy

    Parameters
    ----------
    segy : Pysegy
        input Pysegy class

    Returns
    -------
    List
        possible result
    """
    tracecout = segy.trace_count
    options = [189, 5, 9, 221]
    select = [189, 5, 9, 221]
    force = 4

    for op in options:
        l0 = get_trace_keys(segy, op, 0, force=force)
        ll = get_trace_keys(segy, op, tracecout - 1, force=force)
        l2 = get_trace_keys(segy, op, tracecout // 2, force=force)
        if sum([x >= 0 for x in [l0, ll, l2]]) != 3:
            select.remove(op)
            continue
        if l0 == ll or l0 == l2 or ll == l2:
            select.remove(op)
            continue
        if max([l0, ll, l2]) - min([l0, ll, l2]) > min(tracecout // 10 - 1,
                                                       100000):
            select.remove(op)
            continue
        if (l0 < l2 and l2 > ll) or (l0 > l2 and l2 < ll):
            select.remove(op)
            continue

    if select == []:
        raise RuntimeError("Cannot evaluate inline location")

    return select


def eval_xstep(segy: Pysegy, xline: int) -> int:
    """
    To guess the crossline step

    Parameters
    ----------
    segy : Pysegy
        input Pysegy class
    xline: int
        crossline location

    Returns
    -------
    List
        possible result
    """
    tracecout = segy.trace_count

    d = get_trace_keys(segy, xline, tracecout // 2, tracecout // 2 + 10)
    diff = np.diff(d)
    if diff.min() > 0:
        return diff.min()
    idx = np.where(diff <= 0)[0]
    if len(idx) > 1:
        return 0

    d = get_trace_keys(segy, xline, tracecout // 2 + idx[0] + 1,
                       tracecout // 2 + 10 + idx[0] + 1)
    diff = np.diff(d)
    if diff.min() > 0:
        return diff.min()
    else:
        return 0


def eval_istep(segy: Pysegy, iline: int) -> int:
    """
    To guess the inline step

    Parameters
    ----------
    segy : Pysegy
        input Pysegy class
    iline: int
        inline location

    Returns
    -------
    List
        possible result
    """
    tracecount = segy.trace_count
    i0 = get_trace_keys(segy, iline, tracecount // 2, force=4)
    x1 = tracecount // 2 + 1
    while get_trace_keys(segy, iline, x1, force=4) == i0:
        x1 += 1
    i1 = get_trace_keys(segy, iline, x1, force=4)

    x2 = x1 + 1
    while get_trace_keys(segy, iline, x2, force=4) == i1:
        x2 += 1
    i2 = get_trace_keys(segy, iline, x2, force=4)

    if i2 - i1 != i1 - i0:
        return 0
    else:
        return i2 - i1


def guess(segy_name: Union[str, Pysegy],
          iline=None,
          xline=None,
          istep=None,
          xstep=None,
          xloc=None,
          yloc=None) -> List:
    """
    TODO: scan 3 times when call metaInfo() ?
    guess the locations and steps of inline and crossline

    Parameters
    ----------
    segy_name : str or Pysegy 
        the input segy file

    Returns
    -------
    List
        locations, [loc1, loc2, ...], all possible loctaions,
          each location is like: [iline, xline, istep, xstep]
    """
    if isinstance(segy_name, (str, Path)):
        segy = Pysegy(str(segy_name))
    elif isinstance(segy_name, Pysegy):
        segy = segy_name
    else:
        raise TypeError("Invalid type of `segy_name`")
    
    xlines = eval_xline(segy) if xline is None else [xline]
    ilines = eval_iline(segy) if iline is None else [iline]
    xselect = []
    iselect = []
    xsteps = []
    isteps = []

    if xstep is not None:
        xsteps = [xstep]
        xselect = xlines
    else:
        for xline in xlines:
            xstep = eval_xstep(segy, xline)
            if xstep:
                xselect.append(xline)
                xsteps.append(xstep)
    
    if istep is not None:
        isteps = [istep]
        iselect = ilines
    else:
        for iline in ilines:
            istep = eval_istep(segy, iline)
            if istep:
                iselect.append(iline)
                isteps.append(istep)
    segy.close_file()

    out = []
    for i in range(len(iselect)):
        for x in range(len(xselect)):
            try:
                s = Pysegy(str(segy_name))
                s.setInlineLocation(iselect[i])
                s.setCrosslineLocation(xselect[x])
                s.setSteps(isteps[i], xsteps[x])
                t = s.metaInfo()
            except:
                continue

            out.append([iselect[i], xselect[x], isteps[i], xsteps[x]])

    if out == []:
        raise RuntimeError("cannot define the geometry through the location and steps")

    if xloc is None or yloc is None:
        s = Pysegy(str(segy_name))
        s.setInlineLocation(out[0][0])
        s.setCrosslineLocation(out[0][1])
        s.setSteps(out[0][2], out[0][3])
        s.setXLocation(181)
        s.setYLocation(185)
        s.scan()
        mt = s.get_metaInfo()
        if np.isnan(mt.Z_interval) or np.isnan(mt.Y_interval) or mt.Z_interval == 0 or mt.Y_interval == 0:
            xloc, yloc = 73, 77
        else:
            xloc, yloc = 181, 185
        s.close_file()

    for o in out:
        o += [xloc, yloc]

    return out


def convert_header(arr: np.ndarray,
                   help: Dict[int, Tuple[str, int]],
                   to_str: bool = False,
                   keyoffset: int = 0) -> Dict:
    """
    convert header (in np.ndarray[np.uint8]) to a dict

    Parameters
    ----------
    arr : np.ndarray
        little endian np.uint8 array
    help : Dict
        kBinaryHeaderHelp or kTraceHeaderHelp
    to_str : bool
        convert to str
    keyoffset : int
        if to_str is True, keyoffset is added to each key location,
            i.e., keyoffset = 3200, 1-4 -> 3201-3204

    Returns
    -------
    Dict
        dict as {loc: {helpstr, v}}
    """
    unpackstr = {2: 'h', 4: 'i', 5: 'f', 8: 'q', 9: 'd'}
    result = {}

    for key, (desc, length) in help.items():
        if length == 1:
            result[key] = (desc, arr[key - 1])
        elif length > 9 or length == 6:
            result[key] = (desc, 0)
        else:
            upstr = unpackstr[length]
            if length % 2 == 1:
                length -= 1
            subarr = arr[key - 1:key + length - 1]

            value = struct.unpack(upstr, struct.pack(f'{length}B', *subarr))[0]
            result[key] = (desc, value)

    if to_str:
        out = []
        for key, (desc, value) in result.items():
            length = help[key][1]
            if length > 1 and length % 2 == 1:
                length -= 1
            key += keyoffset
            out.append(
                f'Bytes {key:<4} - {key+length-1:<4}: {value:<8} -- {desc}')
        return out

    return result


def convert_trace(arr: np.ndarray) -> np.ndarray:
    """
    convert np.uint8 byte arr[4*n] to np.float32 array[n],
    n means the number of float points

    Parameters
    ----------
    arr : np.ndarray[np.uint8]
        input trace array in np.uint8 format

    Returns
    -------
    np.ndarray[np.float32]
        output trace array in np.float32 format
    """
    assert arr.shape[0] % 4 == 0
    n = arr.shape[0] // 4
    out = np.zeros((n), np.float32)
    for i in range(n):
        out[i] = struct.unpack('f', struct.pack('4B',
                                                *arr[i * 4:i * 4 + 4]))[0]

    return out


def progress_bar() -> None:
    """
    Disable progress bar in jupyter
    """
    injupyter = 'ipykernel_launcher.py' in sys.argv[0] or 'lab' in sys.argv[0]
    if injupyter:
        disable_progressbar()


def metainfo_to_dict(metainfo: MetaInfo, apply_scalar: bool = False) -> Dict:
    """
    convert MetaInfo to a Dict

    Parameters
    ----------
    metainfo : MetaInfo
        input Metainfo class
    apply_scalar : bool
        apply scalar to inline/crossline interval

    Returns
    -------
    Dict
        meta information in dict format
    """
    out = {}
    out['nt'] = metainfo.sizeX
    out['nx'] = metainfo.sizeY
    out['ni'] = metainfo.sizeZ
    out['trace_count'] = metainfo.trace_count
    out['dt'] = metainfo.sample_interval
    if metainfo.data_format == 1:
        out['dtype'] = '>4f-ibm'
    elif metainfo.data_format == 2:
        out['dtype'] = '>4i'
    elif metainfo.data_format == 3:
        out['dtype'] = '>2i'
    elif metainfo.data_format == 5:
        out['dtype'] = '>4f-ieee'
    elif metainfo.data_format == 8:
        out['dtype'] = '>1i'

    else:
        raise TypeError("don't support this data format")

    out['scalar'] = metainfo.scalar
    zinterval = metainfo.Z_interval
    yinterval = metainfo.Y_interval
    if apply_scalar:
        if metainfo.scalar == 0:
            metainfo.scalar = 1
        scalar = -1 / metainfo.scalar if metainfo.scalar < 0 else metainfo.scalar
        zinterval *= scalar
        yinterval *= scalar
    out['i-interval'] = zinterval
    out['x-interval'] = yinterval

    out['start_time'] = metainfo.start_time

    out['min-iline'] = metainfo.min_inline
    out['max-iline'] = metainfo.max_inline
    out['min-xline'] = metainfo.min_crossline
    out['max-xline'] = metainfo.max_crossline
    out['isnormal'] = metainfo.isNormalSegy

    out['iline'] = metainfo.inline_field
    out['xline'] = metainfo.crossline_field
    out['xloc'] = metainfo.X_field
    out['yloc'] = metainfo.Y_field
    out['istep'] = metainfo.inline_step
    out['xstep'] = metainfo.crossline_step

    out['fills'] = metainfo.fillNoValue

    return out
