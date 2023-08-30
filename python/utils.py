# Copyright (c) 2023 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
#
# github: https://github.com/JintaoLee-Roger

import sys
import struct
import numpy as np
from typing import List, Tuple, Dict, Union
from .cigsegy import (Pysegy, disable_progressbar, MetaInfo, kTraceHeaderHelp)


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

    if isinstance(segy, str):
        segyc = Pysegy(segy)
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

    for i in range_func(beg, end):
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

    if isinstance(segy, str):
        segyc.close_file()

    out = np.array(out)
    if len(keyloc) > 1:  # multi-locations
        out = out.reshape(-1, len(keyloc))

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
        l = []
        l.append(get_trace_keys(segy, op, 0))
        l.append(get_trace_keys(segy, op, 1))
        l.append(get_trace_keys(segy, op, 2))
        l.append(get_trace_keys(segy, op, tracecout // 2))
        l.append(get_trace_keys(segy, op, tracecout - 1))
        if sum([x > 0 for x in l]) != 5:
            select.remove(op)
            continue
        if l[0] == l[1] or l[0] == l[2] or l[1] == l[2]:
            select.remove(op)
            continue
        step = l[1] - l[0]
        if max(l) - min(l) > min((tracecout * step / 10), 10000 * step):
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
        if sum([x > 0 for x in [l0, ll, l2]]) != 3:
            select.remove(op)
            continue
        if l0 == ll or l0 == l2 or ll == l2:
            select.remove(op)
            continue
        if max([l0, ll, l2]) - min([l0, ll, l2]) > tracecout - 1:
            select.remove(op)
            continue
        l1 = get_trace_keys(segy, op, 1, force=force)
        ll2 = get_trace_keys(segy, op, tracecout - 2, force=force)
        if l0 != l1 or ll != ll2:
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
    l0 = get_trace_keys(segy, xline, 0)
    l1 = get_trace_keys(segy, xline, 1)
    l2 = get_trace_keys(segy, xline, 2)

    if l2 - l1 != l1 - l0:
        return 0
    else:
        return l2 - l1


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
    i0 = get_trace_keys(segy, iline, 0)
    x1 = 1
    while get_trace_keys(segy, iline, x1) == i0:
        x1 += 1
    i1 = get_trace_keys(segy, iline, x1)

    x2 = x1 + 1
    while get_trace_keys(segy, iline, x2) == i1:
        x2 += 1
    i2 = get_trace_keys(segy, iline, x2)

    if i2 - i1 != i1 - i0:
        return 0
    else:
        return i2 - i1


def guess(segy_name: str or Pysegy) -> List:
    """
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
    if isinstance(segy_name, str):
        segy = Pysegy(segy_name)
    elif isinstance(segy_name, Pysegy):
        segy = segy_name
    else:
        raise TypeError("Invalid type of `segy_name`")
    xlines = eval_xline(segy)
    ilines = eval_iline(segy)
    xselect = []
    iselect = []
    xsteps = []
    isteps = []

    for xline in xlines:
        xstep = eval_xstep(segy, xline)
        if xstep:
            xselect.append(xline)
            xsteps.append(xstep)

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
                s = Pysegy(segy_name)
                s.setInlineLocation(iselect[i])
                s.setCrosslineLocation(xselect[x])
                s.setSteps(isteps[i], xsteps[x])
                t = s.metaInfo()
            except:
                continue

            out.append([iselect[i], xselect[x], isteps[i], xsteps[x]])

    if out == []:
        raise RuntimeError("cannot guess the location and steps")

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
    elif metainfo.data_format == 5:
        out['dtype'] = '>4f-ieee'
    else:
        raise TypeError("don't support this data format")

    out['scalar'] = metainfo.scalar
    zinterval = metainfo.Z_interval
    yinterval = metainfo.Y_interval
    if apply_scalar:
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
