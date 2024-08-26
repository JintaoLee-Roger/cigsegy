# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from pathlib import Path
from typing import List
import numpy as np
# from segy import get_trace_keys
from segy.cpp._CXX_SEGY import Pysegy
from .constinfo import *


def ebcdic_to_ascii(ebcdic_str: bytes) -> str:
    """Convert EBCDIC encoded string to ASCII."""
    return ''.join(kEBCDICtoASCIImap.get(byte, '?') for byte in ebcdic_str)


def ascii_to_ebcdic(ascii_str: str) -> bytes:
    """Convert ASCII encoded string to EBCDIC."""
    return bytes(kASCIItoEBCDICmap.get(char, 0) for char in ascii_str)


def parse_bheader(bheader: np.ndarray):
    """
    Parse binary header.
    """
    assert bheader.size == 400, "Binary header size must be 400 bytes."
    out, hstring = _parse_header(bheader, kBinaryHeaderHelp)

    return out, hstring


def parse_theader(theader: np.ndarray):
    """
    Parse trace header.
    """
    assert theader.size == 240, "Trace header size must be 240 bytes."
    out, hstring = _parse_header(theader, kTraceHeaderHelp)

    return out, hstring


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
    ntrace = segy.ntrace
    options = [193, 17, 21, 13]
    select = [193, 17, 21, 13]
    for op in options:
        if _get_keys4(segy, op, ntrace - 1) - _get_keys4(segy, op, 0) > ntrace // 2: # yapf: disable
            select.remove(op)
            continue
        d = _get_keys4(segy, op, ntrace // 2, ntrace // 2 + 10)
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
    ntrace = segy.ntrace
    options = [189, 5, 9, 221]
    select = [189, 5, 9, 221]

    for op in options:
        l0 = _get_keys4(segy, op, 0)
        ll = _get_keys4(segy, op, ntrace - 1)
        l2 = _get_keys4(segy, op, ntrace // 2)
        if sum([x >= 0 for x in [l0, ll, l2]]) != 3:
            select.remove(op)
            continue
        if l0 == ll or l0 == l2 or ll == l2:
            select.remove(op)
            continue
        if max([l0, ll, l2]) - min([l0, ll, l2]) > min(ntrace // 10 - 1, 100000): # yapf: disable
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
    ntrace = segy.ntrace

    d = _get_keys4(segy, xline, ntrace // 2, ntrace // 2 + 10)
    diff = np.diff(d)
    if diff.min() > 0:
        return diff.min()
    idx = np.where(diff <= 0)[0]
    if len(idx) > 1:
        return 0

    d = _get_keys4(segy, xline, ntrace // 2 + idx[0] + 1,
                   ntrace // 2 + 10 + idx[0] + 1)
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
    ntrace = segy.ntrace
    i0 = _get_keys4(segy, iline, ntrace // 2)
    x1 = ntrace // 2 + 1
    while _get_keys4(segy, iline, x1) == i0:
        x1 += 1
    i1 = _get_keys4(segy, iline, x1)

    x2 = x1 + 1
    while _get_keys4(segy, iline, x2) == i1:
        x2 += 1
    i2 = _get_keys4(segy, iline, x2)

    if i2 - i1 != i1 - i0:
        return 0
    else:
        return i2 - i1


def guess(segy_name: str,
          iline=None,
          xline=None,
          offset=None,
          istep=None,
          xstep=None,
          ostep=None,
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
        raise RuntimeError(
            "cannot define the geometry through the location and steps")

    if xloc is None or yloc is None:
        s = Pysegy(str(segy_name))
        s.setInlineLocation(out[0][0])
        s.setCrosslineLocation(out[0][1])
        s.setSteps(out[0][2], out[0][3])
        s.setXLocation(181)
        s.setYLocation(185)
        s.scan()
        mt = s.get_metaInfo()
        if np.isnan(mt.Z_interval) or np.isnan(
                mt.Y_interval) or mt.Z_interval == 0 or mt.Y_interval == 0:
            xloc, yloc = 73, 77
        else:
            xloc, yloc = 181, 185
        s.close_file()

    for o in out:
        o += [xloc, yloc]

    return out


def eval_offset(segyname, offset=37):
    if isinstance(segyname, Pysegy):
        segy = segyname
    else:
        segy = Pysegy(str(segyname))
    ntrace = segy.ntrace
    ks = segy.get_trace_keys([offset], [4], ntrace // 3, ntrace // 3 + 100)
    if len(np.unique(ks)) == 1:
        return -1
    # TODO: deal with decreasing offset
    dif = np.diff(ks)
    dif = dif[dif > 0]
    return dif.min()


def parse_metainfo():
    pass


############# Internal functions #############


def _get_keys4(segy: Pysegy, keyloc, beg=0, end=-1):
    if beg < 0:
        beg = 0
        end = segy.ntrace
    if end == 0:
        end = beg + 1
    if end < 0:
        end = segy.ntrace

    if isinstance(keyloc, int):
        keyloc = [keyloc]
    d = segy.get_trace_keys(keyloc, [4] * len(keyloc), beg, end).squeeze()
    if d.size == 1 and d.ndim == 0:
        return int(d)
    elif d.size == 1 and d.ndim == 1:
        return int(d[0])
    else:
        return d


def _to_number(d, loc, ksize, dtype):
    return np.frombuffer(d[loc - 1:loc + ksize - 1].tobytes(), dtype=dtype)[0]


def _parse_header(header: np.ndarray, help_dict: dict) -> dict:
    """
    Parse binary header.
    """
    out = {}
    hstring = []
    for key, (disc, ksize) in help_dict.items():
        if ksize == 1:
            out[key] = header[key - 1]
        elif ksize == 2:
            out[key] = _to_number(header, key, ksize, '>i2')
        elif ksize == 4:
            out[key] = _to_number(header, key, ksize, '>i4')
        elif ksize == 8:
            out[key] = _to_number(header, key, ksize, '>i8')
        else:
            out[key] = 0
        hstring.append(f"{key}-{key+ksize}: {out[key]} - {disc}")

    return out, hstring
