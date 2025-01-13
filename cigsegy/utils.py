# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from pathlib import Path
from typing import List
import warnings
import numpy as np
# from cigsegy import get_trace_keys
from cigsegy.cpp._CXX_SEGY import Pysegy
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


def eval_iline(segy: Pysegy) -> int:
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
    options = [189, 5, 9, 221, 13, 17]
    select = [189, 5, 9, 221, 13, 17]

    for op in options:
        l0 = _get_keys4(segy, op, 0)
        ll = _get_keys4(segy, op, ntrace - 1)
        l2 = _get_keys4(segy, op, ntrace // 2)
        # line number should  not be negative
        if sum([x >= 0 for x in [l0, ll, l2]]) != 3:
            select.remove(op)
            continue
        # line number should be increasing/decreasing
        if l0 == ll or l0 == l2 or ll == l2:
            select.remove(op)
            continue
        # ni is too large
        if max([l0, ll, l2]) - min([l0, ll, l2]) > min(ntrace // 10 - 1, 50000): # yapf: disable
            select.remove(op)
            continue
        # line number should be increasing/decreasing
        if (l0 < l2 and l2 > ll) or (l0 > l2 and l2 < ll):
            select.remove(op)
            continue
        ls = _get_keys4(segy, op, ntrace // 2, ntrace // 2 + 50)
        if len(np.unique(ls)) > 10:
            select.remove(op)
            continue
        return op

    raise RuntimeError("Cannot evaluate inline location")


def is_valid_arr(arr, z_eval=False):
    arr_tuples = [tuple(row) for row in arr]
    if len(arr_tuples) != len(set(arr_tuples)):
        return False

    if z_eval:
        arr = arr[:, :2]
        z = arr[:, 2]
    x = arr[:, 0]
    y = arr[:, 1]

    if z_eval:
        ic = (z, y, x)
        dc = (-z, -y, -x)
    else:
        ic = (y, x)
        dc = (-y, -x)

    idx_asc = np.lexsort(ic)
    arr_asc = arr[idx_asc]

    idx_desc = np.lexsort(dc)
    arr_desc = arr[idx_desc]

    if np.array_equal(arr, arr_asc):
        return True
    elif np.array_equal(arr, arr_desc):
        return True
    else:
        return False


def guess2(segy_name: str,
           iline=None,
           xline=None,
           offset=None,
           istep=None,
           xstep=None,
           ostep=None,
           xloc=None,
           yloc=None):
    """
    guess the locations and steps of inline and crossline
    """
    if isinstance(segy_name, Pysegy):
        segy = segy_name
    else:
        segy = Pysegy(str(segy_name))

    N = segy.ntrace

    # 0. evaluate offset
    offset = 37 if offset is None else offset
    oostep, is4d = eval_offset(segy, offset)
    if ostep is None:
        oostep = ostep
    elif abs(oostep) < abs(ostep):
        warnings.warn(f"We evaluate istep == {istepx}, but you give istep == {istep}") # yapf: disable
    elif oostep * ostep < 0:
        warnings.warn(f"We evaluate istep == {istepx}, but you give istep == {istep}, order maybe wrong") # yapf: disable


    # 1. evaluate iline
    iline = eval_iline(segy) if iline is None else iline
    ni = _get_keys4(segy, iline, N - 1) - _get_keys4(segy, iline, 0)

    # 2. evaluate istep
    jumpi = N // ni // 4 if ni != 0 else 0
    if istep is None and ni == 0:
        istep = 1
    else:
        start = N // 4
        ils = [_get_keys4(segy, iline, start)]
        idx = [start]
        start += jumpi
        while len(ils) < 6 and start < N:
            il = _get_keys4(segy, iline, start)
            start += jumpi
            if il != idx[-1]:
                idx.append(start)
                ils.append(il)
        ilset = set(ils)
        ilset.update(_get_keys4(segy, iline, N // 4, N // 4 + 10).flatten())

        dif = np.diff(np.array(sorted(ilset)))
        istepx = dif.min() if dif[0] > 0 else dif.max()

        if istep is None:
            istep = istepx
        elif abs(istepx) < abs(istep):
            warnings.warn(f"We evaluate istep == {istepx}, but you give istep == {istep}") # yapf: disable
        elif istepx * istep < 0:
            warnings.warn(f"We evaluate istep == {istepx}, but you give istep == {istep}, order maybe wrong") # yapf: disable

    ni = ni / istep + 1
    if ni < 0 or not (isinstance(ni, int) or (isinstance(ni, float) and ni.is_integer())): # yapf: disable
        raise RuntimeError(f"Cannot evaluate iline/istep, becuase ni = {ni}")

    # 3. evaluate xline
    di = None
    if xline is None:
        candidate = [193, 17, 21, 13]
        start, end = N // 2, N // 2 + 50
        keys = candidate + [iline, offset]
        keys = keys if is4d else keys[:-1] # need offset?
        d = _get_keys4(segy, keys, start, end)
        for i, cxl in enumerate(candidate):
            mnx = min(N - 1, N // ni * 6)
            m1 = _get_keys4(segy, cxl, N - 1)
            m2 = _get_keys4(segy, cxl, 0)
            if abs(m1 - m2) > mnx:
                continue

            di = d[:, [4, i, 5]]
            while np.all(di[:, 0] == di[0, 0]) and np.all(di[:, 1] == di[0, 1]) and end < N // 2 + 1100: # yapf: disable
                start += 50
                end += 50
                d_new = _get_keys4(segy, [xline, candidate[i], offset], start, end) # yapf: disable
                di = np.concatenate([di, d_new])
            if end > 1000:
                continue

            if sum(di >= 0) != len(di):
                continue
            if not is_valid_arr(di):
                continue
            xline = cxl
            break

        if xline is None:
            raise RuntimeError("Evaluate xline error")

    # 4. evaluate xstep
    if xstep is None:
        if di is None:
            di = _get_keys4(segy, [iline, xline, offset], start, end)


    # 7. evaluate xloc, yloc
    if xloc is None or yloc is None:
        pass

    return iline, xline, offset, istep, xstep, ostep, xloc, yloc, is4d


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
    guess the locations and steps of inline and crossline
    """
    if isinstance(segy_name, Pysegy):
        segy = segy_name
    else:
        segy = Pysegy(str(segy_name))

    N = segy.ntrace
    offset = 37 if offset is None else offset
    # ostep, is4d = eval_offset(segy, offset)
    is4d = False
    ostep = 1

    iline = eval_iline(segy) if iline is None else iline

    # read 3 lines
    start = int(N // 3)
    lines = set()
    oix = []
    xlines = [193, 17, 21, 13] if xline is None else [xline]
    while True and (start + 400) <= segy.ntrace:
        part = _get_keys4(segy, [offset, iline, *xlines], start, start + 400)
        oix.append(part)
        lines.update(np.unique(part[:, 1]))
        start += 400
        if len(lines) >= 3:
            break

    # eval istep
    lines = np.array(sorted(list(lines)))
    # print(lines)
    dif = np.diff(lines)
    istepx = dif.min() if dif[0] > 0 else dif.max()

    # line: n, n+1, n+2, extract the data of line n+1
    oix = np.concatenate(oix)
    idx = np.where(np.diff(oix[:, 1]) != 0)[0][:2] + 1
    oix = oix[idx[0]:idx[1], :]
    nig = oix.shape[0]

    def _double_check_istep():
        start2 = int(N // 3 * 2)
        oix2 = []
        lines2 = set()
        while True and (start2 + nig) <= segy.ntrace:
            part2 = _get_keys4(segy, iline, start2, start2 + nig)
            oix2.append(part2)
            lines2.update(np.unique(part2))
            start2 += nig
            if len(lines2) >= 3:
                break
        lines2 = np.array(sorted(list(lines2)))
        dif2 = np.diff(lines2)
        istepx2 = dif2.min() if dif2[0] > 0 else dif2.max()
        return istepx2

    if abs(istepx) > 1:  # if istep is not 1, double check
        istepx2 = _double_check_istep()
        if abs(istepx2) < abs(istepx):
            istepx = istepx2

    if istep is not None and istep != istepx:
        warnings.warn(f"You set istep={istep}, but we scan the istep is {istepx}, be careful!") # yapf: disable
    else:
        istep = istepx

    # eval xline
    def _eval_xline(i):
        dif = np.diff(oix[:, i])
        idx = np.where(dif != 0)[0]
        # values, counts = np.unique(idx, return_counts=True)
        dif = dif[dif != 0]
        if len(dif) == 0:
            return 1
        xstepi = dif.min() if dif[0] > 0 else dif.max()
        return xstepi

    if len(xlines) == 1:
        xstepx = _eval_xline(2)
        xidx = 2
    else:
        for i in range(len(xlines)):
            if abs(_get_keys4(segy, xlines[i], N - 1) - _get_keys4(segy, xlines[i], 0)) > N // 10: # yapf: disable
                continue
            xstepx = _eval_xline(i + 2)
            if abs(xstepx) > 0 and abs(xstepx) < 100:
                xline = xlines[i]
                xidx = i + 2
                # oix = oix[:, [0, 1, i + 2]]
                break

    if xstep is not None and xstep != xstepx:
        warnings.warn(f"You set xstep={xstep}, but we scan the xstep is {xstepx}, be careful!") # yapf: disable
    elif xstepx == 0:
        raise RuntimeError("Cannot evaluate xline location, we evaluate xstep is 0") # yapf: disable
    else:
        xstep = xstepx

    # eval offset
    change_points = np.where(np.diff(oix[:, xidx]) != 0)[0] + 1
    segments = np.split(oix[:, xidx], change_points)
    count = sum(1 for seg in segments if len(seg) > 1)
    if count > 15:
        ostep, is4d = eval_offset(segy, offset)

    # eval xloc and yloc
    if xloc is None or yloc is None:
        xys = _get_keys4(segy, [181, 185, 73, 77], 0, 100)
        if len(np.unique(xys[:, 2])) == 0 and len(np.unique(xys[:, 3])) == 0:
            xloc, yloc = 181, 185
        elif len(np.unique(xys[:, 0])) == 0 and len(np.unique(xys[:, 1])) == 0:
            xloc, yloc = 73, 77
        else:
            scalar = segy.keyi2(0, 71)
            if scalar < -1000 or scalar > 1000:
                scalar = 1
            scalar = 1 if scalar == 0 else scalar
            scalar = -1 / scalar if scalar < 0 else scalar
            xys = xys * scalar
            d1 = ((xys[-1, 0] - xys[0, 0])**2 +
                  (xys[-1, 1] - xys[0, 1])**2)**0.5
            d2 = ((xys[-1, 2] - xys[0, 2])**2 +
                  (xys[-1, 3] - xys[0, 3])**2)**0.5
            if d1 < 3 and d2 > 3:
                xloc, yloc = 73, 77
            else:
                xloc, yloc = 181, 185
    return iline, xline, offset, istep, xstep, ostep, xloc, yloc, is4d


def eval_offset(segyname, offset: int = 37) -> int:
    """
    return the ostep, if return -1, the file is not prestack SEG-Y
    """
    if isinstance(segyname, Pysegy):
        segy = segyname
    else:
        segy = Pysegy(str(segyname))

    ntrace = segy.ntrace
    ks = _get_keys4(segy, offset, ntrace // 3, ntrace // 3 + 200)

    if not isinstance(segyname, Pysegy):
        segy.close()

    if len(np.unique(ks)) == 1:
        return 0, False

    dif = np.diff(ks)
    udif, count = np.unique(dif, return_counts=True)
    if (len(udif) > 10):
        warnings.warn("offset is not constant, and is unsorted, cannot evaluate by `scan`. So we treat this file as a 3D.") # yapf: disable
        # raise RuntimeError("offset is not constant, and is unsorted, cannot evaluate by `scan`") # yapf: disable
        return 0, False

    ostep = udif[np.argmax(count)]
    if ostep == 0:
        return 0, False
    return ostep, True


def parse_metainfo(meta: dict):
    out = ""

    # shape information
    shapeinfo = "shape: "
    if meta['ndim'] == 2:
        shapeinfo += f"(n-trace, n-time) = ({meta['ntrace']}, {meta['nt']})"
    elif meta['ndim'] == 3:
        shapeinfo += f"(n-inline, n-crossline, n-time) = ({meta['ni']}, {meta['nx']}, {meta['nt']})"
    else:
        shapeinfo += f"(n-inline, n-crossline, n-offset, n-time) = ({meta['ni']}, {meta['nx']}, {meta['no']}, {meta['nt']})"

    out += shapeinfo + "\n"
    out += f"N traces: {meta['ntrace']}\n"

    # interval
    intervalinfo = "interval: "
    if meta['ndim'] == 2:
        intervalinfo += f"dt = {meta['dt']//1000} ms"
    else:
        intervalinfo += f"di(iline) = {meta['di']:.2f} {meta['unit']}, dx(xline) = {meta['dx']:.2f} {meta['unit']}, dt = {meta['dt']//1000} ms"
    out += intervalinfo + "\n"

    # range
    rangeinfo = "range: "
    end_time = meta['start_time'] + (meta['nt'] - 1) * meta['dt'] / 1000
    timer = f"t: {meta['start_time']} - {end_time} ms"
    if meta['ndim'] > 2:
        rangeinfo += f"inline: {meta['start_iline']} - {meta['end_iline']}, crossline: {meta['start_xline']} - {meta['end_xline']}, "
    if meta['ndim'] == 4:
        rangeinfo += f"offset: {meta['start_offset']} - {meta['end_offset']}, "
    rangeinfo += timer
    out += rangeinfo + "\n"

    tracesort = "trace sorting code: " + kTraceSortingHelp[meta['trace_sorting_code']] # yapf: disable
    out += tracesort + "\n"

    dformat = f"scalar: {meta['scalar']}, data format: " + kDataSampleFormatHelp[meta['dformat']] # yapf: disable
    out += dformat + "\n"

    kinfo = "(key info) "
    stepinfo = "           "
    if meta['ndim'] == 3:
        kinfo += f"iline: {meta['iline']:3}, xline: {meta['xline']:3}"
        stepinfo += f"istep: {meta['istep']:3}, xstep: {meta['xstep']:3}"
    if meta['ndim'] == 4:
        kinfo += f", offset: {meta['offset']:3}"
        stepinfo += f"ostep: {meta['ostep']:3}"
    kinfo += f", xloc: {meta['xloc']:3}, yloc: {meta['yloc']:3}\n"
    out += kinfo
    out += stepinfo + "\n"

    return out


def post_process_meta(segy: Pysegy, meta: dict, apply_scalar=True):
    unit = segy.bkeyi2(55)
    if apply_scalar:
        if meta['scalar'] == 0:
            meta['scalar'] = 1
        scalar = -1 / meta['scalar'] if meta['scalar'] < 0 else meta['scalar']
        meta['di'] *= scalar
        meta['dx'] *= scalar
    if unit == 2:
        meta['unit'] = 'ft'
    else:
        meta['unit'] = 'm'
    return meta


############# Internal functions #############


def _get_keys4(segy: Pysegy, keyloc, beg=-1, end=0):
    if beg < 0:
        beg = 0
        end = segy.ntrace
    if end < 0:
        end = segy.ntrace
    if end == 0:
        end = beg + 1

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
        hstring.append(f"{key:^3} - {key+ksize-1:^3}: {out[key]:<8} - {disc}")

    return out, hstring
