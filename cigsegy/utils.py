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
    options = [189, 5, 9, 221, 13, 17, 41]
    select = [189, 5, 9, 221, 13, 17, 41]

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
    xlines = [193, 17, 21, 13, 45] if xline is None else [xline]
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
            return 0
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
        xloc, yloc = _guess_xy_locations(segy)
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


_PROFILE_KEY_SPECS = [
    ("trace_sequence_line", 1, 4),
    ("trace_sequence_file", 5, 4),
    ("field_record", 9, 4),
    ("trace_in_field_record", 13, 4),
    ("energy_source_point", 17, 4),
    ("ensemble", 21, 4),
    ("trace_in_ensemble", 25, 4),
    ("offset", 37, 4),
    ("elevation_or_shot_z", 41, 4),
    ("elevation_or_receiver_z", 45, 4),
    ("source_x", 73, 4),
    ("source_y", 77, 4),
    ("receiver_x", 81, 4),
    ("receiver_y", 85, 4),
    ("shot_line", 139, 2),
    ("shot_number", 141, 2),
    ("depth", 181, 4),
    ("iline", 189, 4),
    ("xline", 193, 4),
    ("shotpoint", 197, 4),
    ("extension_221", 221, 4),
]


def _normalize_profile_key_specs(key_specs):
    if key_specs is None:
        return list(_PROFILE_KEY_SPECS)

    out = []
    for spec in key_specs:
        if isinstance(spec, dict):
            name = spec.get("name", f"key_{spec['loc']}")
            out.append((name, int(spec["loc"]), int(spec.get("length", 4))))
        elif len(spec) == 2:
            loc, length = spec
            out.append((f"key_{loc}", int(loc), int(length)))
        else:
            name, loc, length = spec
            out.append((str(name), int(loc), int(length)))
    return out


def _mode_value(values):
    values, counts = np.unique(values, return_counts=True)
    return values[np.argmax(counts)], counts.max()


def _dominant_step(values):
    dif = np.diff(values)
    dif = dif[dif != 0]
    if dif.size == 0:
        return None
    step, _ = _mode_value(dif)
    return int(step)


def _range_count(vmin, vmax, step):
    if step is None or step == 0:
        return None
    span = int(vmax) - int(vmin)
    astep = abs(int(step))
    if astep == 0 or span < 0 or span % astep != 0:
        return None
    return span // astep + 1


def _estimate_period(values, step):
    if values.size < 3 or step == 0:
        return None

    repeated = np.flatnonzero(values == values[0])
    repeated = repeated[repeated >= 4]
    if repeated.size > 0:
        return int(repeated[0])

    dif = np.diff(values)
    if step > 0:
        resets = np.flatnonzero(dif < 0) + 1
    else:
        resets = np.flatnonzero(dif > 0) + 1

    if resets.size == 0:
        return None
    if resets.size == 1:
        period = int(resets[0])
        return period if period >= 4 else None
    periods = np.diff(np.r_[0, resets])
    period, _ = _mode_value(periods)
    period = int(period)
    return period if period >= 4 else None


def _score_fast_axis(values, loc):
    dif = np.diff(values)
    nonzero = dif[dif != 0]
    if nonzero.size == 0:
        return {
            "score": 0.0,
            "step": 0,
            "period": None,
            "unique": int(np.unique(values).size),
        }

    step, step_count = _mode_value(nonzero)
    step = int(step)
    step_ratio = float(step_count) / max(1, dif.size)
    unique_ratio = min(1.0, float(np.unique(values).size) / max(1, values.size))
    period = _estimate_period(values, step)
    reset_bonus = 1.0 if period is not None else 0.0

    priors = {
        193: 0.18,
        21: 0.12,
        13: 0.14,
        25: 0.10,
        37: 0.12,
        45: 0.04,
        181: 0.04,
    }
    score = 0.45 * step_ratio + 0.25 * unique_ratio + 0.25 * reset_bonus
    score += priors.get(loc, 0.0)
    if period is None:
        score *= 0.45
    score = min(1.0, score)
    vmin = int(values.min())
    vmax = int(values.max())
    range_values = values
    if period is not None and period <= values.size:
        range_values = values[:period]
        vmin = int(range_values.min())
        vmax = int(range_values.max())

    return {
        "score": float(score),
        "step": step,
        "period": period,
        "unique": int(np.unique(values).size),
        "min": vmin,
        "max": vmax,
        "range_count": _range_count(vmin, vmax, step),
    }


def _score_group_axis(values, period, loc):
    if period is None or period < 1:
        return {"score": 0.0, "unique": int(np.unique(values).size)}

    ngroup = values.size // period
    if ngroup < 2:
        return {"score": 0.0, "unique": int(np.unique(values).size)}

    trimmed = values[:ngroup * period].reshape(ngroup, period)
    within = np.array([np.unique(row).size == 1 for row in trimmed])
    group_values = trimmed[:, 0]
    group_unique = np.unique(group_values).size
    group_step = _dominant_step(group_values)
    change_ratio = min(1.0, float(group_unique - 1) / max(1, ngroup - 1))

    priors = {
        189: 0.18,
        221: 0.16,
        9: 0.14,
        5: 0.08,
        139: 0.08,
        141: 0.08,
        73: 0.04,
        77: 0.04,
    }
    score = 0.65 * float(within.mean()) + 0.30 * change_ratio
    score += priors.get(loc, 0.0)
    score = min(1.0, score)
    return {
        "score": float(score),
        "unique": int(group_unique),
        "step": group_step,
    }


def _relation_with_period(values, period):
    if period is None or period < 2:
        return {
            "constant_within": False,
            "monotonic_within": False,
            "changes_within": False,
        }

    ngroup = values.size // period
    if ngroup < 1:
        return {
            "constant_within": False,
            "monotonic_within": False,
            "changes_within": False,
        }

    rows = values[:ngroup * period].reshape(ngroup, period)
    sample = rows[:min(4, ngroup)]
    constant = [np.unique(row).size == 1 for row in sample]
    changes = [np.unique(row).size > 1 for row in sample]
    monotonic = []
    for row in sample:
        dif = np.diff(row)
        nz = dif[dif != 0]
        monotonic.append(nz.size > 0 and (np.all(nz > 0) or np.all(nz < 0)))

    return {
        "constant_within": bool(np.mean(constant) >= 0.75),
        "monotonic_within": bool(np.mean(monotonic) >= 0.75),
        "changes_within": bool(np.mean(changes) >= 0.75),
    }


def _profile_kind(fast, group, by_name, period):
    fast_loc = fast["loc"]
    group_loc = group["loc"] if group is not None else None

    depth_rel = _relation_with_period(by_name["depth"], period) if "depth" in by_name else {} # yapf: disable
    sx_rel = _relation_with_period(by_name["source_x"], period) if "source_x" in by_name else {} # yapf: disable
    sy_rel = _relation_with_period(by_name["source_y"], period) if "source_y" in by_name else {} # yapf: disable
    rx_rel = _relation_with_period(by_name["receiver_x"], period) if "receiver_x" in by_name else {} # yapf: disable
    ry_rel = _relation_with_period(by_name["receiver_y"], period) if "receiver_y" in by_name else {} # yapf: disable

    source_constant = sx_rel.get("constant_within", False) or sy_rel.get("constant_within", False) # yapf: disable
    receiver_changes = rx_rel.get("changes_within", False) or ry_rel.get("changes_within", False) # yapf: disable
    depth_monotonic = depth_rel.get("monotonic_within", False)

    if depth_monotonic and source_constant and receiver_changes:
        return "das_vsp"
    if fast_loc == 37:
        if group_loc == 21:
            return "cdp_gather"
        return "shot_gather"
    if fast_loc in (193, 21, 13, 17, 45) and group_loc in (189, 221, 9, 5, 41):
        return "volume_3d"
    if period is not None:
        return "gather_3d"
    return "trace_collection"


def guess_profile(segy_name,
                  key_specs=None,
                  sample_size: int = 8192,
                  top: int = 5):
    """
    Experimental profile inference for SEG-Y trace order.

    Unlike :func:`guess`, this function does not assume an inline/crossline
    volume. It first looks for a fast axis that changes within a gather/line
    and resets, then looks for a slow/group axis that is constant within that
    period. The result is a structured dictionary for inspection.

    This function is experimental and is not used by existing APIs.
    """
    if isinstance(segy_name, Pysegy):
        segy = segy_name
        need_close = False
    else:
        segy = Pysegy(str(segy_name))
        need_close = True

    try:
        specs = _normalize_profile_key_specs(key_specs)
        sample_size = max(2, min(int(sample_size), int(segy.ntrace)))
        keys = [loc for _, loc, _ in specs]
        lengths = [length for _, _, length in specs]
        data = segy.get_trace_keys(keys, lengths, 0, sample_size)
        if data.ndim == 1:
            data = data.reshape(-1, 1)

        columns = []
        by_name = {}
        for i, (name, loc, length) in enumerate(specs):
            values = np.asarray(data[:, i], dtype=np.int64)
            item = {
                "name": name,
                "loc": loc,
                "length": length,
                "values": values,
            }
            columns.append(item)
            by_name[name] = values

        fast_candidates = []
        for item in columns:
            stat = _score_fast_axis(item["values"], item["loc"])
            if stat["score"] <= 0:
                continue
            fast_candidates.append({
                "name": item["name"],
                "loc": item["loc"],
                "length": item["length"],
                **stat,
            })
        fast_candidates.sort(key=lambda x: x["score"], reverse=True)

        if not fast_candidates:
            return {
                "profile": "trace_collection",
                "confidence": 0.0,
                "ndim": 2,
                "shape": (int(segy.ntrace), int(segy.nt)),
                "axes": {
                    "trace": {"size": int(segy.ntrace)},
                    "time": {"size": int(segy.nt)},
                },
                "candidates": {"fast": [], "group": []},
            }

        fast = fast_candidates[0]
        period = fast["period"]
        group_candidates = []
        for item in columns:
            if item["loc"] == fast["loc"] and item["length"] == fast["length"]:
                continue
            stat = _score_group_axis(item["values"], period, item["loc"])
            if stat["score"] <= 0:
                continue
            group_candidates.append({
                "name": item["name"],
                "loc": item["loc"],
                "length": item["length"],
                **stat,
            })
        group_candidates.sort(key=lambda x: x["score"], reverse=True)
        group = group_candidates[0] if group_candidates else None

        profile = _profile_kind(fast, group, by_name, period)
        nfast = int(period) if period else None
        ngroup = None
        is_regular = None
        dense_trace_count = None
        shape = (int(segy.ntrace), int(segy.nt))
        ndim = 2
        shape_source = "trace"
        if profile == "volume_3d" and group is not None and fast["range_count"]:
            endpoints = segy.get_trace_keys(
                [group["loc"]],
                [group["length"]],
                np.array([0, segy.ntrace - 1], dtype=np.int32),
            ).reshape(-1)
            g0, g1 = int(endpoints[0]), int(endpoints[-1])
            gmin, gmax = min(g0, g1), max(g0, g1)
            ngroup = _range_count(gmin, gmax, group.get("step"))
            if ngroup is not None:
                nfast = int(fast["range_count"])
                dense_trace_count = int(ngroup * nfast)
                is_regular = dense_trace_count == int(segy.ntrace)
                shape = (ngroup, nfast, int(segy.nt))
                ndim = 3
                shape_source = "range"
                group["start"] = g0
                group["end"] = g1
                group["min"] = gmin
                group["max"] = gmax
        elif nfast and profile != "trace_collection":
            ngroup = int((segy.ntrace + nfast - 1) // nfast)
            dense_trace_count = int(ngroup * nfast)
            is_regular = dense_trace_count == int(segy.ntrace)
            shape = (ngroup, nfast, int(segy.nt))
            ndim = 3
            shape_source = "period"

        confidence = fast["score"]
        if group is not None:
            confidence = min(1.0, 0.55 * fast["score"] + 0.45 * group["score"])

        axes = {
            "fast": {
                k: fast[k]
                for k in ("name", "loc", "length", "step", "period", "unique",
                          "min", "max", "range_count")
            },
            "time": {"size": int(segy.nt), "dt": int(segy.bkeyi2(17))},
        }
        if group is not None:
            axes["group"] = {
                k: group[k]
                for k in ("name", "loc", "length", "unique", "step")
                if k in group
            }
            for k in ("start", "end", "min", "max"):
                if k in group:
                    axes["group"][k] = group[k]
        if ngroup is not None:
            axes["group_size"] = ngroup
            axes["fast_size"] = nfast

        out = {
            "profile": profile,
            "confidence": float(confidence),
            "ndim": ndim,
            "shape": shape,
            "shape_source": shape_source,
            "trace_shape": (int(segy.ntrace), int(segy.nt)),
            "ntrace": int(segy.ntrace),
            "nt": int(segy.nt),
            "axes": axes,
            "candidates": {
                "fast": fast_candidates[:top],
                "group": group_candidates[:top],
            },
        }
        if is_regular is not None:
            out["regular"] = is_regular
            out["dense_trace_count"] = dense_trace_count
            out["missing_or_partial_traces"] = dense_trace_count - int(segy.ntrace)
        return out
    finally:
        if need_close:
            segy.close()


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

    tracesort = "trace sorting code: " + kTraceSortingHelp.get(
        meta.get('trace_sorting_code', 0), "Unknown"
    )
    out += tracesort + "\n"

    dformat = f"scalar: {meta.get('scalar', 1)}, data format: " + kDataSampleFormatHelp.get(  # yapf: disable
        meta.get('dformat', 5), "Unknown"
    )
    out += dformat + "\n"

    kinfo = "(key info) "
    stepinfo = "           "
    if meta['ndim'] == 3:
        kinfo += f"iline: {meta['iline']:3}, xline: {meta['xline']:3}"
        stepinfo += f"istep: {meta['istep']:3}, xstep: {meta['xstep']:3}"
    if meta['ndim'] == 4:
        kinfo += f"iline: {meta['iline']:3}, xline: {meta['xline']:3}, offset: {meta['offset']:3}"
        stepinfo += f"istep: {meta['istep']:3}, xstep: {meta['xstep']:3}, ostep: {meta['ostep']:3}"
    elif meta['ndim'] == 2:
        kinfo += f"xloc: {meta.get('xloc', 181):3}, yloc: {meta.get('yloc', 185):3}"
    if meta['ndim'] > 2:
        kinfo += f", xloc: {meta['xloc']:3}, yloc: {meta['yloc']:3}"
    kinfo += "\n"
    out += kinfo
    out += stepinfo + "\n"

    return out


def post_process_meta(segy: Pysegy, meta: dict, apply_scalar=True):
    unit = segy.bkeyi2(55)
    if apply_scalar and 'di' in meta and 'dx' in meta:
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


def make_2d_meta(segy: Pysegy, xloc: int = None, yloc: int = None):
    if xloc is None or yloc is None:
        xloc, yloc = _guess_xy_locations(segy)

    dt = segy.bkeyi2(17)
    if dt <= 0:
        dt = segy.keyi2(0, 117)
    scalar = segy.keyi2(0, 71)
    if scalar == 0:
        scalar = 1

    meta = {
        'ndim': 2,
        'ntrace': segy.ntrace,
        'nt': segy.nt,
        'dt': dt,
        'start_time': segy.keyi2(0, 105),
        'trace_sorting_code': segy.bkeyi2(29),
        'dformat': segy.bkeyi2(25),
        'scalar': scalar,
        'xloc': xloc,
        'yloc': yloc,
    }
    return post_process_meta(segy, meta, False)


############# Internal functions #############


def _get_keys4(segy: Pysegy, keyloc, beg=-1, end=0):
    if beg < 0:
        beg = 0
        end = segy.ntrace
    if end < 0:
        end = segy.ntrace
    if end == 0:
        end = beg + 1

    if isinstance(keyloc, (int, np.integer)):
        keyloc = [keyloc]
    d = segy.get_trace_keys(keyloc, [4] * len(keyloc), beg, end).squeeze()
    if d.size == 1 and d.ndim == 0:
        return int(d)
    elif d.size == 1 and d.ndim == 1:
        return int(d[0])
    else:
        return d


def _guess_xy_locations(segy: Pysegy):
    sample_end = min(segy.ntrace, 100)
    xys = _get_keys4(segy, [181, 185, 73, 77], 0, sample_end)
    if xys.ndim == 1:
        xys = xys.reshape(1, -1)

    if len(np.unique(xys[:, 2])) == 1 and len(np.unique(xys[:, 3])) == 1:
        return 181, 185
    if len(np.unique(xys[:, 0])) == 1 and len(np.unique(xys[:, 1])) == 1:
        return 73, 77

    scalar = segy.keyi2(0, 71)
    if scalar < -1000 or scalar > 1000:
        scalar = 1
    scalar = 1 if scalar == 0 else scalar
    scalar = -1 / scalar if scalar < 0 else scalar
    xys = xys * scalar
    d1 = ((xys[-1, 0] - xys[0, 0])**2 + (xys[-1, 1] - xys[0, 1])**2)**0.5
    d2 = ((xys[-1, 2] - xys[0, 2])**2 + (xys[-1, 3] - xys[0, 3])**2)**0.5
    if d1 < 3 and d2 > 3:
        return 73, 77
    return 181, 185


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
