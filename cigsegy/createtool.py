# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from datetime import datetime

import numpy as np
from .constinfo import kDataSampleFormatHelp
from cigsegy.cpp._CXX_SEGY import create_segy
from cigsegy import utils


def assemble_metainfo(shape, dformat=5, start=None, interval=None):
    ndim = len(shape)
    assert ndim >= 2 and ndim <= 4
    meta = {}
    meta['ndim'] = ndim
    meta['unit'] = 'm'
    meta['dformat'] = dformat
    # shape
    meta['nt'] = shape[-1]
    if ndim == 2:
        meta['ntrace'] = shape[0]
    else:
        meta['ni'] = shape[0]
        meta['nx'] = shape[1]
        if ndim == 4:
            meta['no'] = shape[2]
        meta['ntarce'] = np.prod(shape[:-1])

    # start
    if ndim == 2:
        if start is None:
            meta['start_time'] = 0
        else:
            assert isinstance(start, (int, np.integer))
            meta['start_time'] = start
    else:
        if start is None:
            start = [2000] * (len(shape) - 1)
            start.append(0)
        assert len(start) == ndim
        meta['start_iline'] = start[0]
        meta['end_iline'] = start[0] + meta['ni'] - 1
        meta['start_xline'] = start[1]
        meta['end_xline'] = start[1] + meta['nx'] - 1
        meta['start_time'] = start[-1]
        if ndim == 4:
            meta['start_offset'] = start[1]
            meta['end_offset'] = start[1] + meta['nx'] - 1

    # interval
    if ndim == 2:
        if interval is None:
            meta['dt'] = 2000
        else:
            assert isinstance(interval, (int, np.integer))
            meta['dt'] = interval
    else:
        if interval is None:
            interval = [25, 25, 2000]
        assert len(interval) == 3
        meta['di'] = interval[0]
        meta['dx'] = interval[1]
        meta['dt'] = interval[2]

    meta['istep'] = 1
    meta['xstep'] = 1
    meta['ostep'] = 1
    meta['scalar'] = 1

    return meta


def generate_textual(meta, textual=None):
    if textual is None:
        textual = parser_textual(meta)
    else:
        if isinstance(textual, list):
            if len(textual) <= 8:
                textual = parser_textual(meta, textual)
            else:
                textual = _parser_textual_list(textual)
        else:
            if len(textual) < 76 * 8:
                textual = parser_textual(meta, textual)
    assert len(textual) == 3200
    textual_header = utils.ascii_to_ebcdic(textual)
    return textual_header


def parser_textual(meta, custom=None):
    """
    """
    t = _gettime()
    ndim = meta['ndim']
    end_time = meta['start_time'] + (meta['nt'] - 1) * meta['dt'] // 1000
    timer = f"t: {meta['start_time']} - {end_time} ms"
    if ndim == 2:
        shape = f"n-trace: {meta['ntrace']}, n-time: {meta['nt']}"
        inter = f"dt = {meta['dt']//1000} ms"
        rge1 = timer
        rge2 = ""
        steps = ""
    elif ndim == 3:
        shape = f"n-iline: {meta['ni']}, n-xline: {meta['nx']}, n-time: {meta['nt']}"
        inter = f"di(iline) = {meta['di']:.2f} {meta['unit']}, dx(xline) = {meta['dx']:.2f} {meta['unit']}, dt = {meta['dt']//1000} ms"
        rge1 = f"inline: {meta['start_iline']} - {meta['end_iline']}, crossline: {meta['start_xline']} - {meta['end_xline']}"
        rge2 = timer
        steps = f"iline: {meta['istep']:3}, xline: {meta['xstep']:3}"
    elif ndim == 4:
        shape = f"n-iline: {meta['ni']}, n-xline: {meta['nx']}, n-offset: {meta['no']}, n-time: {meta['nt']}"
        inter = f"di(iline) = {meta['di']:.2f} {meta['unit']}, dx(xline) = {meta['dx']:.2f} {meta['unit']}, dt = {meta['dt']//1000} ms"
        rge1 = f"inline: {meta['start_iline']} - {meta['end_iline']}, crossline: {meta['start_xline']} - {meta['end_xline']}"
        rge2 = f"offset: {meta['start_offset']} - {meta['end_offset']}, " + timer
        steps = f"iline: {meta['istep']:3}, xline: {meta['xstep']:3}, offset: {meta['ostep']:3}"

    dformat = f"{meta['dformat']} -- " + kDataSampleFormatHelp[meta['dformat']]

    # fmt: off
    if custom is None:
        custom = [
            f"C01 Create By CIGSEGY software (CIG, USTC, 2024), see:",
            f"C02 github: https://github.com/JintaoLee-Roger/cigsegy",
            f"C03 ",
            f"C04 Type: {ndim}D seismic  Created Time: {t}",
            f"C05 ",
            f"C06 ",
            f"C07 ",
            f"C08 ",
        ]
    else:
        if isinstance(custom, list):
            assert len(custom) <= 8
            if len(custom) < 8:
                custom += [''] * (8 - len(custom))
            custom = [f'C{i+1:02} ' + s[:76] for i, s in enumerate(custom)]
        else:
            out = []
            for i in range(8):
                start = i * 76
                if i * 76 >= len(custom):
                    out.append(f'C{i+1:02}')
                else:
                    end = (i+1) * 76
                    end = end if end < len(custom) else len(custom)
                    out.append(f'C{i+1:02} ' + custom[start:end])
            custom = out

    metainfo = [
        f"C09 ",
        f"C10 ",
        f"C11 Shape: {shape}",
        f"C12 Interval: {inter}",
        f"C13 Ranges: {rge1}",
        f"C14         {rge2}",
        f"C15 Line steps: {steps}",
        f"C16 ",
        f"C17 scalar to coordinates: {meta['scalar']}",
        f"C18 Data sample format: {dformat}",
        f"C19 ",
        f"C20 ",
        f"C21 ",
        f"C22 ",
        f"C23",
    ]
    keyloc = [
        f"C24 Binary header locations:",
        f"C25 Sample interval             : bytes 17",
        f"C26 Number of samples per trace : bytes 21",
        f"C27 Data sample format code     : bytes 25",
        f"C28 ",
        f"C29 ",
        f"C30 Trace header locations:",
        f"C31 Inline number               : bytes 189",
        f"C32 Crossline number            : bytes 193",
        f"C33 X coordinate                : bytes 181",
        f"C34 Y coordinate                : bytes 185",
        f"C35 Trace start time/depth      : bytes 105",
        f"C36 Offset                      : bytes 37",
        f"C37 ",
        f"C38 ",
        f"C39 ",
        f"C40 ",
    ]

    textual = custom + metainfo + keyloc

    if len(textual) != 40:
        raise ValueError(f"Header must contain exactly 40 lines, got {len(textual)}")

    textual = ''.join(f'{line.ljust(80)}' for line in textual)

    if len(textual) != 3200:
        raise ValueError(f"Header length must be exactly 3200 bytes, got {len(textual)}")

    return textual
    # fmt: on


class SegyCreate:
    """
    Create a SEG-Y file from an array.

    Parameters
    -----------
    fname : str
        out SEG-Y file name
    data : np.ndarray
        array
    metainfo : None
        meta information

    Examples
    ---------
    ```python
    # data.shape == (1101, 801, 1201)
    segyc = SegyCreate('out.segy', data)

    # 1. default
    segyc.create()

    # 2. set start_time, dt
    segyc.set_dt(4000) # 4 ms
    segyc.set_start_time(200) # 200 ms
    segyc.create()

    # 3. manually set iline, xline for each trace
    iline = np.arange(2000, 2000+data.shape[0])
    xline = np.arange(3000, 3000+data.shape[1])
    il, xl = np.meshgrid(iline, xline, indexing='ij')
    keys = np.c_[il.flatten(), xl.flatten()]
    segyc.create(keys)

    # 4. manually set iline, xline, x, y for each trace
    iline = np.arange(2000, 2000+data.shape[0])
    xline = np.arange(3000, 3000+data.shape[1])
    x = np.arange(605959, 605959 + data.shape[0]*25)
    y = np.arange(6058307, 6058307 + data.shape[1]*12)
    il, xl = np.meshgrid(iline, xline, indexing='ij')
    x, y = np.meshgrid(x, y, indexing='ij')
    keys = np.c_[il.flatten(), xl.flatten(), x.flatten(), y.flatten()]
    segyc.create(keys)
    ```
    """

    def __init__(self, fname: str, data: np.ndarray, metainfo=None):
        self.filename = fname
        self.textual_header = None
        self.binary_header = None
        self.trace_header = None
        assert data.ndim >= 2 and data.ndim <= 4
        self.data = data
        self.shape = data.shape
        if metainfo is None:
            metainfo = assemble_metainfo(data.shape)
        self.metainfo = metainfo
        self._ndim = len(self.shape)

    def set_textual_header(self, textual=None):
        self.textual_header = generate_textual(self.metainfo, textual)

    def set_binary_header(self, bheader: np.ndarray = None):
        if bheader is None:
            bheader = self._init_bheader()
        bheader = bheader.view(np.uint8)
        assert bheader.size == 400
        self._set_keyi2(bheader, 17, self.metainfo['dt'])
        self._set_keyi2(bheader, 21, self.metainfo['nt'])
        self._set_keyi2(bheader, 25, self.metainfo['dformat'])
        self.binary_header = bheader

    def set_trace_header(self, theader: np.ndarray = None):
        if theader is None:
            theader = self._init_theader()
        theader = theader.view(np.uint8)
        assert theader.size == 240
        self._set_keyi2(theader, 69, self.metainfo['scalar'])
        self._set_keyi2(theader, 115, self.metainfo['nt'])
        self._set_keyi2(theader, 117, self.metainfo['dt'])
        self._set_keyi2(theader, 105, self.metainfo['start_time'])
        self.trace_header = theader

    def set_dt(self, dt):
        self.metainfo['dt'] = dt
        if self.binary_header is not None:
            self._set_keyi2(self.binary_header, 17, dt)
            self._set_keyi2(self.trace_header, 117, dt)

    def set_start_time(self, start_time):
        self.metainfo['start_time'] = start_time
        if self.trace_header is not None:
            self._set_keyi2(self.trace_header, 105, start_time)

    def copy_textual_from(self, segyname):
        with open(segyname, 'r') as f:
            self.textual_header = str(f.read(3200))

    def copy_bheader_from(self, segyname):
        self.binary_header = np.fromfile(
            segyname,
            dtype=np.uint8,
            count=400,
            offset=3200,
        )
        self._set_keyi2(self.binary_header, 21, self.shape[-1])

    def copy_theader_from(self, segyname):
        self.trace_header = np.fromfile(
            segyname,
            dtype=np.uint8,
            count=240,
            offset=3600,
        )
        self._set_keyi2(self.trace_header, 115, self.shape[-1])

    def create(self, keys: np.ndarray = None):
        if keys is None:
            keys = self._generate_keys()
        assert keys.ndim == 2
        if (self._ndim < 4 and keys.shape[1] == 2) or (self._ndim == 4
                                                       and keys.shape[1] == 3):
            xy = self._generate_xy()
            assert xy.shape[0] == keys.shape[0]
            keys = np.concatenate([keys, xy], axis=1)
        if self._ndim < 4:
            assert keys.shape[1] == 4
        else:
            assert keys.shape[1] == 5

        assert keys.shape[0] == np.prod(self.shape[:-1])

        if self.textual_header is None:
            self.set_textual_header()
        if self.binary_header is None:
            self.set_binary_header()
        if self.trace_header is None:
            self.set_trace_header()

        create_segy(self.fname, self.data, keys, self.textual_header,
                    self.binary_header, self.trace_header)

    # fmt: off
    def _generate_keys(self):
        if self._ndim == 2:
            x, y = np.meshgrid(np.array([10]), np.arange(100, self.shape[0]+100), indexing='ij')
            keys = np.c_[x.flatten(), y.flatten()]
        else:
            ib = self.metainfo['start_iline']
            ie = ib + self.shape[0]
            xb = self.metainfo['start_xline']
            xe = xb + self.shape[1]
            if self._ndim == 3:
                x, y = np.meshgrid(np.arange(ib, ie, 1), np.arange(xb, xe, 1), indexing='ij')
                keys = np.c_[x.flatten(), y.flatten()]
            else:
                ob = self.metainfo['start_offset']
                oe = ob + self.shape[2]
                x, y, z = np.meshgrid(np.arange(ib, ie, 1), np.arange(xb, xe, 1), np.arange(ob, oe, 1), indexing='ij')
                keys = np.c_[x.flatten(), y.flatten(), z.flatten()]
        return keys

    def _generate_xy(self):
        if self._ndim == 2:
            x, y = np.meshgrid(np.array([]), np.arange(), indexing='ij')
            xy = np.c_[x.flatten(), y.flatten()]
        else:
            di = self.metainfo['di']
            dx = self.metainfo['dx']
            ni = np.arange(108232, 108232+di*self.shape[0], di)
            nx = np.arange(413897, 413897+dx*self.shape[1], dx)
            if self._ndim == 3:
                x, y = np.meshgrid(x, y, indexing='ij')
                xy = np.c_[x.flatten(), y.flatten()]
            else:
                x, y, z = np.meshgrid(x, y, np.arange(10, self.shape[2]+10), indexing='ij')
                xy = np.c_[x.flatten(), y.flatten(), z.flatten()]

        return xy
    # fmt: on

    def _set_keyi2(self, header: np.ndarray, loc: int, value: int):
        loc = loc - 1
        header[loc:loc + 2] = np.array([value], dtype='>i2').view(np.uint8)
        return header

    def _set_keyi4(self, header: np.ndarray, loc: int, value: int):
        loc = loc - 1
        header[loc:loc + 4] = np.array([value], dtype='>i4').view(np.uint8)
        return header

    def _get_keyi2(self, header: np.array, loc: int):
        loc = loc - 1
        return np.frombuffer(header[loc:loc + 2], dtype='>i2')[0]

    def _get_keyi4(self, header: np.array, loc: int):
        loc = loc - 1
        return np.frombuffer(header[loc:loc + 4], dtype='>i4')[0]

    def _init_bheader(self):
        out = np.zeros((400, ), np.uint8)
        self._set_keyi4(out, 5, 1000)
        self._set_keyi2(out, 13, 1)
        self._set_keyi2(out, 55, 1)
        return out

    def _init_theader(self):
        out = np.zeros((240, ), np.uint8)
        if self._ndim == 2:
            self._set_keyi2(out, 29, 0)
        elif self._ndim == 3:
            self._set_keyi2(out, 29, 4)
        else:
            self._set_keyi2(out, 29, 2)

        return out


def _gettime():
    return datetime.now().strftime("%Y-%m-%d %H:%M")


def _parser_textual_list(text):
    textual = [f'C{i:02} ' + t for i, t in enumerate(text)]
    for i in range(len(text), 40):
        textual.append(f'C{i:02} ')

    textual = ''.join(f'{line.ljust(80)}' for line in textual)
    return textual
