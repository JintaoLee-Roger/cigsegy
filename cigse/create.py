# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from datetime import datetime

import numpy as np
from .constinfo import kDataSampleFormatHelp
from cigse.cpp._CXX_SEGY import create_segy


def parser_textual(meta):
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
        steps = "iline: {meta['istep']:3}, xline: {meta['xstep']:3}"
    elif ndim == 4:
        shape = f"n-iline: {meta['ni']}, n-xline: {meta['nx']}, n-offset: {meta['no']}, n-time: {meta['nt']}"
        inter = f"di(iline) = {meta['di']:.2f} {meta['unit']}, dx(xline) = {meta['dx']:.2f} {meta['unit']}, dt = {meta['dt']//1000} ms"
        rge1 = f"inline: {meta['start_iline']} - {meta['end_iline']}, crossline: {meta['start_xline']} - {meta['end_xline']}"
        rge2 = f"offset: {meta['start_offset']} - {meta['end_offset']}, " + timer
        steps = f"iline: {meta['istep']:3}, xline: {meta['xstep']:3}, offset: {meta['ostep']:3}"

    dformat = kDataSampleFormatHelp[meta['dformat']]

    # fmt: off
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
        f"C29",
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

    def __init__(self, fname: str, data: np.ndarray, metainfo=None):
        self.filename = fname
        self.textual_header = None
        self.binary_header = None
        self.trace_header = None
        assert data.ndim >= 2 and data.ndim <= 4
        self.data = data
        self.shape = data.shape
        self.metainfo = metainfo
        self._ndim = len(self.shape)

    def set_textual_header(self, textual=None):
        pass

    def set_binary_header(self, bheader=None):
        pass

    def set_trace_header(self, theader=None):
        pass

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

    def copy_theader_from(self, segyname):
        self.trace_header = np.fromfile(
            segyname,
            dtype=np.uint8,
            count=240,
            offset=3600,
        )

    def create(self, keys: np.ndarray = None):
        if keys is None:
            keys = self._generate_keys()
        assert keys.ndim == 2
        if self._ndim < 4:
            assert keys.shape[1] == 4
        else:
            assert keys.shape[1] == 5

        assert keys.shape[0] == np.prod(self.shape[:-1])

        create_segy(self.fname, self.data, keys, self.textual_header,
                    self.binary_header, self.trace_header)

    def _generate_keys(self):
        pass

    def set_metaInfo(self):
        pass


def _gettime():
    return datetime.now().strftime("%Y-%m-%d %H:%M")


if __name__ == "__main__":
    import time
    print(_gettime())
