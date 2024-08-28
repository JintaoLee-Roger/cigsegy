# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.


def parser_textual():
    """
    """
    # fmt: off
    custom = [
        f"C01 Create By CIGSEGY software (CIG, USTC, 2024), see:",
        f"C02 github: https://github.com/JintaoLee-Roger/cigsegy",
        f"C03 Name: {'segyname'}",
        f"C04 Type: 3D seismic  Created Time: {'now'}",
        f"C05 ",
        f"C06",
        f"C07 ",
        f"C08",
        f"C09",
        f"C10",
        f"C11",
    ]
    metainfo = [
        f"C12 Inline range:    {'meta_info'} - {'meta_info'}",
        f"C13 Crossline range: {'meta_info'} - {'meta_info'}",
        f"C14 Inline step: {'meta_info'}, Crossline step: {'meta_info'}",
        f"C15 shape: (n-time, n-crossline, n-time) = ({'meta_info'}, {'meta_info'}, {'meta_info'})",
        f"C16 Number of samples per data trace: {'meta_info'}",
        f"C17 Sample interval: {'meta_info'} ",
        f"C18 Data sample format: {'dformat'}",
        f"C19 interval of inline: {'Z_interval':.2f} meters, interval of crossline: {'Y_interval':.2f} meters",
        f"C20 Time start: {'meta_info'}",
        f"C21 ",
        f"C22 ",
        f"C23",
    ]
    keyloc = [
        f"C24 Binary header locations:",
        f"C25 Sample interval             : bytes {17}",
        f"C26 Number of samples per trace : bytes {21}",
        f"C27 Data sample format code     : bytes {25}",
        f"C28 ",
        f"C29",
        f"C30 Trace header locations:",
        f"C31 Inline number               : bytes {189}",
        f"C32 Crossline number            : bytes {193}",
        f"C33 X coordinate                : bytes {181}",
        f"C34 Y coordinate                : bytes {185}",
        f"C35 Trace start time/depth      : bytes {105}",
        f"C36 ",
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

    def __init__(self, fname, data=None, shape=None, metainfo=None):
        self.filename = fname
        self.textual_header = None
        self.binary_header = None
        self.trace_header = None
        self.data = data
        self.shape = shape
        self.metainfo = metainfo

    def set_textual_header(self, textual=None):
        pass

    def set_binary_header(self, bheader=None):
        pass

    def set_trace_header(self, theader=None):
        pass

    def copy_textual_from(self, segyname):
        pass

    def copy_bheader_from(self, segyname):
        pass

    def copy_theader_from(self, segyname):
        pass

    def create(self, data, keys):
        pass
