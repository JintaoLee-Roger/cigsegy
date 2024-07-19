# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import Dict, List, Tuple
import numpy as np
from .cigsegy import (Pysegy, disable_progressbar, collect) # type: ignore


class SegyNP:
    """
    A segy file reader that mimics the numpy array style.

    Examples
    ---------
    >>> d = SegyNP('test.sgy')
    >>> d.shape
    # (601, 203, 400) # (n-inline, n-xline, n-time)
    >>> k = d[20:100, 100:200, 100:300] # k is np.array
    >>> k = d[20, :, :] # shape is (203, 400)
    >>> k = d[20, ...] # same as k = d[20, :, :]
    >>> k = d[20] # same as k = d[20, :, :]
    >>> k = d[:, 20, :]
    >>> k = f[20:30]
    >>> print(d.min(), d.max()) # NOTE: min and max are evalated by a small part of data
    """

    def __init__(self, filename, iline, xline, istep, xstep) -> None:
        # TODO: guess
        disable_progressbar()
        self.fname = filename
        self.segy = Pysegy(filename)
        self.segy.setInlineLocation(iline)
        self.segy.setCrosslineLocation(xline)
        self.segy.setSteps(istep, xstep)
        self.segy.scan()

        self.metainfo = self.segy.get_metaInfo()
        self._shape = (self.metainfo.sizeZ, self.metainfo.sizeY, self.metainfo.sizeX)
        self._eval_range()

    def _eval_range(self):
        tracecount = self.metainfo.trace_count
        s, e = tracecount-4000, 4000
        if tracecount < 4000:
            e = tracecount
            s = 0
        p0 = collect(self.fname, 0, e)
        mi, ma = p0.min(), p0.max()

        p0 = collect(self.fname, s, tracecount)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        p0 = collect(self.fname, tracecount//3, tracecount//3+4000)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        self._min = mi
        self._max = ma


    @property
    def shape(self) -> Tuple:
        return self._shape

    def read(self, ib, ie, xb, xe, tb, te) -> np.ndarray:
        shape = [ie - ib, xe - xb, te - tb]
        d = self.segy.read(ib, ie, xb, xe, tb, te)

        if shape.count(1) == 3:
            return d[0, 0, 0]
        return np.squeeze(d)


    def __getitem__(self, slices) -> np.ndarray:
        idx = self._process_keys(slices)

        self._check_bound(*idx)
        data = self.read(*idx)

        return data  # seismic index


    def _process_keys(self, key) -> List:
        if isinstance(key, int):
            if key < 0:
                key += self.shape[0]
            if key < 0 or key >= self.shape[0]:
                raise IndexError("Index out of range")
            return key, key + 1, 0, self.shape[1], 0, self.shape[2]
        elif key is Ellipsis:
            return 0, self.shape[0], 0, self.shape[1], 0, self.shape[2]
        elif isinstance(key, Tuple):
            num_ellipsis = key.count(Ellipsis)
            if num_ellipsis > 1:
                raise ValueError("Only one ellipsis (...) allowed")
            elif num_ellipsis == 1:
                key = (key[0], slice(None, None,
                                     None), slice(None, None, None))

            start_idx = [None] * 3
            end_idx = [None] * 3
            for i, k in enumerate(key):
                if k is None:
                    continue
                if isinstance(k, int):
                    if k < 0:
                        k += self.shape[i]
                    start_idx[i] = k
                    end_idx[i] = k + 1
                elif isinstance(k, slice):
                    if not (k.step is None or k.step == 1):
                        raise IndexError(f"only support step is 1, while got a step {k.step} in the {i}th dimension")

                    start_idx[i] = k.start or 0
                    end_idx[i] = k.stop or self.shape[i]

            for i in range(3):
                if start_idx[i] is None:
                    start_idx[i] = 0
                if end_idx[i] is None:
                    end_idx[i] = self.shape[i]

            ib, xb, tb = start_idx
            ie, xe, te = end_idx

            return ib, ie, xb, xe, tb, te
        elif isinstance(key, slice):
            ib = 0 if key.start is None else key.start
            ie = self._shape[0] if key.stop is None else key.stop
            return ib, ie, 0, self.shape[1], 0, self.shape[2]
        else:
            raise IndexError("Invalid index slices")

    def _check_bound(self, ib, ie, xb, xe, tb, te) -> None:
        assert ib < ie and ie <= self._shape[0]
        assert xb < xe and xe <= self._shape[1]
        assert tb < te and te <= self._shape[2]

    def max(self) -> float:
        return self._max

    def min(self) -> float:
        return self._min

    def close(self) -> None:
        self.segy.close_file()