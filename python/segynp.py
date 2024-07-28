# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import Dict, List, Tuple
import numpy as np
from .cigsegy import (Pysegy, disable_progressbar) # type: ignore
import utils


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

    def __init__(self, filename, iline=None, xline=None, istep=None, xstep=None, as_2d=False) -> None:
        # TODO: guess
        disable_progressbar()
        self.as_3d = not as_2d
        self.fname = filename

        self.segy = Pysegy(filename)
        if self.as_3d:
            [iline, xline, istep, xstep, xloc, yloc] = utils.guess(filename, iline, xline, istep, xstep, 181, 185)
            self.segy.setInlineLocation(iline)
            self.segy.setCrosslineLocation(xline)
            self.segy.setSteps(istep, xstep)
            self.segy.scan()

        self.metainfo = self.segy.get_metaInfo()
        if self.as_3d:
            self._shape = (self.metainfo.sizeZ, self.metainfo.sizeY, self.metainfo.sizeX)
        else:
            self._shape = (self.segy.trace_count, self.metainfo.sizeX)
        self._eval_range()

    def _eval_range(self):
        s, e = self.trace_count-4000, 4000
        if self.trace_count < 4000:
            e = self.trace_count
            s = 0
        p0 = self.segy.collect(0, e)
        mi, ma = p0.min(), p0.max()

        p0 = self.segy.collect(s, self.trace_count)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        p0 = self.segy.collect(self.trace_count//3, self.trace_count//3+4000)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        self._min = mi
        self._max = ma

    @property
    def trace_count(self) -> int:
        return self.segy.trace_count

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
        if self.as_3d:
            return self.read(*idx)
        else:
            data = self.segy.collect(idx[0], idx[1])
            data = np.squeeze(data[:, idx[2]:idx[3]])
            if data.ndim == 0:
                return float(data)
            return data

    def _process_keys(self, key) -> List:
        if isinstance(key, int):
            if key < 0:
                key += self.shape[0]
            if key < 0 or key >= self.shape[0]:
                raise IndexError("Index out of range")
            if self.as_3d:
                return key, key + 1, 0, self.shape[1], 0, self.shape[2]
            else:
                return key, key + 1, 0, self.shape[1]

        elif key is Ellipsis:
            if self.as_3d:
                return 0, self.shape[0], 0, self.shape[1], 0, self.shape[2]
            else:
                return 0, self.shape[0], 0, self.shape[1]

        elif isinstance(key, Tuple):
            N = 3 if self.as_3d else 2
            num_ellipsis = key.count(Ellipsis)
            if num_ellipsis > 1:
                raise ValueError("Only one ellipsis (...) allowed")
            elif num_ellipsis == 1:
                key = (key[0], slice(None, None, None), slice(None, None, None))
                if not self.as_3d:
                    key = (key[0], slice(None, None, None))

            start_idx = [None] * N
            end_idx = [None] * N
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

            for i in range(N):
                if start_idx[i] is None:
                    start_idx[i] = 0
                if end_idx[i] is None:
                    end_idx[i] = self.shape[i]

            if self.as_3d:
                ib, xb, tb = start_idx
                ie, xe, te = end_idx
                return ib, ie, xb, xe, tb, te
            else:
                ib, xb = start_idx
                ie, xe = end_idx
                return ib, ie, xb, xe

        elif isinstance(key, slice):
            ib = 0 if key.start is None else key.start
            ie = self._shape[0] if key.stop is None else key.stop
            if self.as_3d:
                return ib, ie, 0, self.shape[1], 0, self.shape[2]
            else:
                return ib, ie, 0, self.shape[1]
        else:
            raise IndexError("Invalid index slices")


    def _check_bound(self, ib, ie, xb, xe, tb=None, te=None) -> None:
        assert ib < ie and ie <= self._shape[0]
        assert xb < xe and xe <= self._shape[1]
        if self.as_3d:
            assert (tb is not None) and (te is not None)
            assert tb < te and te <= self._shape[2]

    def max(self) -> float:
        return self._max

    def min(self) -> float:
        return self._min

    def close(self) -> None:
        self.segy.close_file()