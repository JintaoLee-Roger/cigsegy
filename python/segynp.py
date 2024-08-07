# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import Dict, List, Tuple
import numpy as np
from .cigsegy import (Pysegy, disable_progressbar)  # type: ignore
from . import utils
from .transform import get_transform_metrix, apply_transform
import warnings


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
    >>> d.to_2d() # as a collection of traces, shape is like (trace_count, n-time)
    >>> k = d[900] # the 900-th traces, output shape is (nt, )
    """

    def __init__(self,
                 filename,
                 iline=None,
                 xline=None,
                 istep=None,
                 xstep=None,
                 xloc=None,
                 yloc=None,
                 as_2d=False) -> None:
        disable_progressbar()
        np.set_printoptions(suppress=True)
        self.as_3d = not as_2d
        self.fname = filename

        self.segy = Pysegy(str(filename))
        self.keylocs = None
        self._shape3d = None
        if self.as_3d:
            self._scan3d(iline, xline, istep, xstep, xloc, yloc)

        self.metainfo = self.segy.get_metaInfo()
        if self.as_3d:
            self._shape3d = (self.metainfo.sizeZ, self.metainfo.sizeY, self.metainfo.sizeX)

        self._shape2d = (self.segy.trace_count, self.metainfo.sizeX)
        self._eval_range()

        if self.metainfo.scalar == 0:
            self._scalar = 1
        elif self.metainfo.scalar < 1:
            self._scalar = -1 / self.metainfo.scalar
        else:
            self._scalar = self.metainfo.scalar

        # for coordinates transform
        self._trans_matrix = None
        self._geomety = None

    def _scan3d(self, iline=None, xline=None, istep=None, xstep=None, xloc=None, yloc=None):
        [iline, xline, istep, xstep, xloc, yloc] = utils.guess(self.fname, iline, xline, istep, xstep, xloc, yloc)[0]
        self.keylocs = [iline, xline, istep, xstep, xloc, yloc]
        self.segy.setInlineLocation(iline)
        self.segy.setCrosslineLocation(xline)
        self.segy.setSteps(istep, xstep)
        self.segy.scan()

    def _eval_range(self):
        s, e = self.trace_count - 4000, 4000
        if self.trace_count < 4000:
            e = self.trace_count
            s = 0
        p0 = self.segy.collect(0, e)
        mi, ma = p0.min(), p0.max()

        p0 = self.segy.collect(s, self.trace_count)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        p0 = self.segy.collect(self.trace_count // 3, self.trace_count // 3 + 4000)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        self._min = mi
        self._max = ma

    @property
    def trace_count(self) -> int:
        return self.segy.trace_count

    @property
    def scalar(self) -> float:
        return self._scalar

    @scalar.setter
    def scalar(self, value):
        self._scalar = value

    @property
    def shape(self) -> Tuple:
        if self.as_3d:
            return self._shape3d
        else:
            return self._shape2d

    def _read(self, ib, ie, xb, xe, tb, te) -> np.ndarray:
        shape = [ie - ib, xe - xb, te - tb]
        d = self.segy.read(ib, ie, xb, xe, tb, te)

        if shape.count(1) == 3:
            return d[0, 0, 0]
        return np.squeeze(d)

    def __getitem__(self, slices) -> np.ndarray:
        idx = self._process_keys(slices)
        self._check_bound(*idx)
        if self.as_3d:
            return self._read(*idx)
        else:
            data = self.segy.collect(idx[0], idx[1])
            data = np.squeeze(data[:, idx[2]:idx[3]])
            if data.ndim == 0:
                return float(data)
            return data

    def _process_keys(self, key) -> List:
        if isinstance(key, (int, np.integer)):
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
                if isinstance(k, (int, np.integer)):
                    if k < 0:
                        k += self.shape[i]
                    start_idx[i] = k
                    end_idx[i] = k + 1
                elif isinstance(k, slice):
                    if not (k.step is None or k.step == 1):
                        raise IndexError(f"only support step is 1, while got a step {k.step} in the {i}th dimension")

                    start_idx[i] = k.start or 0
                    end_idx[i] = k.stop or self.shape[i]
                elif isinstance(k, (List, np.ndarray)):
                    raise NotImplementedError("Not implemented yet: TODO:") 
                else:
                    raise IndexError("Invalid index slices")

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
            ie = self.shape[0] if key.stop is None else key.stop
            if self.as_3d:
                return ib, ie, 0, self.shape[1], 0, self.shape[2]
            else:
                return ib, ie, 0, self.shape[1]
        elif isinstance(key, (List, np.ndarray)):
            raise NotImplementedError("Not implemented yet: TODO:") 
        else:
            raise IndexError("Invalid index slices")

    def _check_bound(self, ib, ie, xb, xe, tb=None, te=None) -> None:
        assert ib < ie and ie <= self.shape[0]
        assert xb < xe and xe <= self.shape[1]
        if self.as_3d:
            assert (tb is not None) and (te is not None)
            assert tb < te and te <= self.shape[2]

    def to_2d(self) -> None:
        """Treat the SEG-Y file as a collection of traces, shape is like (trace_count, nt)"""
        self.as_3d = False

    def to_3d(self, iline=None, xline=None, istep=None, xstep=None) -> None:
        """
        Treat the SEG-Y file as a 3D array. If the SEG-Y file is scanned, the parameters (iline, xline, istep, xstep) will be ignored 
        """
        self.as_3d = True
        if self._shape3d is None:
            self._scan3d(iline, xline, istep, xstep)
            self.metainfo = self.segy.get_metaInfo()
            self._shape3d = (self.metainfo.sizeZ, self.metainfo.sizeY, self.metainfo.sizeX)
        elif iline or xline or istep or xstep:
            warnings.warn(
                "The file has been scanned, so `iline, xline, istep, xstep` will be ignored",
                DeprecationWarning,
                stacklevel=2)

    def max(self, real=False) -> float:
        if real:
            return self[...].min()
        return self._max

    def min(self, real=False) -> float:
        if real:
            return self[...].max()
        return self._min

    def close(self) -> None:
        self.segy.close_file()

    def __array__(self):
        """To support np.array(SegyNP(xxx))"""
        return self[...]

    def to_numpy(self):
        """like pandas"""
        return self[...]

    def __array_function__(self, func, types, args, kwargs):
        if func is np.nanmin:
            return self.min()
        elif func is np.nanmax:
            return self.max()
        raise NotImplementedError(f"Function {func} is not implemented for SegyNP")

    @property
    def iline_range(self):
        if self._shape3d is None:
            raise NotImplementedError("Need scan the file first, please call `to_3d` first") # HACK: need to be optimized when as_2d
        return self.metainfo.min_inline, self.metainfo.max_inline

    @property
    def xline_range(self):
        if self._shape3d is None:
            raise NotImplementedError("Need scan the file first, please call `to_3d` first") # HACK: need to be optimized when as_2d
        return self.metainfo.min_crossline, self.metainfo.max_crossline

    @property
    def time_start(self):
        """
        Get the start time of the SEG-Y file, unit in ms
        """
        return self.metainfo.start_time

    @property
    def dt(self):
        """
        Get the sample interval of the SEG-Y file, unit in ms
        """
        return self.metainfo.sample_interval / 1000

    @property
    def interval(self):
        """
        Get the interval of the SEG-Y file, (inline_interval, xline_interval, dt)
        """
        di = np.round(self.metainfo.Z_interval * self.scalar, 2)
        dx = np.round(self.metainfo.Y_interval * self.scalar, 2)
        return (di, dx, self.dt)

    def __del__(self):
        self.close()

    def lineinfo(self):
        """
        Get the line information of the SEG-Y file,
        each line represents: inline, crossline_start, crossline_end, trace_start, trace_end, count
        """
        if self._shape3d is None:
            raise NotImplementedError("Need scan the file first, please call `to_3d` first")
        return self.segy.get_lineInfo()

    def update_trans_matrix(self):
        if self._shape3d is not None:
            line = self.lineinfo()
            ni = line.shape[0]
            nx = self._shape3d[1]
            xyic = np.zeros((ni*2+nx*2, 4))
            keylocs = [self.keylocs[4], self.keylocs[5], self.keylocs[0], self.keylocs[1]]
            for i in range(ni):
                xyic[i*2, :] = self.segy.get_trace_keys(keylocs, [4]*4, line[i, 3], line[i, 3]+1)
                xyic[i*2+1, :] = self.segy.get_trace_keys(keylocs, [4]*4, line[i, 4], line[i, 4]+1)
            xyic[ni*2:ni*2+nx, :] = self.segy.get_trace_keys(keylocs, [4]*4, 0, nx)
            xyic[ni*2+nx:ni*2+nx*2, :] = self.segy.get_trace_keys(keylocs, [4]*4, self.trace_count-nx, self.trace_count)
        else:
            raise NotImplementedError("Need scan the file first, please call `to_3d` first") # HACK: need to be optimized when as_2d

        self._trans_matrix = get_transform_metrix(xyic[:, 2:], xyic[:, :2])

    def xy_to_ix(self, xy, zero_origin=True):
        """
        Convert the x/y to inline/crossline

        Parameters
        ----------
        xy : np.ndarray
            The x/y array, shape is (n, 2)
        zero_origin : bool, optional
            Whether the x/y is zero-based, by default True,
            if is False, the x/y will be added by the min_x/min_y
        """
        xy = np.array(xy)
        shape = xy.shape
        if xy.ndim == 1:
            xy = xy.reshape(1, -1)
        if self._trans_matrix is None:
            self.update_trans_matrix()
        ic = apply_transform(xy, self._trans_matrix, inv=True)
        if zero_origin:
            ic[:, 0] -= self.iline_range[0]
            ic[:, 1] -= self.xline_range[0]
        return np.round(ic, 2).reshape(shape)

    def ix_to_xy(self, ix, zero_origin=True):
        """
        Convert the inline/crossline to x/y

        Parameters
        ----------
        ix : np.ndarray
            The inline/crossline array, shape is (n, 2)
        zero_origin : bool, optional
            Whether the inline/crossline is zero-based, by default True,
            if is False, the inline/crossline will be added by the min_inline/min_crossline
        """
        ix = np.array(ix)
        shape = ix.shape
        if ix.ndim == 1:
            ix = ix.reshape(1, -1)
        if zero_origin:
            ix[:, 0] += self.iline_range[0]
            ix[:, 1] += self.xline_range[0]
        if self._trans_matrix is None:
            self.update_trans_matrix()
        return np.round(apply_transform(ix, self._trans_matrix), 2).reshape(shape)


    def create_geometry(self):
        if self._shape3d is None:
            raise NotImplementedError("Need scan the file first, please call `to_3d` first")

        raise NotImplementedError("Not implemented yet: TODO:")
        geom = np.zeros(self.shape[:2]) - 1 # (ni, nx), -1 means no trace


    def read_traces_fast(self, index):
        raise NotImplementedError("Not implemented yet: TODO:")
