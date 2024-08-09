# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import Dict, List, Tuple
import numpy as np
from .cigsegy import (Pysegy, disable_progressbar)  # type: ignore
from . import utils
from .transform import get_transform_metrix, apply_transform
from .interp import arbitray_line
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
    >>> d.to_3d()# back to 3d
    >>> line = d.arbitray_line([[0, 0], [100, 100]]) # arbitray line
    >>> d.xy_to_ic([555627, 998321]) # map X-Y to numpy array index
    >>> d.xy_to_ic([555627, 998321], zero_origin=False) # map X-Y to inline/xline index
    >>> d.ic_to_xy([10, 100]) # map numpy array index to X-Y
    >>> d = SegyNP('test.sgy', iline=189, xline=193, as_unsorted=True) # load as unsorted SEG-Y file
    """

    def __init__(
        self,
        filename: str,
        iline: int = None,
        xline: int = None,
        istep: int = None,
        xstep: int = None,
        xloc: int = None,
        yloc: int = None,
        as_2d: bool = False,
        as_unsorted: bool = False,
        ic=None,
    ) -> None:  # TODO: support 4D gather?
        disable_progressbar()
        np.set_printoptions(suppress=True)
        self.as_3d = not as_2d
        self.fname = filename

        self.segy = Pysegy(str(filename))
        self.keylocs = None
        self._shape3d = None
        self.metainfo = None

        # for coordinates transform
        self._trans_matrix = None
        self._geomety = None

        if as_2d and as_unsorted:
            warnings.warn("`as_2d` and `as_unsorted` can't be True at the same time, so as_unsorted will be ignored")
            as_unsorted = False
        self.as_unsorted = as_unsorted

        if self.as_3d:
            if self.as_unsorted:
                self._scan_unsorted(iline, xline, ic)
            else:
                try:
                    self._scan3d(iline, xline, istep, xstep, xloc, yloc)
                except Exception as e:
                    raise RuntimeError(f"{str(e)}\n This SEG-Y file may be unsorted, you can pass `as_unsorted` to view it as unsorted file, but it may be slow") from e
                self.metainfo = self.segy.get_metaInfo()
            self._shape3d = (self.metainfo.sizeZ, self.metainfo.sizeY, self.metainfo.sizeX)
        else:
            self.metainfo = self.segy.get_metaInfo()

        self._shape2d = (self.segy.trace_count, self.metainfo.sizeX)
        self._eval_range()

        if self.metainfo.scalar == 0:
            self._scalar = 1
        elif self.metainfo.scalar < 1:
            self._scalar = -1 / self.metainfo.scalar
        else:
            self._scalar = self.metainfo.scalar

    def _scan_unsorted(self, iline=None, xline=None, ic=None):
        """scan the unsorted SEG-Y file and create the geometry"""
        if ic is not None:
            ic = np.array(ic)
            if iline is not None or xline is not None:
                warnings.warn("'ic' is provided, so 'iline' and 'xline' will be ignored.", UserWarning)
            if ic.shape[0] != self.trace_count:
                raise ValueError("The length of ic should be the same as the trace count")

        elif iline is None or xline is None:
            raise ValueError("When 'as_unsorted' is True, either 'iline' and 'xline' must both be set, or 'ic' must be provided.")

        if ic is None:
            ic = self.segy.get_trace_keys([iline, xline], [4]*2, 0, self.trace_count)

        i0 = ic[:, 0].min()
        ie = ic[:, 0].max()
        diff = np.diff(np.sort(ic[:, 0]))
        diff = diff[diff != 0]
        istep = diff.min()
        if (ie - i0) % istep != 0:
            raise RuntimeError("can not create geomtry (error when determine iline/istep)")

        x0 = ic[:, 1].min()
        xe = ic[:, 1].max()
        diff = np.diff(np.sort(ic[:, 1]))
        diff = diff[diff != 0]
        xstep = diff.min()
        if (xe - x0) % xstep != 0:
            raise RuntimeError("can not create geomtry (error when determine xline/xstep)")

        ni = int((ie - i0) // istep + 1)
        nx = int((xe - x0) // xstep + 1)
        self.segy.setSteps(istep, xstep)
        self.metainfo = self.segy.get_metaInfo()
        self.metainfo.sizeZ = ni
        self.metainfo.sizeY = nx
        self.metainfo.min_inline = i0
        self.metainfo.max_inline = ie
        self.metainfo.min_crossline = x0
        self.metainfo.max_crossline = xe
        self.keylocs = [iline, xline, istep, xstep, 181, 185] # TODO: How to handle it?
        self._shape3d = (ni, nx, self.metainfo.sizeX)
        self.update_geometry(ic)


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

    def _read3d(self, idx):
        num = idx[1::2].count(-1)
        if num == 0 and not self.as_unsorted:
            self._check_bound(*idx)
            # HACK: when the cube is small, is it slow compare with using 2d collect? 
            d = self.segy.read(*idx)
        else:
            if not self.is_create_geometry:
                raise RuntimeError("Need create the geometry first, please call `update_geometry` first")
            ib, ie, xb, xe, tb, te = idx
            x = ib if ie == -1 else np.arange(ib, ie)
            y = xb if xe == -1 else np.arange(xb, xe)
            # TODO: CHECK bound
            x, y = np.meshgrid(x, y, indexing='ij')
            shape = x.shape
            if te == -1:
                d = self.read_traces_with_index(np.c_[x.flatten(), y.flatten()]).reshape(*shape, -1)
                d = d[:, :, tb]
            else:
                d = self.read_traces_with_index(np.c_[x.flatten(), y.flatten()], tb, te).reshape(*shape, -1)
        
        d = np.squeeze(d)
        if d.ndim == 0:
            d = float(d)
        return d

    def _read2d(self, idx):
        tb = 0 if idx[-1] == -1 else idx[2]
        te = self._shape2d[1] if idx[-1] == -1 else idx[3]

        if idx[1] == -1:
            assert idx[0].ndim == 1, "indices must be 1D array"
            data = self.segy.collect(idx[0], tb, te)
        else: # the first dimension is continuous
            assert idx[0] < idx[1] and idx[0] >= 0 and idx[1] <= self.shape[0], "index out of range"
            data = self.segy.collect(idx[0], idx[1], tb, te)

        if idx[-1] == -1:
            data = data[..., idx[-2]]

        return data

    def __getitem__(self, slices) -> np.ndarray:
        idx = self._process_keys(slices)
        if self.as_3d:
            return self._read3d(idx)
        else:
            return self._read2d(idx)

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
            if len(key) > N:
                raise IndexError(f"Too many dimensions: expected at most {N}, got {len(key)}")

            num_ellipsis = key.count(Ellipsis)
            if num_ellipsis > 1:
                raise ValueError("Only one ellipsis (...) allowed")
            if num_ellipsis == 1:
                idx = key.index(Ellipsis)
                n_insert = N - len(key) + 1
                key = key[:idx] + (slice(None),) * n_insert + key[idx + 1:]

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
                    start_idx[i] = np.array(k)
                    end_idx[i] = -1
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

        elif isinstance(key, (slice, List, np.ndarray)):
            if isinstance(key, slice):
                ib = 0 if key.start is None else key.start
                ie = self.shape[0] if key.stop is None else key.stop
            else:
                ib, ie = np.array(key), -1
            if self.as_3d:
                return ib, ie, 0, self.shape[1], 0, self.shape[2]
            else:
                return ib, ie, 0, self.shape[1]
        else:
            raise IndexError("Invalid index slices")

    def _check_bound(self, ib, ie, xb, xe, tb=None, te=None) -> None:
        assert ib < ie and ie <= self.shape[0] and ib >= 0, "index out of range"
        assert xb < xe and xe <= self.shape[1] and xb >= 0, "index out of range"
        if self.as_3d:
            assert (tb is not None) and (te is not None)
            assert tb < te and te <= self.shape[2] and tb >= 0, "index out of range"

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
        return self.to_numpy()

    def to_numpy(self):
        """like pandas"""
        if self.as_3d:
            if self.as_unsorted:
                return self[...]
            else:
                return self.segy.read()
        else:
            return self.segy.collect()

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
            raise RuntimeError("Need scan the file first, please call `to_3d` first")
        return self.segy.get_lineInfo()

    def update_trans_matrix(self, xyic=None):
        """
        Update the transformation matrix for inline/crossline to x/y

        parameters
        ----------
        xyic : np.ndarray, optional
            The inline/crossline and x/y information, shape is (n, 4), default None
        """
        if xyic is not None:
            assert xyic.shape[1] == 4, "The shape of xyic should be (n, 4)"
        elif self._shape3d is not None: # using lineinfo to get the xyic
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
            raise RuntimeError("Need scan the file first, please call `to_3d` first")

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

    @property
    def is_create_geometry(self):
        return self._geomety is not None

    def update_geometry(self, ic=None):
        if self._shape3d is None:
            raise RuntimeError("Need scan the file first, please call `to_3d` first")

        if ic is not None and ic.shape[0] != self.trace_count:
            raise ValueError("The length of ic should be the same as the trace count")

        if ic is None:
            ic = self.segy.get_trace_keys(self.keylocs[:2], [2]*2, 0, self.trace_count)
        if ic.max() >= max(self._shape3d[:2]): # move ic to zero-origin
            ic[:, 0] = (ic[:, 0] - self.iline_range[0]) / self.keylocs[2]
            ic[:, 1] = (ic[:, 1] - self.xline_range[0]) / self.keylocs[3]
            ic = np.round(ic).astype(np.int32)

        # check ic is in the range of (n1, n2)
        if ic[:, 0].min() < 0 or ic[:, 0].max() >= self._shape3d[0]:
            raise IndexError("the first dimension out of range")
        if ic[:, 1].min() < 0 or ic[:, 1].max() >= self._shape3d[1]:
            raise IndexError("the second dimension out of range")

        self._geomety = np.full(self.shape[:2], -1, dtype=np.int32) # (ni, nx), -1 means no trace
        self._geomety[ic[:, 0], ic[:, 1]] = np.arange(self.trace_count)

    def read_traces_with_index(self, index, tbeg=-1, tend=0):
        index = np.array(index).astype(np.int32)

        if index.ndim == 2 and index.shape[1] == 2:
            if not self.as_3d:
                raise RuntimeError("The SEG-Y file is treat as 2D array, index must be 1D array. Call `to_3d()` to view as a 3D")
            if self._geomety is None:
                raise RuntimeError("geometry is not created, Call `update_geometry` first")

            index = self._geomety[index[:, 0], index[:, 1]]

        assert index.ndim == 1, "index must be 1D array or list"
        return self.segy.collect(index, tbeg, tend)

    def map_to_index(self, index):
        if self._geomety is None:
            raise RuntimeError("geometry is not created, Call `update_geometry` first")

        index = np.array(index)
        assert index.ndim == 2 and index.shape[1] == 2, "index shape must be (N, 2)"
        index = self._geomety[index[:, 0], index[:, 1]]
        return index

    def arbitrary_line(self, points, ptype='auto', return_path=True, di=1):
        """
        Extract an arbitrary line from the the SEG-Y file. 
        The path is consisted by points

        Parameters
        -----------
        points : ArrayLike
            shape is (N, 2), it also can be a list
        ptype : str
            one of ['auto', 'zero', 'line', 'xy'], default is 'auto'.
            'zero' means points are taken from a geometry with zero-origin (i.e., numpy array indexes), 
            'line' means points are taken from the inline/xline geometry, 
            'xy' means points are taken from the X-Y geometry.
        return_path : bool
            if true, will return the path
        di : float
            the interval between two points of the path
        """
        if not self.as_3d:
            raise RuntimeError("`arbitrary_line` only support 3D SEG-Y file. Mybe you could call `to_3d` if is a 3D file")

        points = np.array(points)
        # HACK: Need optimize
        if ptype == 'auto':
            if points.max() > 100000:
                ptype = 'xy'
            elif points[:, 0].max() >= self.shape[0] or  points[:, 1].max() >= self.shape[1]:
                ptype = 'line'
            else:
                ptype = 'zero'

        if ptype == 'line':
            points[:, 0] -= self.iline_range[0]
            points[:, 1] -= self.xline_range[0]
        elif ptype == 'xy':
            points = self.xy_to_ix(points)

        if points[:, 0].max() >= self.shape[0] or points[:, 0].min() < 0 or  points[:, 1].max() >= self.shape[1] or points[:, 1].min() < 0:
            raise RuntimeError("`points` out range of the geometry. Maybe it's because the `ptype` is not set correctly (don't set 'auto')")

        out, p = arbitray_line(self, points, di)
        if not return_path:
            return out
        return out, p
