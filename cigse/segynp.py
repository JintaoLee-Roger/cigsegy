# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import List, Tuple
import numpy as np
from cigse.cpp import _CXX_SEGY
from cigse.transform import get_transform_metrix, apply_transform
from cigse.interp import arbitray_line
from cigse import utils
import warnings


class ScanMixin:

    def _scan_unsorted(self, iline=None, xline=None, ic=None):
        """scan the unsorted SEG-Y file and create the geometry"""
        if ic is not None:
            ic = np.array(ic)
            if iline is not None or xline is not None:
                warnings.warn(
                    "'ic' is provided, so 'iline' and 'xline' will be ignored.",
                    UserWarning)
            if ic.shape[0] != self.trace_count:
                raise ValueError(
                    "The length of ic should be the same as the trace count")

        elif iline is None or xline is None:
            raise ValueError(
                "When 'as_unsorted' is True, either 'iline' and 'xline' must both be set, or 'ic' must be provided."
            )

        if ic is None:
            ic = self.segy.get_trace_keys([iline, xline], [4] * 2, 0,
                                          self.trace_count)

        i0 = ic[:, 0].min()
        ie = ic[:, 0].max()
        diff = np.diff(np.sort(ic[:, 0]))
        diff = diff[diff != 0]
        istep = diff.min()
        if (ie - i0) % istep != 0:
            raise RuntimeError(
                "can not create geomtry (error when determine iline/istep)")

        x0 = ic[:, 1].min()
        xe = ic[:, 1].max()
        diff = np.diff(np.sort(ic[:, 1]))
        diff = diff[diff != 0]
        xstep = diff.min()
        if (xe - x0) % xstep != 0:
            raise RuntimeError(
                "can not create geomtry (error when determine xline/xstep)")

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
        self.keylocs = [iline, xline, istep, xstep, 181,
                        185]  # TODO: How to handle it?
        self._shape3d = (ni, nx, self.metainfo.sizeX)
        self.update_geometry(ic)

    def _scan3d(self,
                iline=None,
                xline=None,
                istep=None,
                xstep=None,
                xloc=None,
                yloc=None):
        [iline, xline, istep, xstep, xloc,
         yloc] = utils.guess(self.fname, iline, xline, istep, xstep, xloc,
                             yloc)[0]
        self.keylocs = [iline, xline, istep, xstep, xloc, yloc]
        self.segy.setInlineLocation(iline)
        self.segy.setCrosslineLocation(xline)
        self.segy.setSteps(istep, xstep)
        self.segy.scan()

        if self.metainfo.scalar == 0:
            self._scalar = 1
        elif self.metainfo.scalar < 1:
            self._scalar = -1 / self.metainfo.scalar
        else:
            self._scalar = self.metainfo.scalar

        # for coordinates transform
        self._trans_matrix = None
        self._geomety = None

    def _scan3d(self,
                iline=None,
                xline=None,
                istep=None,
                xstep=None,
                xloc=None,
                yloc=None):
        [iline, xline, istep, xstep, xloc,
         yloc] = utils.guess(self.fname, iline, xline, istep, xstep, xloc,
                             yloc)[0]
        self.keylocs = [iline, xline, istep, xstep, xloc, yloc]
        self.segy.setInlineLocation(iline)
        self.segy.setCrosslineLocation(xline)
        self.segy.setSteps(istep, xstep)
        self.segy.scan()

    def _eval_range(self):
        if self.trace_count < 6000:
            p0 = self.segy.collect(0, self.trace_count)
            self._min = mi
            self._max = ma
            return

        s, e = self.trace_count - 4000, 4000
        p0 = self.segy.collect(0, e)
        mi, ma = p0.min(), p0.max()

        p0 = self.segy.collect(s, self.trace_count)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        p0 = self.segy.collect(self.trace_count // 3,
                               self.trace_count // 3 + 4000)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        self._min = mi
        self._max = ma


class GeometryMixin:

    def map_to_index(self, index):
        if self._geomety is None:
            raise RuntimeError(
                "geometry is not created, Call `update_geometry` first")

        index = np.array(index)
        assert index.ndim == 2 and index.shape[
            1] == 2, "index shape must be (N, 2)"
        index = self._geomety[index[:, 0], index[:, 1]]
        return index

    def update_geometry(self, ic=None):
        if self._shape3d is None:
            raise RuntimeError(
                "Need scan the file first, please call `to_3d` first")

        if ic is not None and ic.shape[0] != self.trace_count:
            raise ValueError(
                "The length of ic should be the same as the trace count")

        if ic is None:
            ic = self.segy.get_trace_keys(self.keylocs[:2], [2] * 2, 0,
                                          self.trace_count)
        if ic.max() >= max(self._shape3d[:2]):  # move ic to zero-origin
            ic[:, 0] = (ic[:, 0] - self.iline_range[0]) / self.keylocs[2]
            ic[:, 1] = (ic[:, 1] - self.xline_range[0]) / self.keylocs[3]
            ic = np.round(ic).astype(np.int32)

        # check ic is in the range of (n1, n2)
        if ic[:, 0].min() < 0 or ic[:, 0].max() >= self._shape3d[0]:
            raise IndexError("the first dimension out of range")
        if ic[:, 1].min() < 0 or ic[:, 1].max() >= self._shape3d[1]:
            raise IndexError("the second dimension out of range")

        self._geomety = np.full(self.shape[:2], -1,
                                dtype=np.int32)  # (ni, nx), -1 means no trace
        self._geomety[ic[:, 0], ic[:, 1]] = np.arange(self.trace_count)

    def lineinfo(self):
        """
        Get the line information of the SEG-Y file,
        each line represents: inline, crossline_start, crossline_end, trace_start, trace_end, count
        """
        if self._shape3d is None:
            raise RuntimeError(
                "Need scan the file first, please call `to_3d` first")
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
        elif self._shape3d is not None:  # using lineinfo to get the xyic
            line = self.lineinfo()
            ni = line.shape[0]
            nx = self._shape3d[1]
            xyic = np.zeros((ni * 2 + nx * 2, 4))
            keylocs = [
                self.keylocs[4], self.keylocs[5], self.keylocs[0],
                self.keylocs[1]
            ]
            for i in range(ni):
                xyic[i * 2, :] = self.segy.get_trace_keys(
                    keylocs, [4] * 4, line[i, 3], line[i, 3] + 1)
                xyic[i * 2 + 1, :] = self.segy.get_trace_keys(
                    keylocs, [4] * 4, line[i, 4], line[i, 4] + 1)
            xyic[ni * 2:ni * 2 + nx, :] = self.segy.get_trace_keys(
                keylocs, [4] * 4, 0, nx)
            xyic[ni * 2 + nx:ni * 2 + nx * 2, :] = self.segy.get_trace_keys(
                keylocs, [4] * 4, self.trace_count - nx, self.trace_count)
        else:
            raise RuntimeError(
                "Need scan the file first, please call `to_3d` first")

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
        return np.round(apply_transform(ix, self._trans_matrix),
                        2).reshape(shape)

    @property
    def north(self):
        if self._north is None:
            xy1 = self.segy.get_trace_keys(self.keylocs[4:], [4]*2, int(self.trace_count//3), int(self.trace_count//3)+1) # yapf: disable
            xy1 = xy1.flatten()
            xy2 = xy1.copy()
            xy2[1] += 1000
            n = self.xy_to_ix([xy1, xy2])
            di = n[1] - n[0]
            self._north = di / np.linalg.norm(di)

        return self._north


class InterpMixin:

    def _check_points_valid(self, points, ptype='auto'):
        if not self.as_3d:
            raise RuntimeError(
                "`arbitrary_line` only support 3D SEG-Y file. Mybe you could call `to_3d` if is a 3D file"
            )

        points = np.array(points)
        # HACK: Need optimize
        if ptype == 'auto':
            if points.max() > 100000:
                ptype = 'xy'
            elif points[:, 0].max() >= self.shape[0] or points[:, 1].max(
            ) >= self.shape[1]:
                ptype = 'line'
            else:
                ptype = 'zero'

        if ptype == 'line':
            points[:, 0] -= self.iline_range[0]
            points[:, 1] -= self.xline_range[0]
        elif ptype == 'xy':
            points = self.xy_to_ix(points)

        if points[:, 0].max() >= self.shape[0] or points[:, 0].min(
        ) < 0 or points[:, 1].max() >= self.shape[1] or points[:, 1].min() < 0:
            raise RuntimeError(
                "`points` out range of the geometry. Maybe it's because the `ptype` is not set correctly (don't set 'auto')"
            )

        return points

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
        points = self._check_points_valid(points, ptype)

        out, p = arbitray_line(self, points, di)
        if not return_path:
            return out
        return out, p

    def interp(self, points, ptype='auto', method='linear'):
        """
        interpolate traces at points

        Parameters
        -----------
        points : ArrayLike
            shape is (..., 2)
        ptype : str
        ptype : str
            one of ['auto', 'zero', 'line', 'xy'], default is 'auto'.
            'zero' means points are taken from a geometry with zero-origin (i.e., numpy array indexes), 
            'line' means points are taken from the inline/xline geometry, 
            'xy' means points are taken from the X-Y geometry.
        method : str
            interpolation method, NOTE: Only support linear method now
        """
        points = self._check_points_valid(points, ptype)


class RWMixin:

    def _read3d(self, idx):
        num = idx[1::2].count(-1)
        if num == 0 and not self.as_unsorted:
            self._check_bound(*idx)
            # HACK: when the cube is small, is it slow compare with using 2d collect?
            d = self.segy.read(*idx)
        else:
            if not self.is_create_geometry:
                raise RuntimeError(
                    "Need create the geometry first, please call `update_geometry` first"
                )
            ib, ie, xb, xe, tb, te = idx
            x = ib if ie == -1 else np.arange(ib, ie)
            y = xb if xe == -1 else np.arange(xb, xe)
            # TODO: CHECK bound
            x, y = np.meshgrid(x, y, indexing='ij')
            shape = x.shape
            if te == -1:
                d = self.read_traces_with_index(np.c_[x.flatten(),
                                                      y.flatten()]).reshape(
                                                          *shape, -1)
                d = d[:, :, tb]
            else:
                d = self.read_traces_with_index(
                    np.c_[x.flatten(), y.flatten()], tb,
                    te).reshape(*shape, -1)

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
        else:  # the first dimension is continuous
            assert idx[0] < idx[1] and idx[0] >= 0 and idx[1] <= self.shape[
                0], "index out of range"
            data = self.segy.collect(idx[0], idx[1], tb, te)

        if idx[-1] == -1:
            data = data[..., idx[-2]]

        return data

    def read_traces_with_index(self, index, tbeg=-1, tend=0):
        index = np.array(index).astype(np.int32)

        if index.ndim == 2 and index.shape[1] == 2:
            if not self.as_3d:
                raise RuntimeError(
                    "The SEG-Y file is treat as 2D array, index must be 1D array. Call `to_3d()` to view as a 3D"
                )
            if self._geomety is None:
                raise RuntimeError(
                    "geometry is not created, Call `update_geometry` first")

            index = self._geomety[index[:, 0], index[:, 1]]

        assert index.ndim == 1, "index must be 1D array or list"
        return self.segy.collect(index, tbeg, tend)

    def to_numpy(self):
        """like pandas"""
        if self.as_3d:
            if self.as_unsorted:
                return self[...]
            else:
                return self.segy.read()
        else:
            return self.segy.collect()

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
                raise IndexError(
                    f"Too many dimensions: expected at most {N}, got {len(key)}"
                )

            num_ellipsis = key.count(Ellipsis)
            if num_ellipsis > 1:
                raise ValueError("Only one ellipsis (...) allowed")
            if num_ellipsis == 1:
                idx = key.index(Ellipsis)
                n_insert = N - len(key) + 1
                key = key[:idx] + (slice(None), ) * n_insert + key[idx + 1:]

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
                        raise IndexError(
                            f"only support step is 1, while got a step {k.step} in the {i}th dimension"
                        )

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
        elif isinstance(key, (List, np.ndarray)):
            raise NotImplementedError("Not implemented yet: TODO:")
        else:
            raise IndexError("Invalid index slices")

    def _check_bound(self, ib, ie, xb, xe, tb=None, te=None) -> None:
        assert ib < ie and ie <= self.shape[0] and ib >= 0, "index out of range"
        assert xb < xe and xe <= self.shape[1] and xb >= 0, "index out of range"
        if self.as_3d:
            assert (tb is not None) and (te is not None)
            assert tb < te and te <= self.shape[
                2] and tb >= 0, "index out of range"

    def tofile(self, fpath: str, load: bool = True):
        """
        save the SEG-Y file to a binary file without headers

        Parameters
        -----------
        fpath : str
            the save path
        load : bool
            if load is true, will load the file into memery first, then write to a file
        """
        if load:
            self.to_numpy().tofile(fpath)
        else:
            self.segy.tofile(fpath, not self.as_3d)


class InnerMixin:

    def __array__(self):
        """To support np.array(SegyNP(xxx))"""
        return self.to_numpy()

    def __getitem__(self, slices) -> np.ndarray:
        idx = self._process_keys(slices)
        if self.as_3d:
            return self._read3d(idx)
        else:
            return self._read2d(idx)

    def __array_function__(self, func, types, args, kwargs):
        if func is np.min:
            return self.min()
        elif func is np.max:
            return self.max()
        elif func is np.nanmin:
            return self.min()
        elif func is np.nanmax:
            return self.max()
        elif func is np.save:
            fpath, obj = args
            np.save(fpath, self.__array__())
            return
        raise NotImplementedError(
            f"Function {func} is not implemented for SegyNP")

    def __del__(self):
        self.close()

    def __repr__(self) -> str:
        out = f"cigsey.SegyNP class, file name: '{self.fname}'\n"
        out += f"shape: {self.shape}, nt is {self.shape[-1]}. Total traces: {self.trace_count}, interval: {self.interval}"
        return out

    def __len__(self) -> int:
        return self.shape[0]


class SegyNP:

    def __init__(self,
                 filename,
                 iline=None,
                 xline=None,
                 offset=None,
                 istep=None,
                 xstep=None,
                 ostep=None,
                 xloc=None,
                 yloc=None,
                 ndim=None,
                 as_unsorted=False,
                 keys=None) -> None:
        np.set_printoptions(suppress=True)
        self._ndim = ndim
        self._fname = filename
        self.segy = _CXX_SEGY.Pysegy(str(filename))
        self.keylocs = None

        self._shape2d = None
        self._shape3d = None
        self.as_unsorted = as_unsorted

        # for coordinates transform
        self._trans_matrix = None
        self._geomety = None
        self._north = None
