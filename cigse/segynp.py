# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import List, Tuple
import numpy as np
from cigse.cpp import _CXX_SEGY
from cigse.transform import get_transform_metrix, apply_transform
from cigse.interp import arbitray_line
from cigse import utils, plot
import warnings


class ScanMixin:
    _segy: _CXX_SEGY.Pysegy

    def _eval_range(self):
        if self.ntrace < 6000:
            p0 = self._segy.collect(0, self.ntrace, 0, self.nt)
            self._min = mi
            self._max = ma
            return

        s, e = self.ntrace - 4000, 4000
        p0 = self._segy.collect(0, e, 0, self.nt)
        mi, ma = p0.min(), p0.max()

        p0 = self._segy.collect(s, self.ntrace, 0, self.nt)
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        p0 = self._segy.collect(
            self.ntrace // 3,
            self.ntrace // 3 + 4000,
            0,
            self.nt,
        )
        mi = min(mi, p0.min())
        ma = max(ma, p0.max())

        self._min = mi
        self._max = ma


class GeometryMixin:
    _segy: _CXX_SEGY.Pysegy

    def map_to_index(self, index):
        pass

    def update_geometry(self, ic=None):
        pass

    def update_trans_matrix(self, xyic=None):
        pass

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


class PlotMixin:
    _segy: _CXX_SEGY.Pysegy

    def plot_region(self, mode='line'):
        """
        plot the region map (x and y axis are inline and crossline)

        Parameters
        -----------
        mode : str
            one of ['line', 'cdpxy', 'xy'], default is 'line'
        """
        plot.plot_region(self._segy, mode)

    def plot_trace_keys(self,
                        keyloc: int,
                        beg: int = 0,
                        end: int = 1000) -> None:
        """
        plot the values (at keyloc in each trace) of the traces 
        range from beg to end .
        """
        plot.plot_trace_keys(self._segy, keyloc, beg, end)

    def plot_trace_2(
        self,
        beg: int = 0,
        end: int = 1000,
        k1: int = None,
        k2: int = None,
    ):
        """
        plot the values of the traces at k1 and k2, range from beg to end.
        If k1 and k2 is None, will plot iline/xline
        """
        if k1 is None:
            k1 = self._keyloc['iline']
        if k2 is None:
            k2 = self._keyloc['xline']
        plot.plot_trace_ix(self._segy, k1, k2, beg, end)

    def plot_trace_3(
        self,
        beg: int = 0,
        end: int = 1000,
        k1: int = None,
        k2: int = None,
        k3: int = None,
    ):
        """
        plot the values of the traces at k1, k2 and k3, range from beg to end.
        If k1, k2 and k3 is None, will plot iline/xline/offset
        """
        if k1 is None:
            k1 = self._keyloc['iline']
        if k2 is None:
            k2 = self._keyloc['xline']
        if k3 is None:
            k3 = self._keyloc['offset']
        plot.plot_trace_ixo(self._segy, k1, k2, k3, beg, end)

    def plot3d(self):
        """
        plot 3d
        """
        raise NotImplementedError("not implement yet")


class InterpMixin:
    _segy: _CXX_SEGY.Pysegy

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

    def extract_arbitrary_line_by_view(self,
                                       bmap: str = 'data',
                                       draw_arb: bool = True,
                                       *,
                                       line: bool = True,
                                       idx: int = 50,
                                       cline='#F3AA3C'):
        """
        extract arbitrary line from seismic data by clicking

        Parameters
        ----------
        - bmap: str
            background map, 'data' or 'blank'
        - line : bool
            whether to draw the broken line 
        - idx: int 
            the slice index of the seismic data if bmap is 'data'
        - cline: str
            color of the line

        Returns
        -------
        - out: np.ndarray
            extracted arbitrary line
        - p: np.ndarray
            extracted arbitrary line path
        - coords: np.ndarray
            the coordinates by clicking
        """

    def align_coordinates(self, fname: str):
        """
        Interpolate the input SEG-Y data to align with the self coordinate system.
        """


class RWMixin:
    _segy: _CXX_SEGY.Pysegy

    def to_numpy(self):
        """like pandas"""
        if self.unsorted and self.ndim > 2:
            return self[...]
        else:
            return self._segy.read()

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
            self._segy.tofile(fpath, self.ndim == 2)

    # fmt: off
    def _read4d(self, idx) -> np.ndarray:
        assert self.ndim == 4, "The data is not 4D"
        self._check_bound2(idx)

        num = idx[1::2].count(-1)

        if num == 0 and not self.unsorted:
            d = self._segy.read4d(*idx)
        else:
            if not self.is_create_geometry:
                raise RuntimeError("Need create the geometry first, please call `update_geometry` first")
            grid, shape = self._create_meshgrid(idx[:-2])
            tidx = self.map_to_index(grid)
            if idx[-1] is None:
                d = self._segy.collect(tidx, 0, self.nt).reshape(*shape, -1)
                d = d[..., idx[-2]]
            else:
                d = self._segy.collect(tidx, idx[-2], idx[-1]).reshape(*shape, -1)

        return self._post_process(d)


    def _read3d(self, idx) -> np.ndarray:
        assert self.ndim == 3, "The data is not 3D"
        self._check_bound2(idx)
        num = idx[1::2].count(-1)
        if num == 0 and not self.unsorted:
            d = self._segy.read3d(*idx)
        else:
            if not self.is_create_geometry:
                raise RuntimeError("Need create the geometry first, please call `update_geometry` first")
            grid, shape = self._create_meshgrid(idx[:-2])
            tidx = self.map_to_index(grid)
            if idx[-1] is None:
                d = self._segy.collect(tidx, 0, self.nt).reshape(*shape, -1)
                d = d[..., idx[-2]]
            else:
                d = self._segy.collect(tidx, idx[-2], idx[-1]).reshape(*shape, -1)

        return self._post_process(d)


    def _read2d(self, idx) -> np.ndarray:
        assert self.ndim == 2, "The data is not 2D"
        self._check_bound2(idx)

        # the first dim is ndarray
        if idx[1] is None:
            # the time dim is ndarray
            if idx[-1] is None:
                data = self._segy.collect(idx[0], 0, self.nt)
                data = data[:, idx[2]]
            else:
                data = self._segy.collect(idx[0], idx[2], idx[3])
        else:  # the first dim is ib, ie
            if idx[-1] is None:
                data = self._segy.collect(idx[0], idx[1], 0, self.nt)
                data = data[:, idx[2]]
            else:
                data = self._segy.collect(idx[0], idx[1], idx[2], idx[3])

        return self._post_process(data)


    def _write4d(self, idx, data) -> None:
        assert self.ndim == 4, "The data is not 4D"
        self._check_bound2(idx)
        self._check_wmode_lastdim(idx)
        self._check_data_shape(data, idx)

        num = idx[1::2].count(-1)

        if num == 0 and not self.unsorted:
            self._segy.write4d(data, *idx)
        else:
            if not self.is_create_geometry:
                raise RuntimeError("Need create the geometry first, please call `update_geometry` first")
            grid, shape = self._create_meshgrid(idx[:-2])
            tidx = self.map_to_index(grid)
            self._segy.write_traces(data, tidx, idx[-2], idx[-1])


    def _write3d(self, idx, data) -> None:
        assert self.ndim == 3, "The data is not 3D"
        self._check_bound2(idx)
        self._check_wmode_lastdim(idx)
        self._check_data_shape(data, idx)

        num = idx[1::2].count(-1)

        if num == 0 and not self.unsorted:
            self._segy.write3d(data, *idx)
        else:
            if not self.is_create_geometry:
                raise RuntimeError("Need create the geometry first, please call `update_geometry` first")
            grid, shape = self._create_meshgrid(idx[:-2])
            tidx = self.map_to_index(grid)
            self._segy.write_traces(data, tidx, idx[-2], idx[-1])


    def _write2d(self, idx, data) -> None:
        assert self.ndim == 2, "The data is not 2D"
        assert data.ndim <= 2, "The data is not 2D"
        self._check_bound2(idx)
        self._check_wmode_lastdim(idx)
        self._check_data_shape(data, idx)

        # the first dim is ndarray
        if idx[1] is None:
            # the time dim is ndarray
            self._segy.write_traces(data, idx[0], idx[2], idx[3])
        else:  # the first dim is ib, ie

            self._segy.write_traces(idx[0], idx[1], idx[2], idx[3])


    def _create_meshgrid(self, idx):
        if idx[1] is None:
            x = idx[0]
        else:
            x = np.arange(idx[0], idx[1])
        if idx[3] is None:
            y = idx[2]
        else:
            y = np.arange(idx[2], idx[3])

        if self.ndim == 3:
            X, Y = np.meshgrid(x, y, indexing='ij')
            shape = X.shape
            return np.c_[X.flatten(), Y.flatten()], shape
        elif self.ndim == 4:
            if idx[5] is None:
                z = idx[4]
            else:
                z = np.arange(idx[4], idx[5])
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            shape = X.shape
            return np.c_[X.flatten(), Y.flatten(), Z.flatten()], shape
        else:
            raise RuntimeError("only support ndim == 3 or ndim == 4")


    def _post_process(self, d: np.ndarray):
        d = d.squeeze()
        if d.ndim == 1 and d.size == 1:
            return d[0]
        elif d.ndim == 0:
            return float(d)
        return d


    def _process_keys(self, key) -> List:
        if isinstance(key, (int, np.integer)):

            if key < 0:
                key += self.shape[0]
            if key < 0 or key >= self.shape[0]:
                raise IndexError("Index out of range")

            if self.ndim == 2:
                return key, key + 1, 0, self.shape[1]
            elif self.ndim == 3:
                return key, key + 1, 0, self.shape[1], 0, self.shape[2]
            else:
                return key, key + 1, 0, self.shape[1], 0, self.shape[2], 0, self.shape[3]

        elif key is Ellipsis:

            if self.ndim == 2:
                return 0, self.shape[0], 0, self.shape[1]
            elif self.ndim == 3:
                return 0, self.shape[0], 0, self.shape[1], 0, self.shape[2]
            else:
                return 0, self.shape[0], 0, self.shape[1], 0, self.shape[2], 0, self.shape[3]

        elif isinstance(key, Tuple):

            return self._process_keys_tuple(key)

        elif isinstance(key, (slice, List, np.ndarray)):

            if isinstance(key, slice):
                ib = 0 if key.start is None else key.start
                ie = self.shape[0] if key.stop is None else key.stop
                assert key.step is None, "step is not supported"
            else:
                ib, ie = np.array(key), None

            if self.ndim == 2:
                return ib, ie, 0, self.shape[1]
            elif self.ndim == 3:
                return ib, ie, 0, self.shape[1], 0, self.shape[2]
            else:
                return ib, ie, 0, self.shape[1], 0, self.shape[2], 0, self.shape[3]

        else:
            raise IndexError("Invalid index slices")


    def _process_keys_tuple(self, key) -> List:
        if len(key) > self.ndim:
            raise IndexError(f"Too many dimensions: expected at most {self.ndim}, got {len(key)}")

        num_ellipsis = key.count(Ellipsis)
        if num_ellipsis > 1:
            raise ValueError("Only one ellipsis (...) allowed")

        if num_ellipsis == 1:
            idx = key.index(Ellipsis)
            n_insert = self.ndim - len(key) + 1
            key = key[:idx] + (slice(None),) * n_insert + key[idx + 1:]

        start_idx = [None] * self.ndim
        end_idx = [None] * self.ndim
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

        out = []
        for i in range(self.ndim):
            if start_idx[i] is None:
                out.append(0)
            else:
                out.append(start_idx[i])

            if end_idx[i] is None:
                out.append(self.shape[i])
            else:
                out.append(end_idx[i])

        return out



class CheckMixin:

    def _check_bound(self, dim, ib, ie):
        if ie == None:
            assert isinstance(ib, np.ndarray) and ib.ndim == 1, f"if ie is 0, ib must be a 1D numpy array"
            assert ib.min() >= 0 and ib.max() <= self.shape[dim], f"index array out of range in dim {dim}"
        else:
            assert ib >= 0 and ib < ie and ib <= self.shape[dim], f"index out of range in dim {dim}"


    def _check_bound2(self, idx):
        assert len(idx) == self.ndim * 2, f"ndim is {self.ndim}, need {self.ndim*2} idx, but got {len(idx)}"
        for i in range(self.ndim):
            self._check_bound(i, idx[i*2], idx[i*2+1])

    def _check_wmode(self):
        if not self._wmode:
            raise RuntimeError("The SEG-Y file is not writable, as you set the `mode` to 'r'. If you want to enable write mode, set to `rw`")

    def _check_wmode_lastdim(self, idx):
        if idx[-1] is None:
            raise TypeError("Indexing with a list or ndarray for the last dimension is not supported when writing.")

    def _check_data_shape(self, data: np.ndarray, idx):
        dstshape = []
        for i in range(self.ndim):
            if idx[i * 2 + 1] is None:
                dstshape.append(idx[i * 2].shape)
            else:
                dstshape.append(idx[i * 2 + 1] - idx[i * 2])
        dstshape = [k for k in dstshape if k != 1]
        assert tuple(dstshape) == data.squeeze().shape, "shape of the input data is not match the shape of the data to write"
# fmt: on


class SegyCMixin:

    def cut(self):
        pass

    def create_by_sharing_header(self):
        pass


# TODO:
class SeGyAccessMixin:

    def __getattr__(self, name):
        if name in ['iline', 'xline', 'offset', 'coordx', 'coordy', 'itrace']:
            return self._SeGyAccessor(self._segy, name)
        return super().__getattribute__(name)

    def __setattr__(self, name: str, value) -> None:
        if name in ['iline', 'xline', 'offset', 'coordx', 'coordy', 'itrace']:
            return self._SeGyAccessor(self._segy, name)
        return super().__setattr__(name)

    class _SeGyAccessor:

        def __init__(self, segy, attribute):
            self._segy = segy
            self.attribute = attribute

        def __getitem__(self, index):
            method = getattr(self._segy, self.attribute)
            return method(index)

        def __setitem__(self, index, value):
            return


class InnerMixin:
    _segy: _CXX_SEGY.Pysegy

    def __array__(self):
        """To support np.array(SegyNP(xxx))"""
        return self.to_numpy()

    def __getitem__(self, slices) -> np.ndarray:
        idx = self._process_keys(slices)
        if self._ndim == 4:
            return self._read4d(idx)
        elif self._ndim == 3:
            return self._read3d(idx)
        else:
            return self._read2d(idx)

    def __setitem__(self, slices, data: np.ndarray) -> None:
        self._check_wmode()
        if data.dtype != np.float32:
            raise TypeError("The data type of the input data must be np.float32") # yapf: disable
        idx = self._process_keys(slices)
        if self._ndim == 4:
            return self._write4d(idx, data)
        elif self._ndim == 3:
            return self._write3d(idx, data)
        else:
            return self._write2d(idx, data)

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
        self._segy.close()

    def __repr__(self) -> str:
        out = f"cigsey.SegyNP class, file name: '{self.fname}'\n\n"
        keys = self._segy.get_keylocs()
        meta = self._segy.get_metainfo()
        meta = {**keys, **meta}
        meta = utils.post_process_meta(meta)
        out += utils.parse_metainfo(meta)
        return out

    def __len__(self) -> int:
        return self.shape[0]


class SegyNP(InnerMixin, RWMixin, InterpMixin, PlotMixin, GeometryMixin,
             ScanMixin, CheckMixin):

    def __init__(self,
                 filename: str,
                 keyloc: dict = None,
                 mode: str = 'r',
                 *,
                 ndim: int = None,
                 as_unsorted: bool = False,
                 keys: dict = None) -> None:
        np.set_printoptions(suppress=True)

        assert mode in ['r', 'rw'], "`mode` only can be 'r' or 'rw'"
        self._ndim = ndim
        self._fname = filename
        self._segy = _CXX_SEGY.Pysegy(str(filename))
        self._wmode = False if mode == 'r' else True
        if self._wmode:
            warnings.warn(
                "Dangerous!!! You are using a writable mode ('rw'), which may **alter** the SEG-Y file. "
                "It is strongly recommended to make a **backup copy** of the file before proceeding "
                "to avoid any potential irreversible changes.", UserWarning)

        self._shape2 = None
        self._shape3 = None
        self._unsorted = as_unsorted

        # for coordinates transform
        self._trans_matrix = None
        self._geomety = None
        self._north = None

        # for values
        self._min = None
        self._max = None

        self._keyloc = None
        self._metainfo = None
        self._lininfo = None

    @property
    def ntrace(self) -> int:
        """
        Number of traces in the SEG-Y file
        """
        return self._segy.ntrace

    @property
    def ndim(self) -> int:
        """
        Number of dimensions of the data
        """
        return self._ndim

    @property
    def shape(self) -> Tuple:
        """
        the shape of the data
        """
        if self.ndim == 2:
            return self._shape2
        else:
            return self._shape3

    @property
    def dtype(self) -> np.dtype:
        """
        the data type of the data, always be np.float32
        """
        return np.float32

    @property
    def nt(self) -> int:
        """
        length of the time axis, i.e., the number of samples for each trace
        """
        return self._segy.nt

    @property
    def file_name(self) -> str:
        """
        the file name of the SEG-Y file
        """
        return self._fname

    @property
    def north(self):
        """
        the north direction of the SEG-Y file, only available for 3D/4D data
        """
        if self._north is None:
            pass

        return self._north

    @property
    def unsorted(self) -> bool:
        """
        whether the SEG-Y file is unsorted
        """
        return self._unsorted

    @property
    def is_create_geometry(self) -> bool:
        """
        whether the geometry is created
        """
        return self._geometry is not None

    def close(self) -> None:
        """
        close the SEG-Y file
        """
        self._segy.close()

    def to_2d(self):
        """
        Treat the SEG-Y file as a collection of traces, shape is like (ntrace, nt)
        """
        self._ndim = 2

    def to_nd(self):
        """
        Treat the SEG-Y file as a 3D/4D array. If the SEG-Y file is scanned, the keylocs will be ignored 
        """
        # TODO:
        self._ndim = self._segy.ndim

    def max(self, real=False) -> float:
        """
        return the maximum value of the data, 
        if real is False, the maximum value is not real, and we just read a part of traces to calculate max
        if real is True, we read all traces to calculate max
        """
        if real:
            return self[...].min()
        return self._max

    def min(self, real=False) -> float:
        """
        return the min value of the data, 
        if real is False, the min value is not real, and we just read a part of traces to calculate min
        if real is True, we read all traces to calculate min
        """
        if real:
            return self[...].max()
        return self._min

    @property
    def lineinfo(self):
        return self._lininfo

    @property
    def metainfo(self):
        return self._metainfo

    @property
    def keylocs(self):
        return self._keyloc
