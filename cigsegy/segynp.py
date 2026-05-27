# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from dataclasses import dataclass
from typing import List, Optional, Tuple
from itertools import product
import numpy as np
from cigsegy.cpp import _CXX_SEGY
from cigsegy.transform import get_transform_metrix, apply_transform
from cigsegy.interp import arbitray_line
from cigsegy import utils, plot, tools, createtool
import warnings
from tqdm import tqdm


@dataclass
class AxisSelector:
    indices: np.ndarray
    is_scalar: bool = False

    def __post_init__(self):
        self.indices = np.asarray(self.indices, dtype=np.int64).reshape(-1)

    @property
    def size(self) -> int:
        return int(self.indices.size)

    def contiguous_bounds(self) -> Optional[Tuple[int, int]]:
        if self.size == 0:
            return (0, 0)
        if self.size == 1:
            idx = int(self.indices[0])
            return (idx, idx + 1)
        steps = np.diff(self.indices)
        if np.all(steps == 1):
            return int(self.indices[0]), int(self.indices[-1]) + 1
        return None

    def constant_step(self) -> Optional[int]:
        if self.size <= 1:
            return 1
        steps = np.diff(self.indices)
        if np.all(steps == steps[0]):
            return int(steps[0])
        return None

    def dense_window(self) -> Tuple[int, int, np.ndarray]:
        if self.size == 0:
            return 0, 0, np.empty(0, dtype=np.int64)
        lo = int(self.indices.min())
        hi = int(self.indices.max()) + 1
        return lo, hi, (self.indices - lo).astype(np.int64, copy=False)

    def full_axis_step(self, axis_size: int) -> Optional[int]:
        if self.is_scalar or self.size == 0:
            return None
        step = self.constant_step()
        if step is None or step <= 0:
            return None
        expected = np.arange(0, axis_size, step, dtype=np.int64)
        if np.array_equal(self.indices, expected):
            return step
        return None


@dataclass
class NormalizedIndex:
    selectors: Tuple[AxisSelector, ...]
    result_tokens: Tuple[Tuple[str, int], ...]

    @property
    def has_newaxis(self) -> bool:
        return any(kind == 'newaxis' for kind, _ in self.result_tokens)

    @property
    def base_shape(self) -> Tuple[int, ...]:
        return tuple(selector.size for selector in self.selectors)

    @property
    def result_shape(self) -> Tuple[int, ...]:
        if not self.result_tokens:
            return ()
        shape = []
        for kind, axis in self.result_tokens:
            if kind == 'newaxis':
                shape.append(1)
            else:
                shape.append(self.selectors[axis].size)
        return tuple(shape)

    def format_result(self, data):
        arr = np.asarray(data)
        if arr.shape != self.base_shape:
            arr = arr.reshape(self.base_shape)
        if not self.result_tokens:
            return arr.reshape(()).item()
        return arr.reshape(self.result_shape)

    def prepare_value(self, value, dtype=np.float32) -> np.ndarray:
        if self.has_newaxis:
            raise IndexError("np.newaxis is not supported in assignment")

        arr = np.asarray(value, dtype=dtype)
        if self.result_shape == ():
            if arr.size != 1:
                raise ValueError(
                    f"cannot assign input with shape {arr.shape} to a scalar selection"
                )
            arr = arr.reshape(())
        else:
            try:
                arr = np.broadcast_to(arr, self.result_shape)
            except ValueError as exc:
                raise ValueError(
                    f"cannot broadcast input shape {arr.shape} to indexed shape {self.result_shape}"
                ) from exc
            arr = np.asarray(arr, dtype=dtype)

        if arr.shape != self.base_shape:
            arr = arr.reshape(self.base_shape)
        return np.ascontiguousarray(arr, dtype=dtype)


class ScanMixin:
    _segy: _CXX_SEGY.Pysegy

    def _eval_range(self):
        if self.ntrace < 6000:
            p0 = self._segy.collect(0, self.ntrace, 0, self.nt)
            self._min = float(p0.min())
            self._max = float(p0.max())
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

    def _scan(self, keylocs=None):
        if keylocs is None:
            keylocs = utils.guess(self._segy)
        elif isinstance(keylocs, dict):
            keylocs = utils.guess(self._segy, **keylocs)
        else:
            keylocs = utils.guess(self._segy, *keylocs)
        self._segy.setLocations(*keylocs[:3])
        self._segy.setSteps(*keylocs[3:6])
        self._segy.setXYLocations(*keylocs[6:8])
        ndim = 4 if keylocs[-1] is True else 3
        self._segy.set_segy_type(ndim)
        self._segy.scan()
        self._keylocs = self._segy.get_keylocs()
        self._metainfo = self._segy.get_metainfo()
        self._metainfo = {**self._keylocs, **self._metainfo}
        utils.post_process_meta(self._segy, self._metainfo, False)
        self._ndim = self._segy.ndim
        self._shape3 = self._segy.shape
        # self._lineinfo = self._segy.get_lineInfo() # TODO:

    def _scan_unsorted(self, keylocs=None, keys=None, ndim=None):
        if keys is not None:
            geom = tools.full_scan(self._segy, keys=keys)
        elif keylocs is None:
            raise ValueError(
                "When 'as_unsorted' is True, either 'keylocs' must be set, or 'keys' must be provided."
            )
        else:
            if ndim is None:
                is4d = None
            else:
                is4d = ndim == 4
            if isinstance(keylocs, (list, tuple)):
                geom = tools.full_scan(self._segy, *keylocs[:3], is4d=is4d)
            else:
                offset = keylocs.get('offset', 37)
                geom = tools.full_scan(self._segy, keylocs['iline'],
                                       keylocs['xline'], offset)
        self.update_geometry(geom['geom'])
        self._parser_unsorted_infos(geom)

    def _parser_unsorted_infos(self, geom):
        self._keylocs = self._segy.get_keylocs()
        self._metainfo = self._segy.get_metainfo()

        self._shape3 = tuple(geom['shape'])
        self._ndim = len(self._shape3)

        loc = geom['location']
        self._keylocs = dict(iline=loc[0], xline=loc[1])
        self._keylocs['istep'] = geom['iline']['istep']
        self._keylocs['xstep'] = geom['xline']['xstep']
        if len(loc) == 3:
            self._keylocs['offset'] = loc[2]
            self._keylocs['ostep'] = geom['offset']['ostep']

        self._metainfo['ni'] = self._shape3[0]
        self._metainfo['start_iline'] = geom['iline']['min_iline']
        self._metainfo['end_iline'] = geom['iline']['max_iline']
        self._metainfo['nx'] = self._shape3[1]
        self._metainfo['start_xline'] = geom['xline']['min_xline']
        self._metainfo['end_xline'] = geom['xline']['max_xline']
        if len(self._shape3) == 4:
            self._metainfo['no'] = self._shape3[2]
            self._metainfo['start_offset'] = geom['offset']['min_offset']
            self._metainfo['end_offset'] = geom['offset']['max_offset']
        self._metainfo['ndim'] = self._ndim


class GeometryMixin:
    _segy: _CXX_SEGY.Pysegy

    def map_to_indices(self, indices):
        if self._geometry is None:
            raise RuntimeError("geometry is not created, Call `update_geometry` first") # yapf: disable
        if self.ndim == 2:
            raise RuntimeError("ndim is 2, unsupport this function") # yapf: disable

        indices = np.array(indices)
        assert indices.ndim == 2
        if self.ndim == 3:
            assert indices.shape[1] == 2
            indices = self._geometry[indices[:, 0], indices[:, 1]]
        elif self.ndim == 4:
            assert indices.shape[1] == 3
            indices = self._geometry[indices[:, 0], indices[:, 1], indices[:, 2]]
        return indices

    def update_geometry(self, geom=None):
        if geom is not None:
            self._geometry = geom
            return
        warnings.warn("This may be slow...")
        geom = tools.full_scan(self._segy, self._keylocs['iline'],
                               self._keylocs['xline'], self._keylocs['offset'])
        self._geometry = geom['geom']

    def update_trans_matrix(self, xyic=None):
        if self.ndim == 2:
            raise RuntimeError("ndim is 2, unsupport this function")
        if xyic is None:
            xyic = tools.get_lineInfo(self._segy, mode='geom')
            xyic = xyic[:, [2, 3, 0, 1]]
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
            ic[:, 0] -= self._metainfo['start_iline']
            ic[:, 1] -= self._metainfo['start_xline']
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
            ix[:, 0] += self._metainfo['start_iline']
            ix[:, 1] += self._metainfo['start_xline']
        if self._trans_matrix is None:
            self.update_trans_matrix()
        return np.round(apply_transform(ix, self._trans_matrix), 2).reshape(shape) # yapf: disable


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

    def plot_trace_keys2(
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
            k1 = self._keylocs['iline']
        if k2 is None:
            k2 = self._keylocs['xline']
        plot.plot_trace_ix(self._segy, k1, k2, beg, end)

    def plot_trace_keys3(
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
            k1 = self._keylocs['iline']
        if k2 is None:
            k2 = self._keylocs['xline']
        if k3 is None:
            k3 = self._keylocs['offset']
        plot.plot_trace_ixo(self._segy, k1, k2, k3, beg, end)

    def plot3d(self, use_viser=False):
        """
        plot 3d
        """
        assert self.ndim == 3, "The data is not 3D"
        if not use_viser:
            try:
                import cigvis
            except ImportError:
                raise ImportError("To use this function, you need to install cigvis")
            nodes = cigvis.create_slices(self)
            cigvis.plot3D(nodes)
        else:
            try:
                from cigvis import viserplot
            except ImportError:
                raise ImportError("To use this function, you need to install cigvis['viser']")
            nodes = viserplot.create_slices(self)
            viserplot.plot3D(nodes)


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
            'zero' means points are taken from a geometry with zero-origin (i.e., numpy array indices), 
            'line' means points are taken from the inline/xline geometry, 
            'xy' means points are taken from the X-Y geometry.
        return_path : bool
            if true, will return the path
        di : float
            the interval between two points of the path
        """
        points = self._process_points(points, ptype)
        out, p, indices = arbitray_line(self, points, di)
        if not return_path:
            return out
        return out, p, indices

    def extract_arbitrary_line_by_view(self,
                                       bmap: str = 'data',
                                       draw_arb: bool = True,
                                       *,
                                       return_values: bool = True,
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
        out = plot.extract_arbitrary_line_by_view(
            self,
            bmap,
            draw_arb,
            line=line,
            idx=idx,
            cline=cline,
        )
        if return_values:
            return out

    def align_coordinates(self, fname: str):
        """
        Interpolate the input SEG-Y data to align with the self coordinate system.
        """
        raise NotImplementedError("Not Implemented yet")

    def _process_points(self, points, ptype='auto'):
        assert self.ndim == 3, "The data is not 3D"
        points = np.array(points)
        # HACK: Need optimize
        if ptype == 'auto':
            if points.max() > 100000:
                ptype = 'xy'
            elif points[:, 0].max() >= self.shape[0] or points[:, 1].max() >= self.shape[1]: # yapf: disable
                ptype = 'line'
            else:
                ptype = 'zero'

        if ptype == 'line':
            points[:, 0] -= self.iline_range[0]
            points[:, 1] -= self.xline_range[0]
        elif ptype == 'xy':
            points = self.xy_to_ix(points)

        self._check_bound(0, points[:, 0].min(), points[:, 0].max())
        self._check_bound(1, points[:, 1].min(), points[:, 1].max())
        return points


class RWMixin:
    _segy: _CXX_SEGY.Pysegy

    def to_numpy(self):
        """like pandas"""
        self._segy.show_progress(True)
        if self._ignore or (self.unsorted and self.ndim > 2):
            d = self[...]
        else:
            d = self._segy.read()
            if self._T:
                d = d.T
        self._segy.show_progress(False)

        return d

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
        self._segy.show_progress(True)
        if load:
            self.to_numpy().tofile(fpath)
        else:
            if self._T:
                raise ValueError("Cannot save the data in .T mode and not load")
            self._segy.tofile(fpath, self.ndim == 2)
        self._segy.show_progress(False)

    def _collect_with_valid_indices(self, tidx, shape, tb: int, te: int):
        tidx = np.asarray(tidx, dtype=np.int32).reshape(-1)
        ns = te - tb
        out = np.zeros((tidx.size, ns), dtype=np.float32)
        valid = tidx >= 0
        if np.any(valid):
            out[valid] = self._segy.collect(tidx[valid], tb, te)
        return out.reshape(*shape, ns)

    def _require_writable_trace_indices(self, tidx):
        tidx = np.asarray(tidx, dtype=np.int32).reshape(-1)
        if np.any(tidx < 0):
            raise IndexError("Cannot write into missing traces in the geometry")
        return tidx

    def _data_shape(self) -> Tuple[int, ...]:
        if self._ndim == 2:
            return tuple(self._shape2)
        return tuple(self._shape3)

    def _normalize_axis_key(self, key, axis_size: int, axis: int) -> AxisSelector:
        if isinstance(key, (int, np.integer)):
            idx = int(key)
            if idx < 0:
                idx += axis_size
            if idx < 0 or idx >= axis_size:
                raise IndexError(
                    f"index {key} is out of bounds for axis {axis} with size {axis_size}"
                )
            return AxisSelector(np.array([idx]), is_scalar=True)

        if isinstance(key, slice):
            try:
                start, stop, step = key.indices(axis_size)
            except ValueError as exc:
                raise ValueError(str(exc)) from exc
            return AxisSelector(np.arange(start, stop, step, dtype=np.int64))

        if isinstance(key, (list, tuple, np.ndarray)):
            arr = np.asarray(key)
            if arr.ndim == 0:
                return self._normalize_axis_key(arr.item(), axis_size, axis)
            if arr.ndim != 1:
                raise IndexError(
                    f"only support 1D index arrays, while got ndim={arr.ndim} in axis {axis}"
                )
            if arr.dtype == np.bool_:
                if arr.size != axis_size:
                    raise IndexError(
                        f"boolean index size {arr.size} does not match axis {axis} size {axis_size}"
                    )
                return AxisSelector(np.flatnonzero(arr))
            if not np.issubdtype(arr.dtype, np.integer):
                raise IndexError(
                    f"index array for axis {axis} must be integer or boolean, got {arr.dtype}"
                )
            arr = arr.astype(np.int64, copy=False)
            arr = np.where(arr < 0, arr + axis_size, arr)
            if arr.size and ((arr < 0).any() or (arr >= axis_size).any()):
                raise IndexError(
                    f"index array out of bounds for axis {axis} with size {axis_size}"
                )
            return AxisSelector(arr)

        raise IndexError("Invalid index slices")

    def _normalize_key(self, key, shape=None) -> NormalizedIndex:
        shape = tuple(self.shape if shape is None else shape)
        ndim = len(shape)

        if not isinstance(key, tuple):
            key = (key, )
        tokens = list(key)

        num_ellipsis = sum(token is Ellipsis for token in tokens)
        if num_ellipsis > 1:
            raise IndexError("Only one ellipsis (...) allowed")

        if num_ellipsis == 1:
            consumed = sum(token is not Ellipsis and token is not None
                           for token in tokens)
            if consumed > ndim:
                raise IndexError(
                    f"Too many dimensions: expected at most {ndim}, got {consumed}"
                )
            fill = ndim - consumed
            ellipsis_idx = next(i for i, token in enumerate(tokens)
                                if token is Ellipsis)
            tokens = tokens[:ellipsis_idx] + [slice(None)] * fill + tokens[
                ellipsis_idx + 1:]

        consumed = sum(token is not None for token in tokens)
        if consumed > ndim:
            raise IndexError(
                f"Too many dimensions: expected at most {ndim}, got {consumed}"
            )

        tokens.extend([slice(None)] * (ndim - consumed))

        selectors = []
        result_tokens = []
        axis = 0
        for token in tokens:
            if token is None:
                result_tokens.append(('newaxis', -1))
                continue

            selector = self._normalize_axis_key(token, shape[axis], axis)
            selectors.append(selector)
            if not selector.is_scalar:
                result_tokens.append(('axis', len(selectors) - 1))
            axis += 1

        return NormalizedIndex(tuple(selectors), tuple(result_tokens))

    def _selectors_to_bounds(self, selectors) -> Optional[List[int]]:
        out = []
        for selector in selectors:
            bounds = selector.contiguous_bounds()
            if bounds is None:
                return None
            out.extend(bounds)
        return out

    def _fast_tslice_params(self, selectors) -> Optional[Tuple[int, int, int]]:
        if self.ndim != 3 or self._ignore or self._view_mode == 'unsorted':
            return None
        if not self._fast_read or not selectors[-1].is_scalar:
            return None

        shape = self._data_shape()
        stepi = selectors[0].full_axis_step(shape[0])
        stepx = selectors[1].full_axis_step(shape[1])
        if stepi is None or stepx is None:
            return None
        return int(selectors[2].indices[0]), stepi, stepx

    def _read_regular(self, idx) -> np.ndarray:
        assert self._ignore, "The data is not regular or shape_hint is not given"

        shp = self._data_shape()
        assert len(shp) in (3, 4), f"shape_hint must be 3D or 4D, got {shp}"

        if self.ndim == 3:
            assert len(shp) == 3, f"ndim=3 but shape_hint={shp}"
            nx, ny, nt = shp
            ab, ae = 0, 1
            xb, xe, yb, ye, tbeg, tend = idx
        elif self.ndim == 4:
            assert len(shp) == 4, f"ndim=4 but shape_hint={shp}"
            _, nx, ny, nt = shp
            ab, ae, xb, xe, yb, ye, tbeg, tend = idx
        else:
            raise ValueError(
                f"_read_regular only supports ndim=3 or 4, got {self.ndim}")

        dt = tend - tbeg

        def base_index(ia, ix):
            return ia * (nx * ny) + ix * ny

        oshape = (ae - ab, xe - xb, ye - yb, dt)
        bytes_total = np.prod(oshape) * 4
        show_bar = bytes_total >= (1 << 30) if self._show_progress else False

        out = np.zeros(oshape, dtype=np.float32)

        iterator = product(range(ab, ae), range(xb, xe))
        total_iters = oshape[0] * oshape[1]
        if show_bar:
            iterator = tqdm(iterator, total=total_iters)

        for ia, ix in iterator:
            beg = base_index(ia, ix) + yb
            end = base_index(ia, ix) + ye
            block = self._segy.collect(beg, end, tbeg, tend)
            out[ia - ab, ix - xb, :, :] = block

        if self.ndim == 3:
            out = out[0]
        return out

    def _spatial_trace_indices(self, selectors) -> np.ndarray:
        shape = tuple(selector.size for selector in selectors)
        if any(size == 0 for size in shape):
            return np.empty(shape, dtype=np.int32)

        grids = np.meshgrid(*[selector.indices for selector in selectors],
                            indexing='ij')
        flat_coords = [grid.reshape(-1) for grid in grids]

        if self._ignore or self._view_mode != 'unsorted':
            tidx = np.ravel_multi_index(flat_coords,
                                        self._data_shape()[:-1]).astype(
                                            np.int32)
        else:
            if not self.is_create_geometry:
                raise RuntimeError(
                    "Need create the geometry first, please call `update_geometry` first"
                )
            coord_grid = np.stack(flat_coords, axis=1)
            tidx = np.asarray(self.map_to_indices(coord_grid), dtype=np.int32)

        return tidx.reshape(shape)

    def _read_trace_window(self, tidx, time_selector: AxisSelector) -> np.ndarray:
        tidx = np.asarray(tidx, dtype=np.int32).reshape(-1)
        ns = time_selector.size
        if tidx.size == 0 or ns == 0:
            return np.empty((tidx.size, ns), dtype=np.float32)

        bounds = time_selector.contiguous_bounds()
        if bounds is not None:
            tb, te = bounds
            if np.any(tidx < 0):
                return self._collect_with_valid_indices(tidx, (tidx.size, ),
                                                        tb, te)
            return self._segy.collect(tidx, tb, te)

        lo, hi, rel = time_selector.dense_window()
        if np.any(tidx < 0):
            block = self._collect_with_valid_indices(tidx, (tidx.size, ), lo,
                                                     hi)
        else:
            block = self._segy.collect(tidx, lo, hi)
        return block[:, rel]

    def _write_trace_window(self, tidx, time_selector: AxisSelector,
                            data: np.ndarray) -> None:
        tidx = self._require_writable_trace_indices(tidx)
        ns = time_selector.size
        if tidx.size == 0 or ns == 0:
            return

        payload = np.ascontiguousarray(data.reshape(tidx.size, ns),
                                       dtype=np.float32)
        trace_bounds = AxisSelector(tidx).contiguous_bounds()
        time_bounds = time_selector.contiguous_bounds()

        if time_bounds is not None:
            tb, te = time_bounds
            if trace_bounds is not None:
                self._segy.write_traces(payload, trace_bounds[0],
                                        trace_bounds[1], tb, te)
            else:
                self._segy.write_traces(payload, tidx, tb, te)
            return

        lo, hi, rel = time_selector.dense_window()
        block = self._segy.collect(tidx, lo, hi)
        block[:, rel] = payload
        if trace_bounds is not None:
            self._segy.write_traces(block, trace_bounds[0], trace_bounds[1],
                                    lo, hi)
        else:
            self._segy.write_traces(block, tidx, lo, hi)

    def _read_base(self, selectors) -> np.ndarray:
        base_shape = tuple(selector.size for selector in selectors)
        if any(size == 0 for size in base_shape):
            return np.empty(base_shape, dtype=np.float32)

        if self.ndim == 2:
            data = self._read_trace_window(selectors[0].indices, selectors[1])
            return data.reshape(base_shape)

        if self._ignore:
            bounds = self._selectors_to_bounds(selectors)
            if bounds is not None:
                return self._read_regular(bounds)

        if self._view_mode != 'unsorted':
            tslice_params = self._fast_tslice_params(selectors)
            if tslice_params is not None:
                t, stepi, stepx = tslice_params
                return self._segy.read_tslice(t, stepi, stepx)[..., None]

            bounds = self._selectors_to_bounds(selectors)
            if bounds is not None:
                if self.ndim == 3:
                    return self._segy.read3d(*bounds)
                return self._segy.read4d(*bounds)

        spatial = self._spatial_trace_indices(selectors[:-1])
        data = self._read_trace_window(spatial.reshape(-1), selectors[-1])
        return data.reshape(*spatial.shape, selectors[-1].size)

    def _write_base(self, selectors, data: np.ndarray) -> None:
        if any(selector.size == 0 for selector in selectors):
            return

        data = np.ascontiguousarray(data, dtype=np.float32)

        if self.ndim == 2:
            self._write_trace_window(selectors[0].indices, selectors[1], data)
            return

        bounds = self._selectors_to_bounds(selectors)
        if bounds is not None and not self._ignore and self._view_mode != 'unsorted':
            if self.ndim == 3:
                self._segy.write3d(data, *bounds)
            else:
                self._segy.write4d(data, *bounds)
            return

        spatial = self._spatial_trace_indices(selectors[:-1]).reshape(-1)
        spatial = self._require_writable_trace_indices(spatial)
        payload = np.ascontiguousarray(
            data.reshape(spatial.size, selectors[-1].size), dtype=np.float32)
        self._write_trace_window(spatial, selectors[-1], payload)


class CheckMixin:

    def _check_bound(self, dim, ib, ie):
        if ie == None:
            assert isinstance(ib, np.ndarray) and ib.ndim == 1, f"if ie is None, ib must be a 1D numpy array"
            assert ib.min() >= 0 and ib.max() < self.shape[dim], f"index array out of range in dim {dim}"
        else:
            assert ib >= 0 and ib < ie and ie <= self.shape[dim], f"index out of range in dim {dim}"


    def _check_bound2(self, idx):
        assert len(idx) == self.ndim * 2, f"ndim is {self.ndim}, need {self.ndim*2} idx, but got {len(idx)}"
        for i in range(self.ndim):
            self._check_bound(i, idx[i*2], idx[i*2+1])


    def _check_bound_idx(self, index, dim):
        assert isinstance(index, (int, np.integer)), "index must be int"
        assert index >= 0 and index < self.shape[dim], f"In dimension {dim}, index {index} out of range"

    def _check_wmode_lastdim(self, idx):
        if idx[-1] is None:
            raise TypeError("Indexing with a list or ndarray for the last dimension is not supported when writing.")

    def _check_data_shape(self, data: np.ndarray, idx):
        dstshape = []
        for i in range(self.ndim):
            if idx[i * 2 + 1] is None:
                dstshape.append(np.asarray(idx[i * 2]).size)
            else:
                dstshape.append(idx[i * 2 + 1] - idx[i * 2])
        dstshape = [k for k in dstshape if k != 1]
        assert tuple(dstshape) == data.squeeze().shape, "shape of the input data is not match the shape of the data to write"


class SegyCMixin:
    _segy: _CXX_SEGY.Pysegy

    def cut(self,
            outname: str,
            ranges: List,
            as2d: bool = False,
            textual: str = '') -> None:
        """
        cut the SEG-Y file into a new SEG-Y file
        """
        if as2d:
            assert len(ranges) == 4, "The length of `ranges` must be 4, as you set the `as2d` to True"
        else:
            assert len(ranges) == self.ndim * 2, f"ndim is {self.ndim}, need {self.ndim*2} idx, but got {len(ranges)}"

        textual = createtool.generate_textual(self._metainfo, textual)

        self._segy.cut(outname, ranges, as2d, textual)

    def create_by_sharing_header(self,
                                 outname: str,
                                 src,
                                 shape=None,
                                 *,
                                 start=None,
                                 as2d=False,
                                 textual='') -> None:
        if isinstance(src, np.ndarray):
            if shape is not None:
                warnings.warn("The shape is ignored, as the src is ndarray")
            shape = src.shape
        elif shape is None:
            raise ValueError("src (filename) is not ndarray, shape is None, need to specify the shape")

        if as2d:
            assert len(shape) == 2, "The shape must be 2D, as you set the `as2d` to True"

        if start is None:
            start = [0] * len(shape)

        if textual not in ("", None):
            textual = createtool.generate_textual(self._metainfo, textual)

        self._segy.show_progress(True)
        if isinstance(src, np.ndarray):
            src = np.ascontiguousarray(src, dtype=np.float32)
            self._segy.create_by_sharing_header(outname, src, start, as2d, textual)
        else:
            self._segy.create_by_sharing_header(outname, src, shape, start, as2d, textual)
        self._segy.show_progress(False)
        


class AccessMixin:
    _segy: _CXX_SEGY.Pysegy

    def __getattr__(self, name):
        if name in ['iline', 'xline', 'offset', 'coordx', 'coordy', 'itrace']:
            return self._SegyAccessor(self._segy, name, getattr(self, '_mode', 'r'))

        return super().__getattr__(name)


    def __setattr__(self, name: str, value) -> None:
        if name in ['iline', 'xline', 'offset', 'coordx', 'coordy', 'itrace']:
            accessor = self._SegyAccessor(self._segy, name, getattr(self, '_mode', 'r'))
            accessor[:] = value
            return
        else:
            super().__setattr__(name, value)

    class _SegyAccessor:
        _HEADER_KEYLOCS = {
            'iline': 'iline',
            'xline': 'xline',
            'offset': 'offset',
            'coordx': 'xloc',
            'coordy': 'yloc',
        }

        def __init__(self, segy, attribute, mode='r'):
            assert mode in ['r', 'rw']
            self._segy = segy
            self.attribute = attribute
            self._mode = mode

        def __getitem__(self, index):
            index = self._process_index(index)
            method = getattr(self._segy, self.attribute)
            if isinstance(index, (int, np.integer)):
                return method(index)
            return self._batch_get(index)

        def __setitem__(self, index, value):
            if self._mode == 'r':
                raise RuntimeError("The SEG-Y file is not writable, as you set the `mode` to 'r'. If you want to enable write mode, set to `rw`") # yapf: disable
            index = self._process_index(index)
            if isinstance(index, (int, np.integer)):
                self._set_scalar(index, value)
                return
            self._set_batch(index, value)

        def _batch_get(self, index: np.ndarray):
            if index.size == 0:
                if self.attribute == 'itrace':
                    return np.empty((0, self._segy.nt), dtype=np.float32)
                return np.empty((0,), dtype=np.int32)
            if self.attribute == 'itrace':
                return self._collect_traces(index)
            keyloc = self._segy.get_keylocs()[self._HEADER_KEYLOCS[self.attribute]]
            return self._get_trace_keys(index, [keyloc], [4]).reshape(-1)

        def _set_scalar(self, index, value):
            if self.attribute == 'itrace':
                trace = np.ascontiguousarray(np.asarray(value, dtype=np.float32))
                self._segy.write_itrace(trace, int(index))
                return

            method = getattr(self._segy, f"set_{self.attribute}")
            method(int(index), value)

        def _set_batch(self, index: np.ndarray, value):
            if index.size == 0:
                return

            if self.attribute == 'itrace':
                traces = np.ascontiguousarray(np.asarray(value, dtype=np.float32))
                if self._is_contiguous_range(index):
                    self._segy.write_traces(
                        traces,
                        int(index[0]),
                        int(index[-1]) + 1,
                        0,
                        self._segy.nt,
                    )
                else:
                    self._segy.write_traces(
                        traces,
                        index.astype(np.int32, copy=False),
                        0,
                        self._segy.nt,
                    )
                return

            values = np.asarray(value)
            if values.ndim == 0:
                values = np.full(index.shape, values.item())
            else:
                values = values.reshape(-1)
            if values.size != index.size:
                raise ValueError(
                    "The number of values must match the number of indices")

            method = getattr(self._segy, f"set_{self.attribute}")
            for i, v in zip(index, values):
                method(int(i), v.item() if isinstance(v, np.generic) else v)

        def _collect_traces(self, index: np.ndarray) -> np.ndarray:
            if self._is_contiguous_range(index):
                return self._segy.collect(
                    int(index[0]),
                    int(index[-1]) + 1,
                    0,
                    self._segy.nt,
                )
            return self._segy.collect(index.astype(np.int32, copy=False), 0,
                                      self._segy.nt)

        def _get_trace_keys(self, index: np.ndarray, keys, length) -> np.ndarray:
            if self._is_contiguous_range(index):
                return self._segy.get_trace_keys(
                    keys,
                    length,
                    int(index[0]),
                    int(index[-1]) + 1,
                )
            return self._segy.get_trace_keys(keys, length,
                                             index.astype(np.int32, copy=False))

        def _is_contiguous_range(self, index: np.ndarray) -> bool:
            return index.size > 0 and np.all(index[1:] == index[:-1] + 1)

        def _process_index(self, keys):
            if isinstance(keys, tuple):
                # Raise an error if more than one dimension is indexed
                raise IndexError("Only 1D indexing is supported.")

            if isinstance(keys, (int, np.integer)):
                if keys < 0:
                    keys += self._segy.ntrace
                return keys

            elif isinstance(keys, slice):
                # Slice index
                start = 0 if keys.start is None else keys.start
                stop = self._segy.ntrace if keys.stop is None else keys.stop
                step = 1 if keys.step is None else keys.step
                if start < 0:
                    start += self._segy.ntrace
                if stop < 0:
                    stop += self._segy.ntrace
                keys = np.arange(start, stop, step)
                if keys.size == 0:
                    return keys
                assert keys.min() >= 0 and keys.max() < self._segy.ntrace, "Index out of range"
                return keys

            elif isinstance(keys, list) or isinstance(keys, np.ndarray):
                # List or numpy array of indices
                keys = np.array(keys).squeeze()
                assert keys.ndim == 1, "Only 1D indexing is supported."
                if keys.size == 0:
                    return keys
                keys = np.where(keys < 0, keys + self._segy.ntrace, keys)
                assert keys.min() >= 0 and keys.max() < self._segy.ntrace, "Index out of range"
                return keys

            else:
                raise TypeError("Unsupported index type.")


        def __repr__(self) -> str:
            d = np.array([self.__getitem__(i) for i in [0, 1, -2, -1]])
            if self.attribute == 'itrace':
                usage = "traces: \n"
                for i in range(4):
                    idx = i
                    if idx > 1:
                        idx = self._segy.ntrace - 4 + i
                    usage += f"trace {idx}: [{d[i, 0]:.4f}, {d[i, 1]:.4f}, {d[i, 2]:.4f}, {d[i, 3]:.4f}, ...]\n"
                    if i == 1:
                        usage += "...\n"
            else:
                usage = f"{self.attribute}: [{d[0]}, {d[1]}, ..., {d[2]}, {d[3]}]"
            return usage

    def bkeyi2(self, loc):
        """get binary header value at loc, view it as int16_t"""
        return self._segy.bkeyi2(loc)

    def bkeyi4(self, loc):
        """get binary header value at loc, view it as int32_t"""
        return self._segy.bkeyi4(loc)

    def keyi2(self, idx, loc):
        """get idx-th trace header value at loc, view it as int16_t"""
        return self._segy.keyi2(idx, loc)

    def keyi4(self, idx, loc):
        """get idx-th trace header value at loc, view it as int16_t"""
        return self._segy.keyi4(idx, loc)

    def set_bkeyi2(self, loc, value):
        """set binary header value at loc, view it as int16_t"""
        return self._segy.set_bkeyi2(loc, value)

    def set_bkeyi4(self, loc, value):
        """set binary header value at loc, view it as int32_t"""
        return self._segy.set_bkeyi4(loc, value)

    def set_keyi2(self, idx, loc, value):
        """set idx-th trace header value at loc, view it as int16_t"""
        return self._segy.set_keyi2(idx, loc, value)

    def set_keyi4(self, idx, loc, value):
        """set idx-th trace header value at loc, view it as int16_t"""
        return self._segy.set_keyi4(idx, loc, value)

    def textual_header(self, code='u', printtext=True):
        """get textual header"""
        t = self._segy.textual_header(code)
        if printtext:
            print(t)
        else:
            return t
# fmt: on


class InnerMixin:
    _segy: _CXX_SEGY.Pysegy

    def __array__(self):
        """To support np.array(SegyNP(xxx))"""
        return self.to_numpy()

    def __getitem__(self, slices) -> np.ndarray:
        visible_index = self._normalize_key(slices)
        selectors = visible_index.selectors

        if self._T:
            base = self._read_base(selectors[::-1])
            axes = tuple(range(base.ndim - 1, -1, -1))
            base = np.transpose(base, axes=axes)
        else:
            base = self._read_base(selectors)

        return visible_index.format_result(base)

    def __setitem__(self, slices, data: np.ndarray) -> None:
        visible_index = self._normalize_key(slices)
        payload = visible_index.prepare_value(data, dtype=np.float32)

        if self._T:
            axes = tuple(range(payload.ndim - 1, -1, -1))
            payload = np.transpose(payload, axes=axes)
            return self._write_base(visible_index.selectors[::-1], payload)

        return self._write_base(visible_index.selectors, payload)

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
        if hasattr(self, '_segy'):
            self._segy.close()

    def __repr__(self) -> str:
        out = f"cigsegy.SegyNP class, file name: '{self.file_name}'\n\n"
        if self.ndim == 2 or self._metainfo is None:
            meta = utils.make_2d_meta(self._segy)
        else:
            meta = {**self._keylocs, **self._metainfo}
            meta = utils.post_process_meta(self._segy, meta)
        out += utils.parse_metainfo(meta)
        return out

    def __len__(self) -> int:
        return self.shape[0]


class SegyNP(InnerMixin, RWMixin, InterpMixin, PlotMixin, GeometryMixin,
             ScanMixin, CheckMixin, AccessMixin, SegyCMixin):
    _VALID_VIEW_MODES = {'scan', 'lazy', '2d', 'unsorted'}

    def __init__(self,
                 filename: str,
                 keylocs: dict = None,
                 mode: str = 'r',
                 *,
                 ndim: int = None,
                 shape_hint: tuple = None,
                 as_unsorted: bool = False,
                 fast_read: bool = False,
                 show_progress: bool = False,
                 keys: dict = None,
                 view_mode: str = None) -> None:
        np.set_printoptions(suppress=True)

        assert mode in ['r', 'rw'], "`mode` only can be 'r' or 'rw'"
        self._ndim = ndim
        self._fname = filename
        self._requested_keylocs = keylocs
        self._requested_keys = keys
        self._requested_ndim = ndim

        self._show_progress = show_progress
        self._segy = _CXX_SEGY.Pysegy(str(filename), mode=='rw')
        self._segy.show_progress(False)
        self._mode = mode
        if self._mode == 'rw':
            warnings.warn(
                "\n**Dangerous!!!** You are using a writable mode ('rw'), which may **alter** the SEG-Y file. \n"
                "It is strongly recommended to make a **backup copy** of the file before proceeding "
                "to avoid any potential irreversible changes.", UserWarning)

        self._fast_read = fast_read  # TODO: set as 'auto'? and how to set step?

        # for values
        self._min = None
        self._max = None

        # for coordinates transform
        self._trans_matrix = None
        self._geometry = None
        self._north = None

        self._keylocs = None
        self._metainfo = None
        self._lineinfo = None
        self._fstep = (2, 2)

        self._shape2 = (self._segy.ntrace, self._segy.nt)
        self._shape3 = None
        self._ignore = False # whether ignore 'scan', only valid when the data is regular and shape_hint is given

        self._T = False

        self._view_mode = self._resolve_view_mode(view_mode, ndim, shape_hint,
                                                  as_unsorted)
        self._unsorted = self._view_mode == 'unsorted'

        if shape_hint is not None:
            if view_mode not in (None, 'lazy'):
                raise ValueError(
                    "shape_hint can only be used with view_mode='lazy' or by leaving view_mode unset"
                )
            if len(shape_hint) < 3 or len(shape_hint) > 4:
                raise ValueError("shape_hint must be of length 3 or 4")
            if shape_hint[-1] != self._shape2[-1]:
                raise ValueError(f"shape_hint's last dimension must be equal to nt ({self._shape2[-1]}), but got {shape_hint[-1]}") # yapf: disable 
            if np.prod(shape_hint) != np.prod(self._shape2):
                raise ValueError(f"The product shape_hint ({np.prod(shape_hint)}) != (ntrace, nt) ({np.prod(self._shape2)})") # yapf: disable
            self._ndim = len(shape_hint)
            self._shape3 = shape_hint
            self._ignore = True
            self._view_mode = 'lazy'
        elif self._view_mode == 'unsorted':
            self._scan_unsorted(keylocs, keys, ndim)
        elif self._view_mode == 'scan':
            try:
                self._scan(keylocs)
            except Exception as e:
                if ndim in (None, 2):
                    self._init_as_2d(str(e))
                else:
                    raise RuntimeError(f"{str(e)}\n This SEG-Y file may be unsorted, you can pass `as_unsorted` to view it as unsorted file, but it may be slow") from e # yapf: disable
            if ndim is not None and self.ndim != ndim:
                raise RuntimeError(f"You set ndim as {ndim}, but the SEG-Y file's ndim is {self.ndim}") # yapf: disable
        else:
            self._init_as_2d()

    def _resolve_view_mode(self, view_mode, ndim, shape_hint, as_unsorted):
        if view_mode is not None:
            if view_mode not in self._VALID_VIEW_MODES:
                raise ValueError(
                    f"view_mode must be one of {sorted(self._VALID_VIEW_MODES)}, got {view_mode!r}"
                )
            if as_unsorted and view_mode != 'unsorted':
                warnings.warn(
                    f"`as_unsorted=True` is ignored because `view_mode='{view_mode}'`"
                )
            resolved = view_mode
        else:
            if as_unsorted:
                resolved = 'unsorted'
            elif shape_hint is not None:
                resolved = 'lazy'
            elif ndim == 2:
                resolved = '2d'
            else:
                resolved = 'scan'

        if resolved == 'unsorted' and ndim == 2:
            warnings.warn("`ndim` is 2, so unsorted view will be ignored")
            return '2d'

        return resolved

    def _init_as_2d(self, reason: str = None):
        if reason:
            warnings.warn(
                f"Failed to create 3D/4D geometry, fallback to 2D trace view: {reason}"
            )
        self._ndim = 2
        self._shape3 = None
        self._metainfo = utils.make_2d_meta(self._segy)
        self._keylocs = {
            'iline': None,
            'xline': None,
            'offset': None,
            'istep': 1,
            'xstep': 1,
            'ostep': 1,
            'xloc': self._metainfo['xloc'],
            'yloc': self._metainfo['yloc'],
        }


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
            out = list(self._shape2)
        else:
            out = list(self._shape3)
        if self._T:
            out = out[::-1]
        return tuple(out)

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
            N = self.ntrace // 3
            xy1 = [self.coordx[N], self.coordy[N]]
            xy2 = xy1.copy()
            xy2[1] += 1000
            n = self.xy_to_ix(np.array([xy1, xy2]))
            di = n[1] - n[0]
            self._north = di / np.linalg.norm(di)

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
        if hasattr(self, '_segy'):
            self._segy.close()

    def to_2d(self):
        """
        Treat the SEG-Y file as a collection of traces, shape is like (ntrace, nt)
        """
        self._ndim = 2
        self._view_mode = '2d'

    def to_nd(self):
        """
        Treat the SEG-Y file as a 3D/4D array. If the SEG-Y file is scanned, the keylocs will be ignored 
        """
        if self._ignore:
            self._ndim = len(self._shape3)
            self._view_mode = 'lazy'
            return

        if self._view_mode == 'unsorted':
            if self._shape3 is None:
                self._scan_unsorted(self._requested_keylocs,
                                    self._requested_keys,
                                    self._requested_ndim)
            self._ndim = len(self._shape3)
            self._view_mode = 'unsorted'
            return

        if self._shape3 is None:
            self._scan(self._requested_keylocs)
        self._ndim = self._segy.ndim
        self._view_mode = 'scan'

    def max(self, real=False) -> float:
        """
        return the maximum value of the data, 
        if real is False, the maximum value is not real, and we just read a part of traces to calculate max
        if real is True, we read all traces to calculate max
        """
        if real:
            return self[...].max()
        if self._max is None:
            self._eval_range()
        return self._max

    def min(self, real=False) -> float:
        """
        return the min value of the data, 
        if real is False, the min value is not real, and we just read a part of traces to calculate min
        if real is True, we read all traces to calculate min
        """
        if real:
            return self[...].min()
        if self._min is None:
            self._eval_range()
        return self._min

    @property
    def lineinfo(self):
        return self._lineinfo

    @property
    def metainfo(self):
        return self._metainfo

    @property
    def keylocs(self):
        return self._keylocs

    @property
    def access_mode(self) -> str:
        return self._mode

    @property
    def view_mode(self) -> str:
        return self._view_mode

    @property
    def fast_read(self) -> bool:
        return self._fast_read

    @fast_read.setter
    def fast_read(self, value: bool) -> None:
        self._fast_read = value

    def set_fast_read_steps(self, istep, xstep):
        self._fstep = (istep, xstep)

    @property
    def T(self):
        self._T = not self._T
        return self
