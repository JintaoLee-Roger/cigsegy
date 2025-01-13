# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from typing import List, Tuple
import numpy as np
from cigsegy.cpp import _CXX_SEGY
from cigsegy.transform import get_transform_metrix, apply_transform
from cigsegy.interp import arbitray_line
from cigsegy import utils, plot, tools, createtool
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
            if isinstance(keylocs, List):
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

    def map_to_index(self, index):
        if self._geometry is None:
            raise RuntimeError("geometry is not created, Call `update_geometry` first") # yapf: disable
        if self.ndim == 2:
            raise RuntimeError("ndim is 2, unsupport this function") # yapf: disable

        index = np.array(index)
        assert index.ndim == 2
        if self.ndim == 3:
            assert index.shape[1] == 2
            index = self._geometry[index[:, 0], index[:, 1]]
        elif self.ndim == 4:
            assert index.shape[1] == 3
            index = self._geometry[index[:, 0], index[:, 1], index[:, 2]]
        return index

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
            'zero' means points are taken from a geometry with zero-origin (i.e., numpy array indexes), 
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
        if self.unsorted and self.ndim > 2:
            d = self[...]
        else:
            d = self._segy.read()
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
            self._segy.tofile(fpath, self.ndim == 2)
        self._segy.show_progress(False)

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
            if self._fast_read and self._is_time_slice(idx):
                d = self._segy.read_tslice(idx[4], self._fstep[0], self._fstep[1])
            else:
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

    def _is_time_slice(self, idx):
        # Only used in _read3d and num==0 and not unsorted
        if idx[1] - idx[0] == self.shape[0] and idx[3] - idx[2] == self.shape[1] and idx[5]-idx[4]==1: # yapf: disable
            return True
        return False


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
                dstshape.append(idx[i * 2].shape)
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

        if len(textual) > 0 and len(textual) != 3200:
            textual = createtool.generate_textual(textual)

        self._segy.show_progress(True)
        if isinstance(src, np.ndarray):
            self._segy.create_by_sharing_header(outname, src, start, as2d, textual)
        else:
            self._segy.create_by_sharing_header(outname, src, shape, start, as2d, textual)
        self._segy.show_progress(False)
        


class AccessMixin:
    _segy: _CXX_SEGY.Pysegy

    def __getattr__(self, name):
        if name in ['iline', 'xline', 'offset', 'coordx', 'coordy', 'itrace']:
            return self._SegyAccessor(self._segy, name)

        return super().__getattr__(name)


    def __setattr__(self, name: str, value) -> None:
        if name in ['iline', 'xline', 'offset', 'coordx', 'coordy', 'itrace']:
            accessor = self._SegyAccessor(self._segy, name)
            accessor[:] = value
            return
        else:
            super().__setattr__(name, value)

    class _SegyAccessor:

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
            else:
                return np.array([method(i) for i in index])

        def __setitem__(self, index, value):
            if self._mode == 'r':
                raise RuntimeError("The SEG-Y file is not writable, as you set the `mode` to 'r'. If you want to enable write mode, set to `rw`") # yapf: disable
            index = self._process_index(index)
            # TODO: check the value size?

            if self.attribute == 'itrace':
                method = getattr(self._segy, f"write_{self.attribute}")
            else:
                method = getattr(self._segy, f"set_{self.attribute}")
            if isinstance(index, (int, np.integer)):
                method(index, value)
            else:
                for i, v in zip(index, value):
                    method(i, v)

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
                start, stop, step = keys.start, keys.stop, keys.step
                if keys.start < 0:
                    start += self._segy.ntrace
                if keys.stop < 0:
                    stop += self._segy.ntrace
                keys = np.arange(start, stop, step)
                assert keys.min() >= 0 and keys.max() < self._segy.ntrace, "Index out of range"
                return keys

            elif isinstance(keys, list) or isinstance(keys, np.ndarray):
                # List or numpy array of indices
                keys = np.array(keys).squeeze()
                assert keys.ndim == 1, "Only 1D indexing is supported."
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
        idx = self._process_keys(slices)
        if self._ndim == 4:
            return self._read4d(idx)
        elif self._ndim == 3:
            return self._read3d(idx)
        else:
            return self._read2d(idx)

    def __setitem__(self, slices, data: np.ndarray) -> None:
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
        if hasattr(self, '_segy'):
            self._segy.close()

    def __repr__(self) -> str:
        out = f"cigsegy.SegyNP class, file name: '{self.file_name}'\n\n"
        keys = self._segy.get_keylocs()
        meta = self._segy.get_metainfo()
        meta = {**keys, **meta}
        meta = utils.post_process_meta(self._segy, meta)
        out += utils.parse_metainfo(meta)
        return out

    def __len__(self) -> int:
        return self.shape[0]


class SegyNP(InnerMixin, RWMixin, InterpMixin, PlotMixin, GeometryMixin,
             ScanMixin, CheckMixin, AccessMixin, SegyCMixin):

    def __init__(self,
                 filename: str,
                 keylocs: dict = None,
                 mode: str = 'r',
                 *,
                 ndim: int = None,
                 as_unsorted: bool = False,
                 fast_read: bool = False,
                 keys: dict = None) -> None:
        np.set_printoptions(suppress=True)

        assert mode in ['r', 'rw'], "`mode` only can be 'r' or 'rw'"
        self._ndim = ndim
        self._fname = filename

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
        self._eval_range()

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

        if ndim == 2 and as_unsorted:
            warnings.warn("`ndim` is 2, so `as_unsorted` will be ignored")
            as_unsorted = False
        self._unsorted = as_unsorted

        if ndim is None or ndim != 2:
            if self._unsorted:
                self._scan_unsorted(keylocs, keys, ndim)
            else:
                try:
                    self._scan(keylocs)
                except Exception as e:
                    raise RuntimeError(f"{str(e)}\n This SEG-Y file may be unsorted, you can pass `as_unsorted` to view it as unsorted file, but it may be slow") from e # yapf: disable
            if ndim is not None and self.ndim != ndim:
                raise RuntimeError(f"You set ndim as {ndim}, but the SEG-Y file's ndim is {self.ndim}") # yapf: disable

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
            return tuple(self._shape2)
        else:
            return tuple(self._shape3)

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

    def to_nd(self):
        """
        Treat the SEG-Y file as a 3D/4D array. If the SEG-Y file is scanned, the keylocs will be ignored 
        """
        if self._shape3 is None:
            self._scan()
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
    def fast_read(self) -> bool:
        return self._fast_read

    @fast_read.setter
    def fast_read(self, value: bool) -> None:
        self._fast_read = value

    def set_fast_read_steps(self, istep, xstep):
        self._fstep = (istep, xstep)
