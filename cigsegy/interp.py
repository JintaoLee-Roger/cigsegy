# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
"""
Interpolation functions for extracting arbitrary lines, resampling data, and data fusion.

This module provides a set of interpolation functions that are essential for handling seismic data 
and other geophysical datasets. These functions enable the extraction of arbitrary lines from data grids, 
resampling of data to different resolutions, and the fusion of multiple datasets.
"""

import numpy as np
from .transform import apply_transform
from cigsegy import ExceptionWrapper

try:
    from numba import njit
except ImportError:

    def njit(func):

        def wrapper(*args, **kwargs):
            print("Warning: Numba is not installed. This function may run slowly without Numba's optimization. To significantly improve performance, please install Numba by running `pip install numba`.") # yapf: disable
            return func(*args, **kwargs)

        return wrapper


try:
    from scipy.spatial import cKDTree
except BaseException as E:
    cKDTree = ExceptionWrapper(E, "run `pip install scipy` to install the the dependency, becuase we need use cKDTree") # yapf: disable


def merge_data(d1, d2, start2, useA=True):
    """
    Merge two data with different start points.
    coordinates range of d1 is (0, n1), (0, n2)
    coordinates range of d2 is (start2[0], start2[0]+n1), (start2[1], start2[1]+n2)
    """
    assert d1.ndim == 3 and d2.ndim == 3
    assert d1.shape[2] == d2.shape[2]

    mask1 = ~np.all(d1 == 0, axis=2)
    mask2 = ~np.all(d2 == 0, axis=2)
    n11, n21, n3 = d1.shape
    n12, n22, _ = d2.shape
    x0, y0 = start2

    x_min = min(0, x0)
    x_max = max(n11, x0 + n12)
    y_min = min(0, y0)
    y_max = max(n21, y0 + n22)

    nx_total = x_max - x_min
    ny_total = y_max - y_min

    x1s, y1s = -x_min, -y_min
    x2s, y2s = x0 - x_min, y0 - y_min

    mmask1 = np.zeros((nx_total, ny_total), dtype=bool)
    mmask2 = np.zeros((nx_total, ny_total), dtype=bool)
    mmask1[x1s:x1s + n11, y1s:y1s + n21] = mask1
    mmask2[x2s:x2s + n12, y2s:y2s + n22] = mask2

    if useA:
        merged_mask = np.logical_not(mmask1) & mmask2
    else:
        merged_mask = mmask2

    out = np.zeros((nx_total, ny_total, n3), dtype=d1.dtype)
    out[x1s:x1s + n11, y1s:y1s + n21] = d1

    # @njit
    # def _merge(out, d2, merged_mask, n12, n22, x2s, y2s):
    #     for i in range(n12):
    #         for j in range(n22):
    #             ii = i + x2s
    #             jj = j + y2s
    #             if merged_mask[ii, jj]:
    #                 out[ii, jj] = d2[i, j]
    #     return out
    # out = _merge(out, d2, merged_mask, n12, n22, x2s, y2s)

    I_grid, J_grid = np.meshgrid(np.arange(n12), np.arange(n22), indexing='ij')
    II = I_grid + x2s
    JJ = J_grid + y2s

    valid_mask = (II >= 0) & (II < out.shape[0]) & (JJ >= 0) & (JJ
                                                                < out.shape[1])

    combined_mask = merged_mask[II, JJ] & valid_mask

    I_selected = I_grid[combined_mask]
    J_selected = J_grid[combined_mask]
    II_selected = II[combined_mask]
    JJ_selected = JJ[combined_mask]

    out[II_selected, JJ_selected, :] = d2[I_selected, J_selected, :]

    return out


def align_coordinates(
    tgt: np.ndarray,
    transr: np.ndarray,
    transt: np.ndarray,
    xyicr: np.ndarray,
    xyict: np.ndarray,
    order: int = 1,
):
    """
    0. information of ref and tgt
        ref: xyr, icr_r, datar
        tgt: xyt, ict_t, datat
        icr_r means the ic coordinates of ref in ref geometry
        ict_r means the ic coordinates of tgt in ref geometry
    
    1. convert xyt to ict_r 
    2. generate gridt_r based on ict_r
    3. convert gridt_r to xy_gt
    4. convert xy_gt to ic_gt_t in tgt geometry
    5. interpolate data_gt_t using datat
    6. reshape data_gt_t to meet the 3d shape

    Parameters
    ----------
    tgt : SegyNP or np.ndarray
        target data
    transr : np.ndarray
        transform matrix from ic to xy of ref
    transt : np.ndarray
        transform matrix from ic to xy of tgt
    xyicr : np.ndarray
        shape is (N, 4), each row is [x, y, iline, xline] in ref geometry
    xyict : np.ndarray
        shape is (N, 4), each row is [x, y, iline, xline] in tgt geometry

    Returns
    -------
    np.ndarray
        shape is (N, n3)
    offset : List
        [x0, y0]
    """
    n1, n2, n3 = tgt.shape
    # 0. information of ref and tgt
    xyr = xyicr[:, :2]
    icr_r = xyicr[:, 2:]
    xyt = xyict[:, :2]
    ict_t = xyict[:, 2:]

    # 1. convert xyt to ict_r (ic of tgt in ref geometry)
    ict_r = apply_transform(xyt, transr, True)

    # 2. generate gridt_r based on ict_r  (grid of tgt in ref geometry)
    start1, end1 = round(ict_r[:, 0].min()), round(ict_r[:, 0].max())
    start2, end2 = round(ict_r[:, 1].min()), round(ict_r[:, 1].max())
    X, Y = np.meshgrid(np.arange(start1, end1 + 1),
                       np.arange(start2, end2 + 1),
                       indexing='ij')
    gridt_r = np.c_[X.flatten(), Y.flatten()]

    # we need know the start point of tgt compared to ref (assume ref start at 0, 0)
    shift_x = gridt_r[:, 0].min() - icr_r[:, 0].min()
    shift_y = gridt_r[:, 1].min() - icr_r[:, 1].min()

    # using kd-tree to filter the grid points (avoid interpolation on the blank area)
    tree = cKDTree(ict_r)
    distances, indices = tree.query(gridt_r, k=1)
    mask = distances <= 1  # using distance to filter the grid points
    to_interp = gridt_r[mask]

    # 3. convert gridt_r to xy_gt
    xy_gt = apply_transform(to_interp, transr, False)

    # 4. convert xy_gt to ic_gt_t in tgt geometry
    ic_gt_t = apply_transform(xy_gt, transt, True)

    # 5. interpolate data_gt_t using datat
    ic_gt_t[:, 0] -= ict_t[:, 0].min()
    ic_gt_t[:, 1] -= ict_t[:, 1].min()

    # another mask to filter the grid points, avoid interpolation out of the boundary
    mask2 = (ic_gt_t[:, 0] >= 0) & (ic_gt_t[:, 0] < n1-1) & \
            (ic_gt_t[:, 1] >= 0) & (ic_gt_t[:, 1] < n2-1)
    ic_gt_t = ic_gt_t[mask2]
    to_interp = to_interp[mask2]

    # interpolate data_gt_t using datat
    tgt_interp = interp3d(tgt, ic_gt_t, order)

    # 6. reshape data_gt_t to meet the 3d shape
    # the shape of output is (nx, ny, n3)
    nx = end1 - start1 + 1
    ny = end2 - start2 + 1

    output = np.zeros((nx, ny, n3), dtype=tgt.dtype)

    ix = (to_interp[:, 0] - start1).astype(int)
    iy = (to_interp[:, 1] - start2).astype(int)

    output[ix, iy, :] = tgt_interp

    return output, [shift_x, shift_y]


@njit
def interp3d(volume: np.ndarray, points: np.ndarray, order: int = 1):
    """
    3D body data is interpolated at a given 2D point.

    Parameters
    -----------
    volume : np.ndarray
        shape is (n1, n2, n3)
    points : np.ndarray
        shape is (N, 2), each row is [x, y]
    order : int
        order of interpolate, can be one of 1, 2, 3。

    Returns
    --------
    out : np.ndarray
        shape is (N, n3)
    """
    n1, n2, n3 = volume.shape
    N = points.shape[0]
    out = np.zeros((N, n3), dtype=volume.dtype)

    for idx in range(N):
        x = points[idx, 0]
        y = points[idx, 1]

        # Ensure that x and y are within the valid range
        # This also means the boundary is interpolated by the repeated values
        x = min(max(x, 0.0), n1 - 1.000001)
        y = min(max(y, 0.0), n2 - 1.000001)

        if order == 1:
            x0 = int(np.floor(x))
            y0 = int(np.floor(y))
            x1 = x0 + 1
            y1 = y0 + 1

            x1 = min(x1, n1 - 1)
            y1 = min(y1, n2 - 1)

            x_frac = x - x0
            y_frac = y - y0

            w00 = (1 - x_frac) * (1 - y_frac)
            w01 = (1 - x_frac) * y_frac
            w10 = x_frac * (1 - y_frac)
            w11 = x_frac * y_frac

            v00 = volume[x0, y0, :]
            v01 = volume[x0, y1, :]
            v10 = volume[x1, y0, :]
            v11 = volume[x1, y1, :]

            out[idx, :] = w00 * v00 + w01 * v01 + w10 * v10 + w11 * v11

        elif order == 2:
            # obtain the 9 grid points around the target point
            x_indices = np.array(
                [int(np.floor(x) - 0.5) + i for i in range(3)])
            y_indices = np.array(
                [int(np.floor(y) - 0.5) + i for i in range(3)])

            # Ensure that the indices are within the valid range
            x_indices = np.clip(x_indices, 0, n1 - 1)
            y_indices = np.clip(y_indices, 0, n2 - 1)

            # Obtain the node coordinates
            x_nodes = x_indices.astype(np.float64)
            y_nodes = y_indices.astype(np.float64)

            # Calculate the Lagrange basis functions
            Lx = np.zeros(3)
            Ly = np.zeros(3)

            for i in range(3):
                xi = x_nodes[i]
                Li = 1.0
                for j in range(3):
                    if i != j:
                        xj = x_nodes[j]
                        Li *= (x - xj) / (xi - xj + 1e-10)
                Lx[i] = Li

            for i in range(3):
                yi = y_nodes[i]
                Li = 1.0
                for j in range(3):
                    if i != j:
                        yj = y_nodes[j]
                        Li *= (y - yj) / (yi - yj + 1e-10)
                Ly[i] = Li

            # Calculate the interpolated value
            value = np.zeros(n3, dtype=volume.dtype)
            for i in range(3):
                xi = x_indices[i]
                for j in range(3):
                    yj = y_indices[j]
                    w = Lx[i] * Ly[j]
                    v = volume[xi, yj, :]
                    value += w * v
            out[idx, :] = value

        elif order == 3:
            # obtain the 16 grid points around the target point
            x_indices = np.array([int(np.floor(x) - 1) + i for i in range(4)])
            y_indices = np.array([int(np.floor(y) - 1) + i for i in range(4)])

            # Ensure that the indices are within the valid range
            x_indices = np.clip(x_indices, 0, n1 - 1)
            y_indices = np.clip(y_indices, 0, n2 - 1)

            # Obtain the node coordinates
            x_nodes = x_indices.astype(np.float64)
            y_nodes = y_indices.astype(np.float64)

            # Calculate the Lagrange basis functions
            Lx = np.zeros(4)
            Ly = np.zeros(4)

            for i in range(4):
                xi = x_nodes[i]
                Li = 1.0
                for j in range(4):
                    if i != j:
                        xj = x_nodes[j]
                        Li *= (x - xj) / (xi - xj + 1e-10)
                Lx[i] = Li

            for i in range(4):
                yi = y_nodes[i]
                Li = 1.0
                for j in range(4):
                    if i != j:
                        yj = y_nodes[j]
                        Li *= (y - yj) / (yi - yj + 1e-10)
                Ly[i] = Li

            # calculate the interpolated value
            value = np.zeros(n3, dtype=volume.dtype)
            for i in range(4):
                xi = x_indices[i]
                for j in range(4):
                    yj = y_indices[j]
                    w = Lx[i] * Ly[j]
                    v = volume[xi, yj, :]
                    value += w * v
            out[idx, :] = value

        else:
            raise ValueError("Order must be 1, 2, or 3")

    return out


def interp(p: np.ndarray, d: np.ndarray) -> np.ndarray:
    """
    Interpolates a 1d series in 3d space by **linear** interpolation.
    
    Parameters
    ----------
    p : ArrayLike
        (N, 2), each row is [x, y], normalized, i.e., 0 <= x, y <= 1
    d : ArrayLike
        shape is (N, 4, n3)
    """
    x = p[:, 0]
    y = p[:, 1]

    # linear interpolation
    # fmt: off
    out = (d[:, 0, :] * (1 - x)[:, None] * (1 - y)[:, None] +
           d[:, 1, :] * (1 - x)[:, None] * y[:, None] +
           d[:, 2, :] * x[:, None] * (1 - y)[:, None] +
           d[:, 3, :] * x[:, None] * y[:, None])

    # fmt: on

    return out


def interpolate_path(points, di=1):
    """
    Interpolates a path from the given points with a given step size.
    """
    points = np.array(points)

    diffs = np.diff(points, axis=0)
    distances = np.sqrt((diffs**2).sum(axis=1))
    cum_distances = np.concatenate(([0], np.cumsum(distances)))
    total_distance = cum_distances[-1]

    distances_interp = np.arange(0, total_distance + di, di)
    x_points = np.interp(distances_interp, cum_distances, points[:, 0])
    y_points = np.interp(distances_interp, cum_distances, points[:, 1])

    indices = np.searchsorted(distances_interp, cum_distances)

    return np.column_stack((x_points, y_points)), indices


def extract_data(data, p):
    """
    Parameters
    ----------
    segy: str
        segyfile
    p : ArrayLike
        shape is (N, 2)
    """
    n1, n2, n3 = data.shape

    N = p.shape[0]
    xq = p[:, 0]
    yq = p[:, 1]

    i = np.clip(np.searchsorted(np.arange(n1), xq, side="right"), 1, n1 - 1)
    j = np.clip(np.searchsorted(np.arange(n2), yq, side="right"), 1, n2 - 1)

    x0 = i - 1
    x1 = i
    y0 = j - 1
    y1 = j

    # normalized path points
    pout = np.zeros_like(p, np.float32)
    pout[:, 0] = p[:, 0] - i + 1
    pout[:, 1] = p[:, 1] - j + 1

    assert pout.min() >= 0 and pout.max() <= 1

    grid_points = np.stack(
        [
            np.array([x0, y0]).T,
            np.array([x0, y1]).T,
            np.array([x1, y0]).T,
            np.array([x1, y1]).T
        ],
        axis=1,
    )

    all_points = grid_points.reshape(-1, 2)
    unique, inverse = np.unique(all_points, axis=0, return_inverse=True)
    unique_data = np.array([data[i, j] for i, j in unique])

    return pout, unique_data[inverse].reshape(N, 4, n3)


def arbitray_line(data, points, di: float = 1):
    """
    data : SegyNP or np.ndarray
    points : the points to interpolate a path
    di : step
    """
    p, indices = interpolate_path(points, di)
    pout, pdata = extract_data(data, p)
    out = interp(pout, pdata)

    return out, p, indices
