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

    f00 = d[:, 0, :]
    f01 = d[:, 1, :]
    f10 = d[:, 2, :]
    f11 = d[:, 3, :]

    # linear interpolation
    out = (f00 * (1 - x)[:, None] * (1 - y)[:, None] + f01 *
           (1 - x)[:, None] * y[:, None] + f10 * x[:, None] *
           (1 - y)[:, None] + f11 * x[:, None] * y[:, None])

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

    return np.column_stack((x_points, y_points))


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
    # TODO: test the speed compared with give geometry
    unique_data = np.array([data[i, j] for i, j in unique])

    return pout, unique_data[inverse].reshape(N, 4, n3)


def arbitray_line(data, points, di: float = 1):
    """
    data : SegyNP or np.ndarray
    points : the points to interpolate a path
    di : step
    """
    p = interpolate_path(points, di)
    pout, pdata = extract_data(data, p)
    out = interp(pout, pdata)

    return out, p
