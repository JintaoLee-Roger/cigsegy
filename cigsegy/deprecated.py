# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.


def fromfile_ignore_header(*args, **kwargs):
    raise RuntimeError(
        "`fromfile_ignore_header` is deprecated and removed. Please use `collect` instead."
        " e.g., `d = cigsegy.collect(segy_name).reshape(ni, nx, nt)`")


def tofile_ignore_header(*args, **kwargs):
    raise RuntimeError(
        "`tofile_ignore_header` is deprecated and removed. Please use `collect` instead."
        " e.g., `d = cigsegy.collect(segy_name).reshape(ni, nx, nt)`")
