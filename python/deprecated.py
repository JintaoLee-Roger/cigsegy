import warnings
from .tools import collect


def fromfile_ignore_header(segy_name: str,
                           sizeZ: int,
                           sizeY: int,
                           sizeX: int,
                           format: int = 5):
    # TODO: to be removed in v1.1.9
    warnings.warn(
        "`fromfile_ignore_header` is deprecated and will be removed in a future version. Please use `collect` instead."
        " e.g., `d = cigsegy.collect(segy_name).reshape(ni, nx, nt)`",
        DeprecationWarning,
        stacklevel=2)

    return collect(str(segy_name)).reshape(sizeZ, sizeY, sizeX)


def tofile_ignore_header(segy_name: str,
                         out_name: str,
                         sizeX: int,
                         sizeY: int,
                         sizeZ: int,
                         format: int = 5):
    # TODO: to be removed in v1.1.9
    warnings.warn(
        "`tofile_ignore_header` is deprecated and will be removed in a future version. Please use `collect` instead."
        " e.g., `d = cigsegy.collect(segy_name).reshape(ni, nx, nt)`",
        DeprecationWarning,
        stacklevel=2)

    d =  collect(str(segy_name)).reshape(sizeZ, sizeY, sizeX)
    d.tofile(out_name)
