class ExceptionWrapper:
    """
    Copy from `trimesh.exceptions.ExceptionWrapper`

    Create a dummy object which will raise an exception when attributes
    are accessed (i.e. when used as a module) or when called (i.e.
    when used like a function)

    For soft dependencies we want to survive failing to import but
    we would like to raise an appropriate error when the functionality is
    actually requested so the user gets an easily debuggable message.
    """

    def __init__(self, exception):
        self.exception = exception

    def __getattribute__(self, *args, **kwargs):
        if args[0] == "__class__":
            return None.__class__
        raise super().__getattribute__("exception")

    def __call__(self, *args, **kwargs):
        raise super().__getattribute__("exception")


from .cigsegy import (  # type: ignore
    Pysegy, fromfile, fromfile_without_scan, tofile, create_by_sharing_header,
    modify_trace_key, modify_bin_key)
from .tools import (create, textual_header, metaInfo, scan_prestack,
                    load_prestack3D, collect, scan_unsorted3D, load_unsorted3D)
from .utils import (progress_bar, get_trace_keys, get_trace_keys2)

try:
    from . import plot
except BaseException as e:
    plot = ExceptionWrapper(e)

from .segynp import SegyNP

from .deprecated import *

progress_bar()

__all__ = [
    "Pysegy", "fromfile", "fromfile_without_scan", "fromfile_ignore_header",
    "tofile", "tofile_ignore_header", "create", "collect", "textual_header",
    "metaInfo", "create_by_sharing_header", "get_trace_keys",
    "get_trace_keys2", "scan_prestack", "load_prestack3D", "scan_unsorted3D",
    "load_unsorted3D", "modify_bin_key", "modify_trace_key", "SegyNP"
]
