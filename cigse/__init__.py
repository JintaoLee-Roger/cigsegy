# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
"""
**cigsegy** is a tool for exchanging data between **SEG-Y** format and 
**NumPy** array inside Python environment.

Core Features
--------------

- Fast (Implemented in c++)
- python wraping and **numpy** array supports
- dealing with normal and **irregular** SEG-Y volume [1]_.
- creating a SEG-Y file using the **existed header** of a SEG-Y


Project url: [https://github.com/JintaoLee-Roger/cigsegy](https://github.com/JintaoLee-Roger/cigsegy) \n
Documents: [https://cigsegy.readthedocs.io](https://cigsegy.readthedocs.io) \n
pypi: [https://pypi.org/project/cigsegy](https://pypi.org/project/cigsegy) 
"""


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

    def __init__(self, e, custom=''):
        if custom:
            self.exception = type(e)(f"{e.args[0]}\n\t{custom}", *e.args[1:])
        else:
            self.exception = e

    def __getattribute__(self, *args, **kwargs):
        if args[0] == "__class__":
            return None.__class__
        raise super().__getattribute__("exception")

    def __call__(self, *args, **kwargs):
        raise super().__getattribute__("exception")


from .deprecated import *
from .factories import *
from .segynp import SegyNP
from .cpp._CXX_SEGY import Pysegy
from . import plot
from . import tools
from . import utils
from . import interp
from . import transform
from . import constinfo
from . import createtool
from .createtool import SegyCreate
