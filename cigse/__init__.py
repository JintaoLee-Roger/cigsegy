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
