"""
test read functions

1. fromfile
2. ...

involving the following ablities:
1. data types: IEEE float, IBM float, (int16, int32)
2. different arbitry geometry
    1. normal one
    2. irregular one
    3. missing lines
    4. 2D/3D/4D
    5. one dimension is 1
    6. ...
"""

import numpy as np
import cigse
from pathlib import Path

root = Path('/disk/d1/share_vm/cigsegy_data')

def test_rogan():
    d1 = cigse.fromfile(root / 'rogan_2039_psdm_full_padded.sgy')
    d2 = np.fromfile(root / 'rogan_h1001x769x663.dat', np.float32).reshape(663, 769, 1001)
    assert np.allclose(d1, d2)

def test_d15():
    d1 = cigse.fromfile(root / 'Z3FIN1994A-2.sgy')
    d2 = np.fromfile(root / 'Z3FIN1994A.dat', np.float32).reshape()
    assert np.allclose(d1, d2)

def test_gr08():
    d1 = cigse.fromfile(root / 'GR_08.sgy')
    d2 = np.fromfile(root / 'GR_08.dat', np.float32).reshape()
    assert np.allclose(d1, d2)

