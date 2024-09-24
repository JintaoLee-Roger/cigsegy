"""
Test the scan function.

Test SEG-Y files:
1. rogan_2039_psdm_full_padded.sgy
    3D, istep = 2, IBM float, 1.6GB, h1001x769x663

2. nlog.nl/D/D15/Z3FIN1994A-2.sgy
    3D, istep=10, xstep=2, int16, 3.9G

3. B10: line missing

4. GR_08.sgy: descending
    istep=1, xstep=-1, IBM, 7.5G
"""

import numpy as np
import cigse