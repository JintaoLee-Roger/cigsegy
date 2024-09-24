"""
test the following functions:

1. create_by_sharing_header
2. cut

"""

import cigse 
from cigse import Pysegy
import numpy as np

######### test create_by_sharing_header #########

# full fill
def test_cbsh1():
    cigse.create_by_sharing_header('out.segy', 'rogan_2039_psdm_full_padded.sgy', )


# with offset
def test_cbsh2():
    cigse.create_by_sharing_header('out.segy', 'rogan_2039_psdm_full_padded.sgy', )

######### test cut #########

def test_cut():
    segy = Pysegy('rogan_2039_psdm_full_padded.sgy')
    segy.cut('cut.segy', [0, 100, 0, 100, 0, 100])
    segy.close()
    