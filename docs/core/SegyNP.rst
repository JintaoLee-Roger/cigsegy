SegyNP class
###################

This class can help you treat SEG-Y files as 3D/2D numpy arrays, without reading the whole file into memory.
This means the reading process only occurs when you access the data, such as a part of data.


3D array
==========

If the SEG-Y file is a 3D post-stack seismic file:

.. code-block:: python

    >>> from cigsegy import SegyNP
    >>> d = SegyNP('3Dpoststack.sgy', iline=189, xline=193)
    >>> d.shape # (ni, nx, nt), use as a numpy array, 3D geometry
    >> sx = d[100] # the 100-th inline profile
    >> sx = d[100:200] # return a 3D array with shape (100, nx, nt)
    >> sx = d[:, 200, :] # the 200-th crossline profile
    >> sx = d[:, :, 100] # the 100-th time slice, note, it may be slow if the file is large
    >> sx.min(), sx.max() 
    # get the min and max value, but they are evaluated from a part of data, 
    # so they may not be the real min and max value
    >> sx.trace_cout # get the number of traces for the file
    >> sx.close() # close the file


2D array
==========

If you don't want to create a 3D geometry, just treat the SEG-Y file as a collection of 
1D traces, i.e., a 2D array, you can use the following code:

.. code-block:: python

    >>> from cigsegy import SegyNP
    >>> d = SegyNP('2Dseismic.sgy', as_2d=True)
    >>> d.shape # (trace_count, nt), use as a numpy array, collection of 1D traces
    >> sx = d[100] # the 100-th trace
    >> sx = d[100:200] # return a 2D array with shape (100, nt)
    >> sx.min(), sx.max() 
    # get the min and max value, but they are evaluated from a part of data, 
    # so they may not be the real min and max value
    >> sx.trace_cout # the number of traces for the file