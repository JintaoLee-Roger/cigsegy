Create a SEG-Y file
###################

Create from numpy array (or binary file) and SEG-Y headers
==========================================================

There is often such a workflow:
    a. Display SEG-Y format data ``orig.segy`` in specialized software, such as Petrel.
    b. Use Python code to process this data and obtain new data ``afterprocess``, which is in NumPy array format
    c. To display this processed data in specialized software, it needs to be converted back to SEG-Y format and use the headers from the original data, i.e., using the NumPy array ``afterprocess`` and the header of ``orig.segy`` to create a new SEG-Y file ``out.segy``.

You can use ``cigsegy.create_by_sharing_header`` to create a new 
SEG-Y file ``out.segy`` whose headers are same as ``orig.segy``.

.. code-block:: python

    # read orig.segy to numpy array
    >>> orig = cigsegy.fromfile('orig.segy', 9, 21)
    >>> orig.shape 
    # (589, 762, 1001)

    # do some prosess
    >>> afterprocess = do_some_process(orig)
    >>> afterprocess.shape 
    # (589, 762, 1001)

    # assume the iline/xline/istep/xstep of **orig.segy** are 9/21/1/1
    >>> cigsegy.create_by_sharing_header('out.segy', 'orig.segy', afterprocess, \
        keylocs=[9, 21, 1, 1])

    # using binary file instead of numpy array
    # assume shape = (589, 762, 1001) = (n-inline, n-crossline, n-time)
    >>> cigsegy.create_by_sharing_header('out.segy', 'orig.segy', 'afterprocess.dat', \
        (589, 762, 1001), keylocs=[9, 21, 1, 1])


If ``afterprocess`` is a sub-volume, i.e., only process a part of ``orig.segy``, you 
can set ``start=(iline_start, xline_start, time_start)`` to create a new SEG-Y file.

For example:

.. code-block:: python

    # read orig.segy to numpy array
    >>> orig = cigsegy.fromfile('orig.segy', 9, 21)
    >>> orig.shape 
    # (589, 762, 1001)

    # cut orig, only process a part of orig array
    # offset is (100, 400, 300)
    >>> toprocess = orig[100:520, 400:700, 300:900]
    >>> afterprocess = do_some_process(toprocess)
    >>> afterprocess.shape 
    # (420, 300, 600)

    # create 'out.segy' using afterprocess (data) and a part
    # of header in 'orig.segy'
    >>> cigsegy.create_by_sharing_header('out.segy', 'orig.segy', afterprocess, \
        keylocs=[9, 21, 1, 1], start=(100, 400, 300))


.. Note::

    For all ``create``-like functions, you can speacify the ``textual`` to 
    custom textual header, e.g.,

    .. code-block:: python
        
        # you can special the first 12 lines of textual header
        >>> textual12 = "xxxxx"
        >>> cigsegy.create_by_sharing_header('out.segy', 'orig.segy', afterprocess, \
            keylocs=[9, 21, 1, 1], start=(100, 400, 300), textual=my_info)
        >>> cigsegy.textual_header('out.segy')
        # xxx



Create from numpy array and some parameters
===========================================

.. code-block:: python

    # d is a numpy array, d.shape == (n-inlines, n-crosslines, n-time)
    >>> cigsegy.create('out.segy', d, format=5, start_time=0, iline_interval=15, ...)


