Read Poststack data
###################


Meta Information
================

Use ``cigsegy.metaInfo`` to scan the SEG-Y file and obtain some meta information
of the file, such as shape, dt, data format ...

From textual header, we may get some key information, such as inline/crossline 
location, inline/crossline step, cdp x/y location, ... See :ref:`textual_header`.

.. code-block:: python

    # inline/xline/istep/xstep/xloc/yloc = 189/193/2/1/73/77
    >>> cigsegy.metaInfo('rogan.sgy', iline=189, xline=193, istep=2, xstep=1, xloc=73, yloc=77)
    # shape: (n-inline, n-crossline, n-time) = (663, 769, 1001)
    # N traces: 380762
    # interval: di(iline) = 35.00 ft, dx(xline) = 17.50 ft, dt = 4 ms
    # range: inline: 360 - 1684, crossline: 1764 - 2532, t: 0 - 4000.0 ms
    # trace sorting code: Unknown
    # scalar: 1, data format: 4-byte IBM floating-point
    # (key info) iline: 189, xline: 193, xloc:  73, yloc:  77
    #         istep:   2, xstep:   1

    >>> cigsegy.metaInfo('rogan.sgy')
    # shape: (n-inline, n-crossline, n-time) = (663, 769, 1001)
    # N traces: 380762
    # interval: di(iline) = 35.00 ft, dx(xline) = 17.50 ft, dt = 4 ms
    # range: inline: 360 - 1684, crossline: 1764 - 2532, t: 0 - 4000.0 ms
    # trace sorting code: Unknown
    # scalar: 1, data format: 4-byte IBM floating-point
    # (key info) iline: 189, xline: 193, xloc:  73, yloc:  77
    #         istep:   2, xstep:   1


In some SEG-Y files, we cannot get useful information from textual header, i.e., 
don't know iline/xline/istep/xstep. You can just ignore them, and cigsegy 
will automatically guess the locations and steps of inline and crossline.

.. Note::

    For the previous version, you may need to set ``use_guess=True``. For the lasted version,
    you don't need do this, just ignore them.

.. code-block:: python

    >>> cigsegy.metaInfo('rogan.sgy', 189, 193) # ignore istep and xstep
    # shape: (n-inline, n-crossline, n-time) = (663, 769, 1001)
    # N traces: 380762
    # interval: di(iline) = 35.00 ft, dx(xline) = 17.50 ft, dt = 4 ms
    # range: inline: 360 - 1684, crossline: 1764 - 2532, t: 0 - 4000.0 ms
    # trace sorting code: Unknown
    # scalar: 1, data format: 4-byte IBM floating-point
    # (key info) iline: 189, xline: 193, xloc:  73, yloc:  77
    #         istep:   2, xstep:   1


To get the meta information in ``dict`` format, use ``cigsegy.tools.get_metaInfo``:

.. code-block:: python

    >>> meta = cigsegy.tools.get_metaInfo('rogan.sgy')
    >>> print(meta)
    # {'iline': 189, 'istep': 2, 'offset': 37, 'ostep': 1, 'xline': 193, 
    # 'xloc': 73, 'xstep': 1, 'yloc': 77, 'dformat': 1, 'di': 34.99929428100586, 
    # 'dt': 4000, 'dx': 17.499540328979492, 'end_iline': 1684, 'end_offset': 0, 
    # 'end_xline': 2532, 'esize': 4, 'fillNoValue': 0.0, 'ndim': 3, 'ni': 663, 
    # 'no': 1, 'nt': 1001, 'ntrace': 380762, 'nx': 769, 'scalar': 1, 
    # 'start_iline': 360, 'start_offset': 0, 'start_time': 0, 'start_xline': 1764, 
    # 'trace_sorting_code': 0, 'tracesize': 4244, 'unit': 'ft'}


.. Note::

    You can use ``cigsegy.tools.trace_count('rogan.sgy')`` to get the trace number.


Read 3D poststack data
======================

Use ``cigsegy.fromfile`` to read data as ``numpy.ndarray``:

.. code-block:: python

    >>> d = cigsegy.fromfile('rogan.sgy', iline=9, xline=21, istep=2, xstep=1)
    >>> d.shape 
    # (663, 769, 1001) # (n-inline, n-crossline, n-time)


Use ``cigsegy.tofile`` to convert SEG-Y file to a binary file (without any headers).
When using ``cigsegy.tofile()``, you **don't** have to worry about 
running out of memory. Therefore, this function is very useful when 
dealing with **huge** files.

.. code-block:: python

    >>> cigsegy.tofile('rogan.sgy', 'out.dat', iline=9, xline=21, istep=2, xstep=1)




Read unsorted 3D poststack data
==================================

If the SEG-Y file is unsorted, you can use ``cigsegy.tools.full_scan`` to scan the 
geometry of the file, and then use ``cigsegy.tools.load_by_geom`` to read the data.

But, please note that ``cigsegy.tools.full_scan`` is slow, because it needs to scan the whole file.
Besides, ``iline`` and ``xline`` are required to be specified.

.. code-block:: python

    >>> geom = cigsegy.tools.full_scan('rogan.sgy', 189, 193) # must pass iline and xline
    >>> d = cigsegy.tools.load_by_geom('rogan.sgy', geom)




Read 3D poststack data by ignoring header
==========================================

You can use ``cigsegy.collect('rogan.sgy').reshape(ni, nx, nt)`` to do it.



Use plot tools you will see like:

.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/rogan3d.png
    :alt: rogan3d
    :align: center


Read 2D poststack data
======================

Use ``cigsegy.collect`` to read all traces as a 2D array:

.. code-block:: python

    >>> d = cigsegy.collect('L03_MIG_CB.sgy')
    >>> d.shape
    # (5038, 6000) # (n-traces, n-time)


.. note::

    ``collect`` function can read traces from ``beg`` to ``end``.

    - default, collect all traces

    .. code-block:: python

        >>> d = cigsegy.collect('L03_MIG_CB.sgy')
        >>> d.shape
        # (5038, 6000) # (n-traces, n-time)

    - set ``beg`` and ``end`` (``0 <= beg < end <= trace_count``) to collect traces from ``beg`` to ``end``, ``end`` is not included.

    .. code-block:: python

        >>> d = cigsegy.collect('L03_MIG_CB.sgy', 10, 100)
        >>> d.shape
        # (90, 6000)


    - read one trace with a trace index, this is equivalent to setting ``beg=index`` and ``end=0``.

    .. code-block:: python

        >>> d = cigsegy.collect('L03_MIG_CB.sgy', 100)
        >>> d.shape
        # (1, 6000), the 100-th trace

        # equivalent to
        >>> d = cigsegy.collect('L03_MIG_CB.sgy', beg=100, end=0)


    - set ``end=-1`` to collect traces from ``beg`` to ``trace_count``.

    .. code-block:: python

        >>> d = cigsegy.collect('L03_MIG_CB.sgy', 1000, -1)
        >>> d.shape
        # (4038, 6000), range like [1000:trace_count]


Arbitrary slicing and extration
===============================

Use ``cigsegy.SegyNP`` class is a more efficient way, which treats the SEG-Y file as a 3D/2D numpy array.
Please see ``SegyNP`` for more details. (From version 1.1.7)


Plot the slices, you will see:

.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/slice.png
    :alt: slices
    :align: center


Cut a sub SEG-Y
===============

Use ``cigsegy.SegyNP`` class is a more efficient way, which treats the SEG-Y file as a 3D/2D numpy array.


Plot region map 
===============

Use ``cigsegy.plot.plot_region`` to plot the region where the segy file was located
(x/y axis is inline/crossline).

.. code-block:: python

    # loc: [iline, xline, istep, xstep]
    >>> cigsegy.plot.plot_region('rogan.sgy', iline=9, xline=21, xxx)

You will see:

.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/rogan.png
    :alt: rogan
    :align: center

``rogan.sgy`` file is a **irregular** SEG-Y file which missing some traces.

If you want to plot the region in CDP X and CDP Y axis, set ``mode='cdpxy'``, and set 
``xloc=cdpx, yloc=cdpy`` if nessesary.

.. code-block:: python

    # loc: [iline, xline, istep, xstep]
    >>> cigsegy.plot.plot_region('rogan.sgy', mode='cdpxy', iline=9, xline=21, xloc=73, yloc=77, xxx)


.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/roganxy.png
    :alt: roganxy
    :align: center
