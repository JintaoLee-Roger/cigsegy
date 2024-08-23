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
    # In python, the shape is (n-inline, n-crossline, n-time) = (663, 769, 1001).

    # shape: (n-time, n-crossline, n-inline) = (1001, 769, 663)
    # sample interval: 4000, data format code: 4-bytes IBM floating-point
    # inline range: 360 - 1684, crossline range: 1764 - 2532
    # interval of inline: 35.0, interval of crossline: 17.5, time start: 0
    # inline field: 189, crossline field: 193
    # inline step: 2, crossline step: 1
    # Is regular file (no missing traces): false

    >>> cigsegy.metaInfo('fx.segy', iline=9, xline=21)
    # In python, the shape is (n-inline, n-crossline, n-time) = (1001, 120, 1500).

    # shape: (n-time, n-crossline, n-inline) = (1500, 120, 1001)
    # sample interval: 4000, data format code: 4-bytes IBM floating-point
    # inline range: 3 - 1003, crossline range: 1 - 120
    # interval of inline: 25.0, interval of crossline: 25.0, time start: 0
    # inline field: 9, crossline field: 21
    # inline step: 1, crossline step: 1
    # Is regular file (no missing traces): true


In some SEG-Y files, we cannot get useful information from textual header, i.e., 
don't know iline/xline/istep/xstep. You can just ignore them, and cigsegy 
will automatically guess the locations and steps of inline and crossline.

.. Note::

    For the previous version, you may need to set ``use_guess=True``. For the lasted version,
    you don't need do this, just ignore them.

.. code-block:: python

    >>> cigsegy.metaInfo('rogan.sgy', 189, 193) # ignore istep and xstep
    # In python, the shape is (n-inline, n-crossline, n-time) = (663, 769, 1001).

    # shape: (n-time, n-crossline, n-inline) = (1001, 769, 663)
    # sample interval: 4000, data format code: 4-bytes IBM floating-point
    # inline range: 360 - 1684, crossline range: 1764 - 2532
    # interval of inline: 35.0, interval of crossline: 17.5, time start: 0
    # inline field: 189, crossline field: 193
    # inline step: 2, crossline step: 1
    # Is regular file (no missing traces): false


To get the meta information in ``dict`` format, use ``cigsegy.tools.get_metaInfo``:

.. code-block:: python

    >>> meta = cigsegy.tools.get_metaInfo('rogan.sgy')
    >>> print(meta)
    # {'nt': 1001, 'nx': 769, 'ni': 663, 'trace_count': 380762, 
    # 'dt': 4000, 'dtype': '>4f-ibm', 'scalar': 1, 'i-interval': 35.01677322387695, 
    # 'x-interval': 17.499217987060547, 'start_time': 0, 'min-iline': 360, 
    # 'max-iline': 1684, 'min-xline': 1764, 'max-xline': 2532, 'isnormal': False, 
    # 'iline': 189, 'xline': 193, 'xloc': 73, 'yloc': 77, 'istep': 2, 'xstep': 1, 
    # 'fills': 0.0}


.. Note::

    You can use ``cigsegy.tools.trace_count('rogan.sgy')`` to get the trace number,
    and use ``cigsegy.tools.nt('rogan.sgy')`` to get the number of time samples for one trace.


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

If the SEG-Y file is unsorted, you can use ``cigsegy.scan_unsorted3D`` to scan the 
geometry of the file, and then use ``cigsegy.load_unsorted3D`` to read the data.

But, please note that ``cigsegy.scan_unsorted3D`` is slow, because it needs to scan the whole file.
Besides, ``iline`` and ``xline`` are required to be specified.

.. code-block:: python

    >>> geom = cigsegy.scan_unsorted3D('rogan.sgy', 189, 193) # must pass iline and xline
    >>> d = cigsegy.load_unsorted3D('rogan.sgy', geom)




Read 3D poststack data by ignoring header
==========================================

.. Note::

    This feature will be deseperated in the future version. 
    You can use ``cigsegy.collect('rogan.sgy').reshape(ni, nx, nt)`` to do the same thing.


If the header is broken and the shape and data format is already known, 
you can ignore header by specify the shape by using ``cigsegy.fromfile_ignore_header``:

.. code-block:: python

    # format: 1 for 4 bytes IBM float, 5 for 4 bytes IEEE float
    >>> d = cigsegy.fromfile_ignore_header('rogan.sgy', 663, 769, 1001, format=1)
    >>> d.shape 
    # (663, 769, 1001) # (n-inline, n-crossline, n-time)

    # tofile mode
    >>> cigsegy.tofile_ignore_header('rogan.segy', 'out.dat', 663, 769, 1001, format=1)


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

.. Note::

    Use ``cigsegy.SegyNP`` class is a more efficient way, which treats the SEG-Y file as a 3D/2D numpy array.
    Please see ``SegyNP`` for more details. (From version 1.1.7)


Use ``Pysegy`` class to read arbitrary sub-volumes, slices, traces.

.. code-block:: python

    import cigsegy
    from cigsegy import Pysegy
    import numpy
    # assume rogan.segy's iline/xline/istep/xstep is 189/193/2/1
    segy = Pysegy('rogan.segy')
    segy.setSteps(2, 1) # set inline and crossline steps
    segy.setInlineLocation(189)
    segy.setCrosslineLocation(193)
    segy.scan()
    iline89 = segy.read_inline_slice(89) # inline 89
    xline101 = segy.read_cross_slice(101) # xline 101
    time200 = segy.read_time_slice(200) # time 200
    trace28_93 = segy.read_trace(500, 300) # inline 500, xline 300
    # inline 10:100, xline 20:200, time 30:100
    subvol = segy.read(10, 100, 20, 200, 30, 100)
    segy.close_file()


.. note::
    If you want to read a trace with trace number index (range 0-N, 
    N is the total number of trace) rather than inline index and 
    crossline index, you can use ``cigsegy.collect`` function:

    .. code-block:: python

        # read the 500-th trace
        >>> cigsegy.collect('rogan.segy', 500)


Plot the slices, you will see:

.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/slice.png
    :alt: slices
    :align: center


Cut a sub SEG-Y
===============

Use ``Pysegy.cut`` to cut a sub-volume and keep its trace headers,
save as a sub SEG-Y file.

.. code-block:: python

    import cigsegy
    from cigsegy import Pysegy
    import numpy
    # assume rogan.segy's iline/xline/istep/xstep is 189/193/2/1
    segy = Pysegy('rogan.segy')
    segy.setSteps(2, 1) # set inline and crossline steps
    segy.setInlineLocation(189)
    segy.setCrosslineLocation(193)
    segy.scan()

    # volume: inline 10:100, xline 20:200, time 30:100
    segy.cut('out1.segy', 10, 100, 20, 200, 30, 100)

    # volume: inline 10:100, xline 20:200, time : (all)
    segy.cut('out2.segy', 10, 100, 20, 200)

    # volume: inline : (all), xline : (all), time 30:100
    segy.cut('out3.segy', 30, 100)

    segy.close_file()




Plot region map 
===============

Use ``cigsegy.plot.plot_region`` to plot the region where the segy file was located
(x/y axis is inline/crossline).

.. code-block:: python

    # loc: [iline, xline, istep, xstep]
    >>> cigsegy.plot.plot_region('rogan.sgy', loc=[9, 21, 2, 1])

You will see:

.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/rogan.png
    :alt: rogan
    :align: center

``rogan.sgy`` file is a **irregular** SEG-Y file which missing some traces.

If you want to plot the region in CDP X and CDP Y axis, set ``mode='cdpxy'``, and set 
``cdpxy_loc=[cdpx, cdpy]`` if nessesary.

.. code-block:: python

    # loc: [iline, xline, istep, xstep]
    >>> cigsegy.plot.plot_region('rogan.sgy', mode='cdpxy', loc=[9, 21, 2, 1], cdpxy_loc=[73, 7])


.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/roganxy.png
    :alt: roganxy
    :align: center
