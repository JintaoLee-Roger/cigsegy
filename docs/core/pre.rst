Read Prestack data
###################


.. note::
    cigsegy-**v1.1.6** doesn't support read a **unsorted** prestack gather in a convenient
    way. Prestack seismic gather will be supported in the next version of cigsegy.
    
    However, this **doesn't** mean that you cannot use cigsegy to read the **unsorted** prestack data.
    You can use ``get_trace_keys()`` function to analyze SEG-Y file's
    trace headers and ``collect()`` to access data traces.


In SEG-Y file, the term CDP (common depth point) as used in this 
document is used as a synonym for the term CMP (common midpoint).

The differences between several gathers are shown in the following figure.

.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/gather.png 
    :alt: gather
    :align: center

In general, there are three keys to evaluate 3D prestack SEG-Y file: ``inline``, ``xline`` or ``crossline``, and ``offset``.


Load 3D prestack SEG-Y (4D array)
----------------------------------

Use ``cigsegy.scan_prestack`` to scan the prestack file.

.. code-block:: python

    >>> geom = cigsegy.scan_prestack('3Dgather.sgy', iline=5, xline=9, offset=37)
    >>> geom
    # {'shape': [41, 482, 61, 1501],
    # 'iline': {'min_iline': 1000, 'max_iline': 1040, 'istep': 1},
    # 'xline': {'min_xline': 1000, 'max_xline': 1481, 'xstep': 1},
    # 'offset': {'min_offset': 175, 'max_offset': 3175, 'ostep': 50}}


Use ``cigsegy.load_prestack3D`` to load 3D prestack file.

.. code-block:: python

    # geom can be obtained by `scan_prestack`
    >>> d = cigsegy.load_prestack3D('3Dgather.sgy', iline=5, xline=9, offset=37, geom=geom)
    >>> d.shape
    # (41, 482, 61, 1501) = (ni, nx, no, nt), no means number of offset

    # ignore geom
    >>> d = cigsegy.load_prestack3D('3Dgather.sgy', iline=5, xline=9, offset=37)


.. note::
    Unlike 3D post-stack seismic, ``scan_prestack`` scans the whole file instead of 
    a part of file, so ``scan_prestack`` may be slow.


Load 2D prestack SEG-Y (3D array)
----------------------------------

Load 2D prestack SEG-Y using the same approach as for 3D poststack data, i.e., 
use ``fromfile``.

Set the cdp/line location as ``iline``, and set offset location as ``xline``.

.. code-block:: python

    # get step of offset, which is 25
    >>> np.diff(cigsegy.get_trace_keys('2Dgather.segy', 37, 0, 10))
    # array([25, 25, 25, 25, 25, 25, 25, 25, 25])

    >>> d = cigsegy.fromfile('2Dgather.segy', iline=9, xline=37, istep=1, xstep=25)
    >>> d.shape
    # (1001, 120, 1500) = (n-line, no, nt), no is number of offset