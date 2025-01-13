Read Prestack data
###################


.. note::
    cigsegy-**v1.1.7** doesn't support read a **unsorted** prestack gather in a convenient
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

Use ``cigsegy.metaInfo`` to scan the prestack file.

.. code-block:: python

    >>> cigsegy.metaInfo('3Dgather.sgy', iline=9, xline=21, offset=37)


Use ``cigsegy.fromfile`` to load 3D prestack file.

.. code-block:: python

    # geom can be obtained by `scan_prestack`
    >>> d = cigsegy.fromfile('3Dgather.sgy', iline=9, xline=21, offset=37)
    >>> d.shape
    # (41, 482, 61, 1501) = (ni, nx, no, nt), no means number of offset


Load 2D prestack SEG-Y (3D array)
----------------------------------

Load 2D prestack SEG-Y using the same approach as for 3D poststack data, i.e., 
use ``fromfile``.

Set the cdp/line location as ``iline``, and set offset location as ``xline``.

.. code-block:: python

    >>> cigsegy.metaInfo('2Dgather.sgy', iline=9, xline=37)

    >>> d = cigsegy.fromfile('2Dgather.segy', iline=9, xline=37)
    >>> d.shape
    # (1001, 120, 1500) = (n-line, no, nt), no is number of offset