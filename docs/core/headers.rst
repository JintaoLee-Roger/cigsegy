Headers
##############

There are 3 kind of headers in a SEG-Y file: Textual Header, Binary Header, Trace Header.
For most SEG-Y files, each file contains 1 Textual Header, 1 Binary Header and 
several Trace Headers (one trace has one trace header).

.. _textual_header:

Textual Header
==============

The first **3200**-bytes, Textual File Header record contains 40 lines 
of textual information, providing a human-readable description of the 
seismic data in the SEG-Y file. It originally was encoded in the EBCDIC 
character set but ASCII is now allowed for all Textual File Headers 
(since `v1 <https://library.seg.org/pb-assets/technical-standards/seg_y_rev1-1686080991247.pdf>`_).

The textual header contains some key information descripted the SEG-Y file 
and some key information to access the file, such as the inline location,
crossline location, cdp x, cdp y ...

Use ``cigsegy.textual_header`` to print this 3200 bytes header:

.. code-block:: python

    >>> cigsegy.textual_header('rogan.sgy')
    # C01 CLIENT: STUART PETROLEUM LTD    AREA:COOPER BASIN   SOUTH AUSTRALIA
    # ...
    # C06 INLINE RANGE: 360 - 1684(2)  CROSSLINE RANGE 1764 - 2532(1)
    # C07 -------------PROCESSING FLOW---------------
    # ...
    # C35  DESC                   BYTE LOCATION       FORMAT 
    # C36  3D INLINE NUMBER        9- 12             32 BIT INTEGER 
    # C37  3D CROSSLINE  NUMBER   21- 24             32 BIT INTEGER 
    # C38  CDP_X                  73- 76             32 BIT INTEGER 
    # C39  CDP_Y                  77- 80             32 BIT INTEGER 
    # C40  

You can get some key information to read the SEG-Y file, such as inline location 
is 9 (C36), crossline location is 21 (C37), X location is 73 (C38), Y location 
is 77 (C39), inline step is 2 (C06), crossline step is 1 (C06).

You can input ``coding='a'`` or ``coding='e'`` to force the 3200 bytes as ASCII or EBCDIC 
coding, such as :

.. code-block:: python

    # 'e' for EBCDIC, 'a' for ASCII
    >>> cigsegy.textual_header('rogan.sgy', coding='e')


Binary Header
=============

The **400** bytes Binary File Header record contains binary values relevant to the whole SEG-Y file.

Use ``cigsegy.tools.read_header`` to get all values of the binary header.

.. code-block:: python

    # just print all keys and their values
    >>> cigsegy.tools.read_header('rogan.sgy', type='bh')
    # Bytes 3201 - 3204: 0        -- Job ID
    # ...
    # Bytes 3217 - 3218: 4000     -- Sample interval(dt)
    # Bytes 3219 - 3220: 0        -- dt of original
    # Bytes 3221 - 3222: 1001     -- N samples per traces(ns)
    # Bytes 3223 - 3224: 0        -- ns of orignal
    # Bytes 3225 - 3226: 1        -- Data sample format code (1-IBM, 5-IEEE)
    # ...

    # get a dict: key: (desc, value), key is the location
    >>> bh = cigsegy.tools.read_header('rogan.sgy', type='bh', printstr=False)
    >>> print(bh)
    # {1: ('Job ID', 0), 5: ('Line number', 0), 9: ('Reel Number', 0), 
    # 13: ('N traces per ensembel', 4734), 15: ('N auxiliary traces per ensembel', 0), 
    # 17: ('Sample interval(dt)', 4000), 19: ('dt of original', 0), 
    # 21: ('N samples per traces(ns)', 1001)
    # ...


For more detail description of each location, pelease read `SEG-Y v2 (2017) <https://library.seg.org/pb-assets/technical-standards/seg_y_rev2_0-mar2017-1686080998003.pdf>`_


Trace Header
============

The SEG-Y trace header (240 bytes for each trace header) contains trace attributes, 
most of which are defined with two-byte or four-byte, two's complement integers.

Use ``cigsegy.tools.read_header`` to get all values of the n-th trace header.

.. code-block:: python

    # print 10-th trace header
    >>> cigsegy.tools.read_header('rogan.sgy', type='th', n=10)
    # Bytes 1    - 4   : 11       -- Trace sequence number within line
    # Bytes 5    - 8   : 0        -- Trace sequence number within SEG-Y file
    # Bytes 9    - 12  : 360      -- Original field record number
    # Bytes 13   - 16  : 0        -- Trace number within the original field record

    # like binary header, set printstr=False to get a dict
    >>> th = cigsegy.tools.read_header('rogan.sgy', type='th', n=10, printstr=False)


To get multiple values (N traces, N values) in a special location, use ``cigsegy.get_trace_keys``,
for example, get the values at location 21 for the first 5 traces:

.. code-block:: python

    >>> d = cigsegy.get_trace_keys('rogan.sgy', keyloc=21, beg=0, end=5)
    >>> d
    # array([2101, 2102, 2103, 2104, 2105])


Some case:

- set ``beg`` as a negative number, e.g., ``beg=-1`` to read the values of all the traces, just like: ``value[:]``;
- set ``end`` as a negative number, e.g., ``end=-1`` to read the values of the traces from ``beg`` to ``trace_count``, just like: ``value[beg:]``;
- set ``end=0`` or ignore it to read just one value of the trace from the ``beg``-th trace, just like: ``value[beg]``.

.. code-block:: python

    >>> d = cigsegy.get_trace_keys('rogan.sgy', keyloc=21, beg=-1)
    >>> d.shape
    # (380762, ) all traces

    >>> d = cigsegy.get_trace_keys('rogan.sgy', keyloc=21, beg=10000, end=-1)
    >>> d.shape 
    # (370762,), 10000 - trace_count

    >>> d = cigsegy.get_trace_keys('rogan.sgy', keyloc=21, beg=100)
    >>> d
    # 2201, the 100-th trace, just a int number


Set ``force`` as the length of the key, force to read the value even if the ``keyloc`` 
not in the standard SEG-Y trace header's keys.

.. code-block:: python

    # 221 is not a key in the standard SEG-Y trace header
    # read 221-224 4 bytes as a int
    >>> cigsegy.get_trace_keys('mx.sgy', keyloc=221, beg=100, force=4)
    # 1986

    # if read a float32, set force=5, read float64, set force=9



Set ``keyloc`` as a ``List`` or 1D ``np.ndarray`` to get keys in multiple locations.

.. code-block:: python

    >>> d = cigsegy.get_trace_keys('rogan.sgy', keyloc=[9, 21], beg=0, end=10)
    >>> d.shape
    # (10, 2)


Use ``cigsegy.plot.plot_trace_keys`` to plot the keys change along trace numbers:

.. code-block:: python

    >>> cigsegy.plot.plot_trace_keys('rogan.sgy', keyloc=21, beg=0, end=1000)


.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/ptk.png
    :alt: plot_trace_keys
    :align: center


Use ``cigsegy.plot.plot_trace_ix`` to plot the inline and crossline.

.. code-block:: python

    >>> cigsegy.plot.plot_trace_ix('rogan.sgy', iline=9, xline=21, beg=1000, end=2000, figsize=(8, 4))


.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/ptix.png
    :alt: plot_trace_keys
    :align: center


Use ``cigsegy.plot.plot_trace_ixo`` to plot the inline, crossline and offset

.. code-block:: python

    >>> cigsegy.plot.plot_trace_ixo('3Dgather.sgy', iline=5, xline=9, offset=37, beg=1000, end=1100, figsize=(10, 4))


.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/ptixo.png
    :alt: plot_trace_keys
    :align: center