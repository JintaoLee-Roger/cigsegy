Headers and Metadata
####################

SEG-Y files usually contain:

- one 3200-byte textual header;
- one 400-byte binary header;
- one 240-byte trace header before each trace;
- encoded sample bytes after each trace header.

``cigsegy`` keeps these concepts visible because reliable geometry depends on
the trace header fields used for inline, crossline, offset, and coordinates.


Textual Header
==============

Use ``textual_header`` to inspect the first 3200 bytes.

.. code-block:: python

   import cigsegy

   cigsegy.textual_header("input.sgy")
   text = cigsegy.textual_header("input.sgy", printtext=False)

Set ``coding`` to force decoding:

.. code-block:: python

   cigsegy.textual_header("input.sgy", coding="a")  # ASCII
   cigsegy.textual_header("input.sgy", coding="e")  # EBCDIC

Textual headers often describe byte locations for inline, crossline, CDP X, and
CDP Y, but they are not always reliable.  The usual workflow is to let
``cigsegy`` guess first, then pass explicit locations only if the result needs
correction.


Metadata Scan
=============

Use ``metaInfo`` for a readable summary:

.. code-block:: python

   cigsegy.metaInfo("input.sgy")

If the inferred geometry is wrong or ambiguous:

.. code-block:: python

   cigsegy.metaInfo("input.sgy", iline=189, xline=193)

Use ``tools.get_metaInfo`` when code needs a dictionary:

.. code-block:: python

   meta = cigsegy.tools.get_metaInfo("input.sgy")
   print(meta["ni"], meta["nx"], meta["nt"], meta["dt"])

Important keys include:

- ``ni``, ``nx``, ``no``, ``nt``: logical shape;
- ``dt``: sample interval in microseconds;
- ``dformat``: SEG-Y sample format code;
- ``start_iline``, ``end_iline``, ``start_xline``, ``end_xline``;
- ``iline``, ``xline``, ``offset``, ``xloc``, ``yloc``: byte locations.


Binary Header
=============

Use ``tools.read_header(..., type="bh")`` to inspect all known binary header
fields.

.. code-block:: python

   cigsegy.tools.read_header("input.sgy", type="bh")

   binary = cigsegy.tools.read_header("input.sgy", type="bh", printstr=False)
   print(binary[17])  # sample interval
   print(binary[21])  # samples per trace
   print(binary[25])  # sample format code

When you need raw bytes, use the C++ object directly:

.. code-block:: python

   segy = cigsegy.Pysegy("input.sgy")
   try:
       binary_bytes = segy.get_binary_header()
   finally:
       segy.close()


Trace Header
============

Read a decoded trace header:

.. code-block:: python

   cigsegy.tools.read_header("input.sgy", type="th", n=0)

   trace0 = cigsegy.tools.read_header("input.sgy", type="th", n=0, printstr=False)
   print(trace0[189])  # inline, if byte 189 is used by this file

Read selected trace header keys for many traces:

.. code-block:: python

   ix = cigsegy.get_trace_keys("input.sgy", keyloc=[189, 193], beg=0, end=1000)

Use ``indices`` for non-contiguous trace access:

.. code-block:: python

   import numpy as np

   indices = np.array([0, 10, 20], dtype=np.int32)
   keys = cigsegy.get_trace_keys("input.sgy", keyloc=[189, 193], indices=indices)

Use ``force`` for non-standard byte locations:

.. code-block:: python

   value = cigsegy.get_trace_keys("input.sgy", keyloc=221, beg=0, force=4)


Geometry Fields
===============

For a regular 3D post-stack file, common locations are:

- inline: 189 or 9;
- crossline: 193 or 21;
- X coordinate: 181 or 73;
- Y coordinate: 185 or 77.

For a 4D/pre-stack file, add an offset field, often byte 37.

These are conventions, not guarantees.  Always confirm with ``textual_header``,
``plot_trace_keys``, or a metadata scan when the file is unfamiliar.
