About SEG-Y
###########

SEG-Y is the common exchange format for seismic data.  A practical reader or
writer must handle two things at the same time:

- the byte-level SEG-Y structure;
- the survey geometry encoded in trace headers.

``cigsegy`` focuses on making those two layers accessible from Python while
using C++ for the heavy file I/O and sample encoding work.


SEG-Y Revisions
===============

Useful references:

- `SEG-Y rev 0 (1975) <https://library.seg.org/pb-assets/technical-standards/seg_y_rev0-1686080980707.pdf>`_
- `SEG-Y rev 1 (2002) <https://library.seg.org/pb-assets/technical-standards/seg_y_rev1-1686080991247.pdf>`_
- `SEG-Y rev 2.0 (2017) <https://library.seg.org/pb-assets/technical-standards/seg_y_rev2_0-mar2017-1686080998003.pdf>`_

Many real files mix conventions from different revisions.  ``cigsegy`` reads
the binary header and trace headers directly, so you can override geometry byte
locations when a file does not follow the expected convention.


File Structure
==============

The common layout is:

1. 3200-byte textual header;
2. 400-byte binary header;
3. optional extended textual headers;
4. repeated trace records:

   - 240-byte trace header;
   - sample bytes for that trace;

5. optional data trailer.

``SegyWriter`` can write textual headers, binary headers, extended textual
headers, trace headers, sample bytes, and an optional data trailer.  Most common
reading workflows still assume one textual header, one binary header, and
file-order trace records.


Geometry
========

SEG-Y does not have a universal array model.  A 3D post-stack cube becomes an
array only after you decide which trace header fields represent inline and
crossline.  A 4D/pre-stack volume also needs an offset field.

Common byte locations are:

- inline: 189 or 9;
- crossline: 193 or 21;
- offset: 37;
- X coordinate: 181 or 73;
- Y coordinate: 185 or 77.

These are conventions, not guarantees.  Use ``textual_header``,
``metaInfo``, ``get_trace_keys``, or plotting helpers to confirm a file.


Sample Formats
==============

The binary header field at bytes 3225-3226 stores the SEG-Y sample format code.
Common values include:

- ``1``: 4-byte IBM floating point;
- ``2``: 4-byte signed integer;
- ``3``: 2-byte signed integer;
- ``5``: 4-byte IEEE floating point;
- ``8``: 1-byte signed integer.

``SegyWriter.write_trace_block`` accepts float32-compatible data and lets the
C++ backend encode it according to the selected sample format.
``write_raw_trace_block`` bypasses conversion and copies already encoded sample
bytes exactly.


Post-Stack and Pre-Stack
========================

Post-stack data is usually represented as:

.. code-block:: text

   (n_inline, n_xline, n_sample)

Pre-stack or gather data often becomes:

.. code-block:: text

   (n_inline, n_xline, n_offset, n_sample)

2D lines can be treated as:

.. code-block:: text

   (n_trace, n_sample)
