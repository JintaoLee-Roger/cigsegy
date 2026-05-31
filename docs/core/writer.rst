Write with SegyWriter
#####################

``cigsegy.SegyWriter`` is the recommended writing API for new code.  It has
three modes:

- ``from_template``: reuse headers from an existing SEG-Y;
- ``from_headers``: write file-order trace headers that you already have;
- ``create``: generate regular headers from shape and geometry parameters.

Configure a builder first, then open it:

.. code-block:: python

   b = cigsegy.SegyWriter.from_template("input.sgy", "out.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(data)


Template Mode
=============

Use template mode when the output trace headers should come from an existing
SEG-Y file.  If ``keylocs`` is omitted, ``cigsegy`` tries to infer the geometry
from the template trace headers.  Use explicit ``keylocs`` only when the
template is ambiguous or uses non-standard header fields.

.. code-block:: python

   b.keylocs(iline=189, xline=193)


Whole Volume
------------

.. code-block:: python

   b = cigsegy.SegyWriter.from_template("input.sgy", "processed.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(processed)

Without ``select``, this path delegates to the stable
``create_by_sharing_header`` implementation.


Continuous Sub-Volume
---------------------

``start`` is the zero-based logical start in the template volume.

.. code-block:: python

   sub = processed[100:300, 40:200, 0:800]

   b = cigsegy.SegyWriter.from_template("input.sgy", "sub.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(sub, start=(100, 40, 0))


Spatial Thinning
----------------

Use ``select`` when the output should copy headers from a non-contiguous subset
of existing traces.  This does not create new spatial positions.

.. code-block:: python

   thin = data[::2, ::3, :]

   b = cigsegy.SegyWriter.from_template("input.sgy", "thin.sgy")
   b.select(iline=slice(None, None, 2), xline=slice(None, None, 3))
   b.overwrite(True)

   with b.open() as w:
       w.write(thin)

Equivalent shorthand:

.. code-block:: python

   b.spatial_stride(iline=2, xline=3)

If any selected spatial trace cannot be mapped back to a real template trace,
the writer raises an error instead of inventing headers.


Spatial Thinning and Time Downsampling
--------------------------------------

When ``sample`` has a constant step and ``sample_interval_us`` is not set
explicitly, the output sample interval is derived from the template interval.

.. code-block:: python

   down = data[::2, ::3, ::4]

   b = cigsegy.SegyWriter.from_template("input_1ms.sgy", "thin_4ms.sgy")
   b.select(
       iline=slice(None, None, 2),
       xline=slice(None, None, 3),
       sample=slice(None, None, 4),
   )
   b.overwrite(True)

   with b.open() as w:
       w.write(down)

Here the output ``dt`` becomes ``input_dt * 4``.


Time Super-Resolution
---------------------

Time super-resolution changes the sample axis but not the spatial trace
locations.  Leave ``sample`` unset and set the new interval explicitly.

.. code-block:: python

   b = cigsegy.SegyWriter.from_template("input_2ms.sgy", "output_1ms.sgy")
   b.strict(False).sample_interval_us(1000).overwrite(True)

   with b.open() as w:
       w.write(super_res_data)

For template writes without ``select``, ``sample_interval_us`` is passed to the
legacy path as ``dt_new``.  Use ``strict(False)`` when the output sample count or
sample interval is not a strict sub-window of the template.


Binary Sample File
------------------

.. code-block:: python

   b = cigsegy.SegyWriter.from_template("input.sgy", "out.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write_file("processed.dat", shape=(589, 762, 1001))

``.npy`` files are loaded with ``numpy.load(..., mmap_mode="r")``.  Raw binary
files are interpreted as float32.


Stored-Header Mode
==================

Use ``from_headers`` when textual, binary, and trace headers already exist in
file order.  This is the natural mode for chunked containers.

.. code-block:: python

   b = cigsegy.SegyWriter.from_headers("out.sgy")
   b.textual(textual)                  # 3200 bytes
   b.binary(binary)                    # 400 bytes
   b.sample_format(1)                  # IBM float
   b.sample_count(1001)
   b.overwrite(True)

   with b.open() as w:
       for trace_headers, samples in blocks:
           w.write_trace_block(trace_headers, samples)

``write_trace_block`` accepts float32-compatible samples and the C++ backend
encodes them according to the SEG-Y sample format code.

Set ``sample_interval_us`` to patch output ``dt`` while preserving the rest of
the trace header bytes:

.. code-block:: python

   b.sample_interval_us(4000)


Raw Sample Bytes
----------------

Use ``write_raw_trace_block`` when samples are already encoded in SEG-Y byte
format and should be copied without conversion.

.. code-block:: python

   b = cigsegy.SegyWriter.from_headers("raw_out.sgy")
   b.textual(textual).binary(binary)
   b.sample_format(1).sample_count(1001)
   b.overwrite(True)

   with b.open() as w:
       for trace_headers, sample_bytes in raw_blocks:
           w.write_raw_trace_block(trace_headers, sample_bytes)


Create Mode
===========

Use ``create`` when no source SEG-Y header is available and a regular SEG-Y is
enough.


3D Post-Stack
-------------

.. code-block:: python

   b = cigsegy.SegyWriter.create(
       "created.sgy",
       shape=(589, 762, 1001),
       sample_format=5,
       sample_interval_us=2000,
       overwrite=True,
   )
   b.grid(iline_start=1000, xline_start=2000, iline_step=1, xline_step=1)
   b.origin(x_start=600000, y_start=4100000, x_step=25, y_step=25)

   with b.open() as w:
       w.write(data)


Block Writing
-------------

``start`` is the zero-based logical block origin in the target volume.

.. code-block:: python

   b = cigsegy.SegyWriter.create("created.sgy", shape=(589, 762, 1001))
   b.sample_format(5).sample_interval_us(2000).overwrite(True)

   with b.open() as w:
       for ibeg, block in inline_blocks:
           w.write_block(block, start=(ibeg, 0, 0))


2D Trace Collection
-------------------

.. code-block:: python

   b = cigsegy.SegyWriter.create(
       "line.sgy",
       shape=(12000, 1500),
       sample_format=5,
       sample_interval_us=2000,
       overwrite=True,
   )

   with b.open() as w:
       w.write(line_data)


4D Pre-Stack
------------

.. code-block:: python

   b = cigsegy.SegyWriter.create(
       "gather.sgy",
       shape=(120, 200, 48, 1500),
       sample_format=5,
       sample_interval_us=2000,
       overwrite=True,
   )
   b.as_4d()
   b.grid(iline_start=1000, xline_start=2000, offset_start=100)

   with b.open() as w:
       w.write(gather_data)


Custom Geometry
---------------

Use ``geometry`` when regular starts and steps are not enough.  Values can be
scalar, 1D file-order arrays, or arrays matching the trace grid shape.

.. code-block:: python

   b = cigsegy.SegyWriter.create("created.sgy", shape=data.shape)
   b.geometry(iline=iline_grid, xline=xline_grid, x=x_grid, y=y_grid)

   with b.open() as w:
       w.write(data)


Choosing a Mode
===============

Use ``from_template`` when output traces correspond to source traces.

Use ``from_headers`` when another container already stores file-order headers.

Use ``create`` when regular generated headers are good enough.

Use ``write_raw_trace_block`` when byte-preserving sample output matters.


Low-Level Backend
=================

``cigsegy.cpp._CXX_SEGY.SegyBlockWriter`` is the low-level C++ backend used by
``SegyWriter``.  It expects file-order trace headers and sample blocks.  Most
user code should use ``SegyWriter`` instead.
