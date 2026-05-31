Read, Process, and Write with NumPy
###################################

This page is for the common script workflow: inspect a SEG-Y file, read samples
into NumPy, process the array, and write or export the result.

Use this path when:

- the data fits in memory and you want a plain ``numpy.ndarray``;
- you are building a conversion or processing script;
- you want explicit control over file creation;
- you need quick header inspection without keeping a SEG-Y object around.

Use :doc:`SegyNP` instead when the file should stay on disk and you need lazy
slicing, repeated random access, coordinate transforms, arbitrary lines, or
in-place editing through an array-like object.

Most functions try to infer inline, crossline, offset, step, and coordinate byte
locations from trace headers.  Pass explicit locations only when the inferred
geometry is wrong or ambiguous.


Inspect Headers and Metadata
============================

Start with the textual header:

.. code-block:: python

   import cigsegy

   cigsegy.textual_header("input.sgy")

Get a geometry summary:

.. code-block:: python

   cigsegy.metaInfo("input.sgy")

If the automatic scan needs help:

.. code-block:: python

   cigsegy.metaInfo("input.sgy", iline=189, xline=193)
   cigsegy.metaInfo("gather.sgy", iline=189, xline=193, offset=37, is4d=True)

Use ``tools.get_metaInfo`` when code needs the metadata dictionary:

.. code-block:: python

   meta = cigsegy.tools.get_metaInfo("input.sgy")
   print(meta["ni"], meta["nx"], meta["nt"], meta["dt"])

Read decoded binary or trace headers:

.. code-block:: python

   binary = cigsegy.tools.read_header("input.sgy", type="bh", printstr=False)
   trace0 = cigsegy.tools.read_header("input.sgy", type="th", n=0, printstr=False)

   print(binary[17])   # sample interval
   print(binary[21])   # samples per trace
   print(trace0[189])  # inline if byte 189 is used by this file


Read Samples
============

``fromfile`` scans geometry and returns a NumPy array:

.. code-block:: python

   data = cigsegy.fromfile("poststack.sgy")
   print(data.shape)  # (n_inline, n_xline, n_sample)

For prestack data, try the same automatic path first:

.. code-block:: python

   gathers = cigsegy.fromfile("prestack.sgy")
   print(gathers.shape)  # (n_inline, n_xline, n_offset, n_sample)

If the geometry needs explicit fields:

.. code-block:: python

   data = cigsegy.fromfile("poststack.sgy", iline=189, xline=193)
   gathers = cigsegy.fromfile(
       "prestack.sgy",
       iline=189,
       xline=193,
       offset=37,
       is4d=True,
   )

``collect`` reads file-order traces as ``(n_trace, n_sample)``.  It is useful
for 2D lines, irregular files, quick sampling, and index-based reads:

.. code-block:: python

   traces = cigsegy.collect("line.sgy")
   part = cigsegy.collect("line.sgy", beg=1000, end=2000)
   trace = cigsegy.collect("line.sgy", beg=100)
   window = cigsegy.collect("line.sgy", beg=1000, end=2000, tbeg=200, tend=800)

Read arbitrary trace indices:

.. code-block:: python

   import numpy as np

   indices = np.array([0, 50, 100], dtype=np.int32)
   traces = cigsegy.collect("line.sgy", indices=indices)


Read Trace Header Keys
======================

``get_trace_keys`` reads selected trace header fields over many traces:

.. code-block:: python

   keys = cigsegy.get_trace_keys("input.sgy", keyloc=[189, 193], beg=0, end=1000)
   ilines = keys[:, 0]
   xlines = keys[:, 1]

Read all trace values for one key:

.. code-block:: python

   ilines = cigsegy.get_trace_keys("input.sgy", keyloc=189)

For non-standard byte locations, pass ``force`` as the byte width:

.. code-block:: python

   values = cigsegy.get_trace_keys("input.sgy", keyloc=221, beg=0, end=1000, force=4)


Export Raw Samples
==================

``tofile`` writes sample values only, without textual, binary, or trace headers.
The output is raw little-endian IEEE float32 samples:

.. code-block:: python

   cigsegy.tofile("poststack.sgy", "poststack.dat")

The main reason to use ``tofile`` is memory pressure.  ``fromfile`` returns a
NumPy array and therefore needs enough RAM for the whole volume.  ``tofile``
streams through the SEG-Y and writes the samples directly to disk, which is
useful when the machine cannot hold the full volume in memory or when another
program expects raw float32 data.

The raw output has no shape metadata.  Record the shape from ``metaInfo`` or
``tools.get_metaInfo`` and reopen it with ``numpy.memmap`` when needed:

.. code-block:: python

   meta = cigsegy.tools.get_metaInfo("poststack.sgy")
   shape = (meta["ni"], meta["nx"], meta["nt"])

   data = np.memmap("poststack.dat", dtype="<f4", mode="r", shape=shape)

Use ``as2d=True`` when you want a simple trace-by-sample dump without geometry:

.. code-block:: python

   cigsegy.tofile("input.sgy", "traces.dat", as2d=True)

If you want a NumPy file with shape and dtype metadata, use ``to_npy``.  It also
streams data and does not materialize the full volume in memory:

.. code-block:: python

   shape = cigsegy.to_npy("poststack.sgy", "poststack.npy")
   data = np.load("poststack.npy", mmap_mode="r")


Write SEG-Y
===========

For new writing code, use :doc:`writer`.  The most common path is
``SegyWriter.from_template``:

.. code-block:: python

   processed = process(data)

   b = cigsegy.SegyWriter.from_template("poststack.sgy", "processed.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(processed)

Write a continuous sub-volume by passing its logical start:

.. code-block:: python

   sub = data[100:300, 40:200, 0:800]

   b = cigsegy.SegyWriter.from_template("poststack.sgy", "sub.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(sub, start=(100, 40, 0))

For time resampling, set the output sample interval explicitly:

.. code-block:: python

   b = cigsegy.SegyWriter.from_template(
       "input_2ms.sgy",
       "output_1ms.sgy",
       sample_interval_us=1000,
       overwrite=True,
   )
   b.strict(False)

   with b.open() as w:
       w.write(super_res_data)

Older shortcuts such as ``create_by_sharing_header`` and ``create`` remain
available for compact scripts and legacy code.  See :doc:`create`.


Edit Header Values
==================

Header editing changes the SEG-Y file in place.  Work on a copy unless you are
intentionally modifying the original file.

Modify a binary header key:

.. code-block:: python

   cigsegy.modify_bin_key("copy.sgy", loc=17, value=2000)

Modify one trace header key:

.. code-block:: python

   cigsegy.modify_trace_key("copy.sgy", loc=189, value=1024, idx=0)

Use ``idx=-1`` to modify all traces:

.. code-block:: python

   cigsegy.modify_trace_key("copy.sgy", loc=117, value=2000, idx=-1)


Irregular or Unsorted Geometry
==============================

For convenient random access, prefer ``SegyNP(..., as_unsorted=True)``.  The
lower-level tools are still available when you want direct control:

.. code-block:: python

   geom = cigsegy.tools.full_scan("input.sgy", 189, 193)
   data = cigsegy.tools.load_by_geom("input.sgy", geom)

Full scans read all trace headers, so they are slower than regular geometry
scans.


Plot Helpers
============

Plotting helpers are useful while confirming header locations and geometry:

.. code-block:: python

   cigsegy.plot.plot_trace_keys("input.sgy", keyloc=193, beg=0, end=2000)
   cigsegy.plot.plot_trace_ix("input.sgy")
   cigsegy.plot.plot_region("input.sgy")

Pass explicit locations when checking a suspected geometry:

.. code-block:: python

   cigsegy.plot.plot_trace_ix("input.sgy", iline=189, xline=193)
   cigsegy.plot.plot_region("input.sgy", iline=189, xline=193)


Sample Format Helpers
=====================

For low-level conversion work, ``cigsegy`` exposes IBM/IEEE floating-point
helpers:

.. code-block:: python

   ieee = cigsegy.ibm_to_ieee(ibm_values, is_big_endian=True)
   ibm = cigsegy.ieee_to_ibm(
       ieee_values,
       is_little_endian_input=True,
       is_big_endian_output=True,
   )

Most users do not need these functions directly; SEG-Y readers and writers use
the file sample format code automatically.
