.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/logo.svg
   :alt: logo

**cigsegy** is a Python and C++ toolkit for reading, inspecting, converting,
and writing SEG-Y seismic data.

The common workflow is intentionally simple:

1. inspect the textual header;
2. scan metadata and let ``cigsegy`` infer geometry byte locations;
3. read samples into a NumPy array with ``fromfile``;
4. process the array;
5. write a SEG-Y back with ``SegyWriter``.

``SegyNP`` is the array-like reader to use when you need partial reads, random
slicing, visualization, or files that should stay on disk.

Source code is available at
`github.com/JintaoLee-Roger/cigsegy <https://github.com/JintaoLee-Roger/cigsegy>`_.


Quick Start
===========

Install
-------

.. code-block:: bash

   pip install cigsegy


Inspect the Header
------------------

Start with the textual header.  It often tells you where inline, crossline, and
coordinate fields are stored.

.. code-block:: python

   import cigsegy

   cigsegy.textual_header("input.sgy")

Then scan metadata.  ``cigsegy`` will try to infer inline, crossline, offset,
step, and coordinate byte locations from the trace headers.

.. code-block:: python

   cigsegy.metaInfo("input.sgy")

Pass explicit byte locations only when the guess is ambiguous or the SEG-Y uses
non-standard headers:

.. code-block:: python

   cigsegy.metaInfo("input.sgy", iline=189, xline=193)

Use ``tools.read_header`` when you need decoded binary or trace header fields:

.. code-block:: python

   binary = cigsegy.tools.read_header("input.sgy", type="bh", printstr=False)
   trace0 = cigsegy.tools.read_header("input.sgy", type="th", n=0, printstr=False)


Read to NumPy
-------------

For many processing scripts, the most direct path is reading the whole volume
into a NumPy array.

.. code-block:: python

   data = cigsegy.fromfile("input.sgy")
   print(data.shape)  # (n_inline, n_xline, n_sample)

For 4D/pre-stack data, try the same automatic path first:

.. code-block:: python

   gathers = cigsegy.fromfile("gather.sgy")

If the inferred geometry is not what you expect, pass the known fields:

.. code-block:: python

   data = cigsegy.fromfile("input.sgy", iline=189, xline=193)
   gathers = cigsegy.fromfile("gather.sgy", iline=189, xline=193, offset=37)

For 2D trace collections:

.. code-block:: python

   traces = cigsegy.collect("line.sgy")

For volumes that are too large for memory, stream samples directly to disk:

.. code-block:: python

   cigsegy.tofile("input.sgy", "samples.dat")  # raw float32, no shape metadata
   cigsegy.to_npy("input.sgy", "samples.npy")  # .npy, can be memory-mapped


Process and Write Back
----------------------

When the output keeps the same trace geometry as the input, write with
``SegyWriter.from_template``.  This is the modern replacement for the old
``create_by_sharing_header`` workflow in new code.

.. code-block:: python

   processed = process(data)

   b = cigsegy.SegyWriter.from_template("input.sgy", "processed.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(processed)

For a continuous sub-volume:

.. code-block:: python

   sub = data[100:300, 40:200, 0:800]

   b = cigsegy.SegyWriter.from_template("input.sgy", "sub.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(sub, start=(100, 40, 0))

For time super-resolution or downsampling where the sample interval changes,
set the output interval explicitly:

.. code-block:: python

   b = cigsegy.SegyWriter.from_template("input_2ms.sgy", "output_1ms.sgy")
   b.strict(False).sample_interval_us(1000).overwrite(True)

   with b.open() as w:
       w.write(super_res_data)

For spatial thinning, copy headers only from real traces in the template:

.. code-block:: python

   thin = data[::2, ::3, :]

   b = cigsegy.SegyWriter.from_template("input.sgy", "thin.sgy")
   b.select(iline=slice(None, None, 2), xline=slice(None, None, 3))
   b.overwrite(True)

   with b.open() as w:
       w.write(thin)


Read Lazily with SegyNP
-----------------------

Use ``SegyNP`` when the file is too large to load, or when you need interactive
random access.

.. code-block:: python

   vol = cigsegy.SegyNP("input.sgy")

   iline = vol[100, :, :]
   xline = vol[:, 200, :]
   time_slice = vol[:, :, 300]
   small_cube = vol[100:140, 200:260, 300:700]


Write from Stored Headers
-------------------------

Use this mode for chunked containers or pipelines that already have textual,
binary, and trace headers.

.. code-block:: python

   b = cigsegy.SegyWriter.from_headers("out.sgy")
   b.textual(textual).binary(binary)
   b.sample_format(1).sample_count(1001)
   b.overwrite(True)

   with b.open() as w:
       for trace_headers, samples in blocks:
           w.write_trace_block(trace_headers, samples)

``write_raw_trace_block`` copies already encoded sample bytes without IBM/IEEE
or integer conversion.


Create Regular Headers from Scratch
-----------------------------------

.. code-block:: python

   b = cigsegy.SegyWriter.create(
       "created.sgy",
       shape=(589, 762, 1001),
       sample_format=5,
       sample_interval_us=2000,
       overwrite=True,
   )
   b.grid(iline_start=1000, xline_start=2000)
   b.origin(x_start=600000, y_start=4100000, x_step=25, y_step=25)

   with b.open() as w:
       w.write(data)


Legacy Shortcuts
================

The following functions are still available:

- ``cigsegy.fromfile``: read a 3D or 4D volume into memory;
- ``cigsegy.collect``: read traces as a 2D array;
- ``cigsegy.tofile``: dump sample data to a raw binary file;
- ``cigsegy.create_by_sharing_header``: create a SEG-Y using headers from an
  existing SEG-Y file;
- ``cigsegy.create``: older 3D regular-volume writer.

For new writing code, prefer ``SegyWriter`` because it covers template writes,
stored-header block writes, raw-byte writes, and generated headers with one API.


License
=======

cigsegy is distributed under the MIT license.


Citation
========

.. code-block:: text

   Li, Jintao. "CIGSEGY: A tool for exchanging data between SEG-Y format and NumPy array inside Python environment".
   URL: https://github.com/JintaoLee-Roger/cigsegy
