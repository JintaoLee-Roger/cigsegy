Quickstart
##########

This page is the shortest path from installation to a complete SEG-Y workflow:
inspect the file, read samples into NumPy, process them, and write a new SEG-Y.


Install
=======

Install from PyPI:

.. code-block:: bash

   pip install cigsegy

For local development, use an editable install from the repository root:

.. code-block:: bash

   pip install -U pip
   pip install -e . --config-settings editable_mode=strict

The package includes the compiled ``cigsegy.cpp._CXX_SEGY`` extension and the
Python command-line entry point ``cigsegy``.  A C++ compiler and ``pybind11`` are
needed when building from source.


Choose a Workflow
=================

Most users work with SEG-Y in one of two ways:

- :doc:`read the data into NumPy <core/post>` when you want a simple script:
  inspect the file, load samples, process an array, export raw samples, or write
  a new SEG-Y;
- :doc:`work with the SEG-Y like an array <core/SegyNP>` when the file should
  stay on disk and you need slicing, random access, visualization, coordinate
  transforms, arbitrary lines, or careful in-place edits.

When unsure, start by reading into NumPy.  Switch to ``SegyNP`` when repeated
partial reads or geometry-aware operations become more important than loading
the whole volume at once.


Inspect the SEG-Y
=================

Start with the textual header.  It often records byte locations for inline,
crossline, coordinate, and offset fields.

.. code-block:: python

   import cigsegy

   cigsegy.textual_header("input.sgy")

Then scan the geometry:

.. code-block:: python

   cigsegy.metaInfo("input.sgy")

By default, ``cigsegy`` tries to infer inline, crossline, offset, step, and
coordinate byte locations from trace headers.  If the guess is ambiguous or the
SEG-Y uses non-standard headers, pass the known byte locations explicitly:

.. code-block:: python

   cigsegy.metaInfo("input.sgy", iline=189, xline=193)
   cigsegy.metaInfo("gather.sgy", iline=189, xline=193, offset=37, is4d=True)


Read to NumPy
=============

For a regular post-stack volume, read the samples into a plain NumPy array:

.. code-block:: python

   data = cigsegy.fromfile("input.sgy")

The returned array is usually shaped as ``(n_inline, n_crossline, n_sample)``.
For prestack data:

.. code-block:: python

   gathers = cigsegy.fromfile("gather.sgy")

If the inferred geometry is not what you expect, provide the fields you know:

.. code-block:: python

   data = cigsegy.fromfile("input.sgy", iline=189, xline=193)
   gathers = cigsegy.fromfile("gather.sgy", iline=189, xline=193, offset=37, is4d=True)

When geometry is irregular or file-order traces are all you need, use
``collect``:

.. code-block:: python

   traces = cigsegy.collect("line.sgy")
   window = cigsegy.collect("line.sgy", beg=0, end=1000, tbeg=100, tend=800)

When the volume is too large for memory, stream it to disk instead of calling
``fromfile``:

.. code-block:: python

   cigsegy.tofile("input.sgy", "samples.dat")       # raw float32, no shape metadata
   cigsegy.to_npy("input.sgy", "samples.npy")       # .npy with shape metadata


Process and Write Back
======================

After processing with NumPy, write a new SEG-Y with ``SegyWriter`` while sharing
headers from a template SEG-Y:

.. code-block:: python

   processed = data * 2.0

   builder = cigsegy.SegyWriter.from_template("input.sgy", "processed.sgy")
   builder.overwrite(True)

   with builder.open() as writer:
       writer.write(processed)

To change the sample interval while keeping the same spatial geometry:

.. code-block:: python

   builder = cigsegy.SegyWriter.from_template(
       "input_2ms.sgy",
       "output_1ms.sgy",
       sample_interval_us=1000,
       overwrite=True,
   )
   builder.strict(False)

   with builder.open() as writer:
       writer.write(resampled)


Use the CLI
===========

The installed ``cigsegy`` command mirrors the same common workflow:

.. code-block:: bash

   cigsegy textual input.sgy
   cigsegy meta input.sgy
   cigsegy fromfile input.sgy volume.npy
   cigsegy to-npy input.sgy volume.npy
   cigsegy create input.sgy processed.npy processed.sgy --overwrite

The CLI is useful for inspection and simple conversions.  For production
processing pipelines, Python code is usually clearer.  See :doc:`cli` for all
subcommands and the optional native C++ tool.


Use Lazy Reading When Needed
============================

Use :class:`cigsegy.SegyNP` when the file is too large to load fully, or when
you need repeated slicing and visualization-oriented reads:

.. code-block:: python

   vol = cigsegy.SegyNP("input.sgy")
   inline = vol[100]


.. toctree::
   :caption: Core Workflows
   :maxdepth: 2
   :hidden:

   core/post
   core/writer
   core/SegyNP
   core/pre
   core/create
   cli

.. toctree::
   :caption: Reference
   :maxdepth: 2
   :hidden:

   aboutsegy
   core/headers
   api/pyapi
   changelog
   ends
