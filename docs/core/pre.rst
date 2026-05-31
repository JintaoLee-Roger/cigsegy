Pre-Stack and 4D Data
#####################

Pre-stack SEG-Y data usually needs three geometry fields:

- inline;
- crossline or CDP;
- offset.

In ``cigsegy`` this maps naturally to a 4D array:
``(n_inline, n_xline, n_offset, n_sample)``.


Scan Metadata
=============

.. code-block:: python

   import cigsegy

   cigsegy.metaInfo("gather.sgy")

``cigsegy`` will try to infer inline, crossline, offset, and geometry steps.
If the guess is wrong or ambiguous, pass explicit fields:

.. code-block:: python

   cigsegy.metaInfo("gather.sgy", iline=189, xline=193, offset=37, is4d=True)

If the file is a 2D line gather, use the line/CDP field as ``iline`` and the
offset field as ``xline``:

.. code-block:: python

   cigsegy.metaInfo("line_gather.sgy", iline=9, xline=37, is4d=False)


Read with SegyNP
================

.. code-block:: python

   gathers = cigsegy.SegyNP("gather.sgy", keylocs=[189, 193, 37])

   one_cmp = gathers[20, 30, :, :]
   near = gathers[:, :, :8, :]
   window = gathers[10:20, 30:50, :, 200:800]


Read into Memory
================

.. code-block:: python

   data = cigsegy.fromfile("gather.sgy")


Unsorted Pre-Stack Files
========================

If trace order is not regular, request an unsorted geometry map:

.. code-block:: python

   gathers = cigsegy.SegyNP(
       "unsorted_gather.sgy",
       keylocs={"iline": 189, "xline": 193, "offset": 37},
       as_unsorted=True,
   )

This is slower to open because all trace headers must be scanned, but later
array indexing uses the geometry map.


Write 4D Data
=============

Use ``SegyWriter.from_template`` when output traces correspond to existing
template traces:

.. code-block:: python

   b = cigsegy.SegyWriter.from_template("gather.sgy", "processed.sgy")
   b.as_4d().overwrite(True)

   with b.open() as w:
       w.write(processed)

Use ``SegyWriter.create(...).as_4d()`` when generating a regular 4D SEG-Y from
scratch:

.. code-block:: python

   b = cigsegy.SegyWriter.create(
       "created_gather.sgy",
       shape=(120, 200, 48, 1500),
       sample_interval_us=2000,
       overwrite=True,
   )
   b.as_4d().grid(iline_start=1000, xline_start=2000, offset_start=100)

   with b.open() as w:
       w.write(gathers)
