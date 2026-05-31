Work with SEG-Y as an Array
###########################

``SegyNP`` presents a SEG-Y file as a NumPy-like array while keeping samples on
disk.  Use it when you want to slice a SEG-Y repeatedly, inspect geometry,
convert between survey coordinates and inline/crossline positions, extract
arbitrary lines, or edit an existing file carefully in writable mode.

If your workflow is simply "read the whole SEG-Y, process it as a NumPy array,
and write it back", start with :doc:`post` instead.


Open a SEG-Y
============

For a regular 3D post-stack volume:

.. code-block:: python

   import cigsegy

   vol = cigsegy.SegyNP("poststack.sgy")
   print(vol.shape)  # (n_inline, n_xline, n_sample)
   print(vol.ndim)

``SegyNP`` tries to infer geometry automatically.  If the inferred geometry is
wrong or ambiguous, pass explicit key locations:

.. code-block:: python

   vol = cigsegy.SegyNP("poststack.sgy", keylocs=[189, 193])

For 4D/prestack data, try automatic inference first:

.. code-block:: python

   gathers = cigsegy.SegyNP("prestack.sgy")
   print(gathers.shape)  # (n_inline, n_xline, n_offset, n_sample)

If needed, provide inline, crossline, and offset byte locations:

.. code-block:: python

   gathers = cigsegy.SegyNP("prestack.sgy", keylocs=[189, 193, 37])

Use ``ndim=2`` when the file should be treated as file-order traces by samples:

.. code-block:: python

   traces = cigsegy.SegyNP("line.sgy", ndim=2)
   print(traces.shape)  # (n_trace, n_sample)


Choose a View Mode
==================

``SegyNP`` can open a file in several ways:

- default scan mode: scan regular 3D/4D geometry and expose logical axes;
- ``ndim=2``: expose the file as ``(n_trace, n_sample)``;
- ``as_unsorted=True``: build a map for unsorted geometry;
- ``shape_hint=(...)``: skip geometry scanning when the regular shape is known.

Examples:

.. code-block:: python

   regular = cigsegy.SegyNP("regular.sgy")
   traces = cigsegy.SegyNP("line.sgy", ndim=2)
   lazy = cigsegy.SegyNP("regular.sgy", shape_hint=(589, 762, 1001))

   unsorted = cigsegy.SegyNP(
       "unsorted.sgy",
       keylocs={"iline": 189, "xline": 193},
       as_unsorted=True,
   )

``as_unsorted=True`` scans trace headers and builds a geometry map, so opening
can be slower than regular scan mode.


Array-Like Reading
==================

Use NumPy-style indexing:

.. code-block:: python

   iline = vol[100, :, :]
   xline = vol[:, 200, :]
   time_slice = vol[:, :, 300]
   cube = vol[100:160, 200:260, 300:900]

For 2D trace views:

.. code-block:: python

   trace = traces[100]
   part = traces[1000:1200, :]

For prestack data:

.. code-block:: python

   one_gather = gathers[20, 30, :, :]
   near_offsets = gathers[:, :, :8, :]

Convert to a regular NumPy array when needed:

.. code-block:: python

   data = vol.to_numpy()

or rely on NumPy conversion:

.. code-block:: python

   data = np.asarray(vol)

For whole-volume reads that should always be materialized immediately,
:func:`cigsegy.fromfile` is usually simpler.


Metadata and Header Access
==========================

``SegyNP`` keeps common metadata on the object:

.. code-block:: python

   print(vol.metainfo)
   print(vol.keylocs)
   print(vol.ntrace, vol.nt, vol.dtype)

Read the textual header:

.. code-block:: python

   vol.textual_header()
   text = vol.textual_header(printtext=False)

Read binary and trace header values:

.. code-block:: python

   dt = vol.bkeyi2(17)
   ns = vol.bkeyi2(21)
   iline0 = vol.keyi4(0, 189)

Convenience header arrays are available for common trace keys:

.. code-block:: python

   ilines = vol.iline[:]
   xlines = vol.xline[:]
   xs = vol.coordx[:]
   ys = vol.coordy[:]

The available convenience names include ``iline``, ``xline``, ``offset``,
``coordx``, ``coordy``, and ``itrace``.


Coordinate Transforms
=====================

For regular 3D or 4D data, ``SegyNP`` can convert between logical
inline/crossline coordinates and survey ``x/y`` coordinates:

.. code-block:: python

   xy = vol.ix_to_xy([[100, 200], [120, 220]])
   ix = vol.xy_to_ix(xy)

By default, ``ix_to_xy`` and ``xy_to_ix`` treat inline/crossline values as
zero-based logical indices.  Pass ``zero_origin=False`` when using real inline
and crossline numbers from the trace headers:

.. code-block:: python

   xy = vol.ix_to_xy([[1024, 2300]], zero_origin=False)
   ix = vol.xy_to_ix(xy, zero_origin=False)

The two directions are fitted separately from survey control points:

- ``ix_to_xy`` fits an inline/crossline-to-``x/y`` affine transform;
- ``xy_to_ix`` fits an ``x/y``-to-inline/crossline affine transform.

This is intentionally not implemented as one matrix plus its inverse.  Real
SEG-Y coordinates are often rounded, scaled, or slightly noisy, and fitting the
two directions independently is usually more stable for practical use.


Geometry Maps
=============

For unsorted files or custom trace mapping, create or update a geometry map:

.. code-block:: python

   vol = cigsegy.SegyNP("unsorted.sgy", keylocs=[189, 193], as_unsorted=True)
   trace_indices = vol.map_to_indices([[10, 20], [11, 20]])

You can also provide a geometry array yourself:

.. code-block:: python

   vol.update_geometry(geom)


Plotting and Arbitrary Lines
============================

``SegyNP`` exposes plotting helpers for geometry checks and interactive work:

.. code-block:: python

   vol.plot_region()
   vol.plot_trace_keys(keyloc=193, beg=0, end=2000)
   vol.plot3d()

Extract an arbitrary line from zero-origin indices, inline/crossline numbers,
or survey ``x/y`` coordinates:

.. code-block:: python

   points = [[100, 200], [120, 260], [180, 300]]
   line, path, indices = vol.arbitrary_line(points, ptype="zero")

   xy_points = [[443210.0, 3378120.0], [444020.0, 3378800.0]]
   line, path, indices = vol.arbitrary_line(xy_points, ptype="xy")

For interactive picking:

.. code-block:: python

   out = vol.extract_arbitrary_line_by_view()


Write Through SegyNP
====================

Open with ``mode="rw"`` only when you intentionally want to modify an existing
SEG-Y in place.  Make a copy first if the source file matters.

.. code-block:: python

   vol = cigsegy.SegyNP("editable.sgy", mode="rw")
   vol[100:110, 200:220, :] = processed
   vol.close()

You can also edit header values through the same object:

.. code-block:: python

   vol = cigsegy.SegyNP("editable.sgy", mode="rw")
   vol.set_bkeyi2(17, 2000)
   vol.set_keyi4(0, 189, 1024)
   vol.close()

For creating a new SEG-Y, use :doc:`writer` instead of modifying a source file.


Export from SegyNP
==================

Dump sample values to a raw binary file:

.. code-block:: python

   vol.tofile("samples.dat")

This is useful when the data is too large for memory and a downstream tool can
read raw float32 samples.  Raw files do not store shape metadata.  If you need a
``.npy`` file that can be memory-mapped later, use ``cigsegy.to_npy``:

.. code-block:: python

   cigsegy.to_npy("input.sgy", "samples.npy")
   data = np.load("samples.npy", mmap_mode="r")

For a new SEG-Y file, prefer ``SegyWriter``:

.. code-block:: python

   b = cigsegy.SegyWriter.from_template("input.sgy", "out.sgy")
   b.overwrite(True)

   with b.open() as w:
       w.write(vol[...])


Choosing a Reader
=================

Use ``SegyNP`` when:

- you need slices, traces, gathers, or time slices;
- you do not want to load the whole file into memory;
- you need coordinate transforms or arbitrary lines;
- the file may be used by an interactive viewer;
- you intentionally want in-place editing through ``mode="rw"``.

Read into NumPy first when:

- you want a whole 3D or 4D NumPy array immediately;
- the data comfortably fits in memory;
- you are writing a simple conversion or batch-processing script;
- you need one-off header inspection.
