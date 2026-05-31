Legacy Creation APIs
####################

For new writing code, use :doc:`writer`.  This page documents the shorter
factory functions that remain useful for old scripts and compact workflows.


create_by_sharing_header
========================

``create_by_sharing_header`` creates a SEG-Y file by copying headers from an
existing SEG-Y file and writing new sample values.

.. code-block:: python

   import cigsegy

   processed = process(data)

   cigsegy.create_by_sharing_header(
       "processed.sgy",
       "template.sgy",
       processed,
       keylocs=[189, 193, 1, 1],
   )

For a raw float32 binary file, pass ``shape``:

.. code-block:: python

   cigsegy.create_by_sharing_header(
       "processed.sgy",
       "template.sgy",
       "processed.dat",
       shape=(589, 762, 1001),
       keylocs=[189, 193, 1, 1],
   )

For a sub-volume, pass the zero-based logical start:

.. code-block:: python

   cigsegy.create_by_sharing_header(
       "sub.sgy",
       "template.sgy",
       sub,
       keylocs=[189, 193, 1, 1],
       start=[100, 40, 0],
   )

For time super-resolution or downsampling where the output sample count or
sample interval no longer matches the template, set ``strict=False`` and pass
``dt_new`` in microseconds:

.. code-block:: python

   cigsegy.create_by_sharing_header(
       "output_1ms.sgy",
       "template_2ms.sgy",
       super_res_data,
       keylocs=[189, 193, 1, 1],
       strict=False,
       dt_new=1000,
   )

This API handles continuous sub-volumes.  For spatial thinning with strides,
use ``SegyWriter.from_template(...).select(...)``.


Textual Header Override
=======================

The ``textual`` argument can be:

- ``""`` or ``None``: keep/generate the default behavior;
- 3200 bytes;
- a 3200-character string;
- a list of strings used to generate a standard 40-line textual header.

.. code-block:: python

   text_lines = [
       "Processed with custom workflow",
       "Input: template.sgy",
   ]

   cigsegy.create_by_sharing_header(
       "processed.sgy",
       "template.sgy",
       processed,
       keylocs=[189, 193, 1, 1],
       textual=text_lines,
   )


create
======

``cigsegy.create`` is an older convenience function for 3D regular volumes.  It
is still available, but ``SegyWriter.create`` is preferred for new code.

.. code-block:: python

   cigsegy.create(
       "created.sgy",
       data,
       format=5,
       dt=2000,
       start_time=0,
       iline_interval=25,
       xline_interval=25,
   )


SegyCreate
==========

``cigsegy.SegyCreate`` is a lower-level builder used by the legacy create path.
Use it only when you need direct access to those older construction steps.
