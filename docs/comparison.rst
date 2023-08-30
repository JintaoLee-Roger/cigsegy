Comparison
##########

Here we comparison ``cigsegy`` 
to `segyio <https://github.com/equinor/segyio>`_ 
and `cegysak <https://github.com/trhallam/segysak>`_.

Compared with ``segysak``, our implementation is much more faster.
``cigsegy`` is a little slower than ``segyio`` when reading a segy file, but the gap is small. However,
``cigsegy`` is faster than ``segyio`` when creating a segy file.

``segyio`` assumes that the file is a regular, sorted 3D volume.
It also supports files (unstrict mode) that are just a collection of traces too,
but lots of features are now disabled and will raise errors in this mode.

However, there are losts of segy files that are sorted, but missing some traces.
``segyio`` doesn't support these files although they are easy dealed with.
``cigsegy`` supports these files just using the same methods, 
e.g. ``cigsegy.fromfile('miss.segy')``. 
``cigsegy`` can also handle the cases that the inline/crossline step is not 1.

For some reasons (confidentiality requirement, mistakes when recording, ...), 
the file headers are broken. If you remember the volume size and sample format (1 for IBM, 5 for IEEE), 
``cigsegy`` can read the files too. Just using 

.. code-block::

    d = cigsegy.fromfile_ignore_header('miss.segy', inline_size, crossline_size, time, dformat)

``segymat`` is implememted in MATLAB, it runs very slowly.

``shape = (651, 951, 462) 1.3G``

+--------+---------+--------+---------+---------+
|  mode  | cigsegy | segyio | segysak | segymat |
+========+=========+========+=========+=========+
|  read  | 1.212s  | 0.944s | 134.8s  | 151.99s |
+--------+---------+--------+---------+---------+
| create |  3.01s  | 14.03s |   \--   |   \--   |
+--------+---------+--------+---------+---------+



Run ``read`` 3 times for ``cigsegy`` and ``segyio``, 
once for ``segysak``, 
``(1062, 2005, 2401) 20G``


+--------+-----------------+---------------+---------+---------+
|  mode  |     cigsegy     |    segyio     | segysak | segymat |
+========+=================+===============+=========+=========+
|  read  | 45.6s+45.5s+45s |  65s+20s+19s  | 612.45s | >1500s  |
+--------+-----------------+---------------+---------+---------+
| create |  41.14s+43.35s  | 78.74s+82.74s |   \--   |   \--   |
+--------+-----------------+---------------+---------+---------+

The second reading by ``segyio`` is faster than the first reading by ``segyio``, 
while the time required for reading three times by ``cigsegy`` is very close.