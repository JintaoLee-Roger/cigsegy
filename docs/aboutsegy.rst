About SEG-Y
###########

The SEG-Y (sometimes SEG Y or SEGY) file format is one of several data 
standards developed by the Society of Exploration Geophysicists (SEG) 
for the exchange of geophysical data. It is an open standard, and is 
controlled by the SEG Technical Standards Committee, a non-profit organization.

The format was originally developed in 1973 to store single-line seismic 
reflection digital data on magnetic tapes. The specification was 
published in 1975.

The format and its name evolved from the SEG "Ex" or Exchange 
Tape Format. However, since its release, there have been significant 
advancements in geophysical data acquisition, such as 3-dimensional 
seismic techniques and high speed, high capacity recording.

The most recent revision of the SEG-Y format was published in 2002, 
named the rev 1 specification. It still features certain legacies 
of the original format (referred as rev 0), such as an optional SEG-Y 
tape label, the main 3200 byte textual file header and a 400 byte 
binary file header [1]_.


Version
=======

There are several version of SEG-Y, you can find the documents at `SEG Technical Standards <https://library.seg.org/seg-technical-standards>`_.

The marjor versions include:

- `SEG-Y v0 (1975) <https://library.seg.org/pb-assets/technical-standards/seg_y_rev0-1686080980707.pdf>`_
- `SEG-Y v1 (2002) <https://library.seg.org/pb-assets/technical-standards/seg_y_rev1-1686080991247.pdf>`_
- `SEG-Y v2 (2017) <https://library.seg.org/pb-assets/technical-standards/seg_y_rev2_0-mar2017-1686080998003.pdf>`_

All version compression (click `here <https://wiki.seg.org/w/images/4/42/SEG-Y_bytestream_all_revisions.pdf>`_ for source pdf version):

.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/comparison.png
    :alt: segy comparison
    :align: center


Structure
=========

By comparing several versions of the SEG-Y files, we can observe that 
the SEG-Y file is primarily composed of the following **mandatory parts**: 
``3200 bytes Textual Header``, ``400 bytes Binary Header``, ``240 bytes Trace Headers`` 
and corresponding ``data traces``. And there are several **optional parts** in 
version 1 and version 2: ``Optional 128 byte SEG-Y Tape Label``, ``Extended
Textual File Header``, ``Data Trailer``.

**CIGSEGY** only supports these **mandatory parts**, 
and other **optional parts** are currently not supported 
(this is because the vast majority of SEG-Y files only include these **mandatory parts**), 
i.e., **cigsegy** considers that a SEG-Y file contains **one** ``Textual Header``, 
**one** ``Binary Header``, **N** ``Trace Headers``, **N** ``data traces``, 
where N is the total number of traces.

Although the SEG-Y file structure and Version 0 closely resemble 
what cigsegy considers, this doesn't imply that cigsegy can only 
handle Version 0 SEG-Y files. 
CIGSEGY takes into account the information from 
the binary header and trace headers of Version 1 and Version 2 as well.


Type of data
============

The seismic data can be devided into two categories: Poststack seismic data
and Prestack seismic gather. 

Inside the binary trace header of a SEG-Y file, there is a 2-byte integer 
stored at bytes 3229-3230. This integer describes the category to 
which this SEG-Y file belongs.

+-------+--------------------------------+-------------------+
| value | data type                      | post or pre stack |
+=======+================================+===================+
| -1    | Other                          | \-                |
+-------+--------------------------------+-------------------+
| 0     | Unknown                        | \-                |
+-------+--------------------------------+-------------------+
| 1     | As recorded (no sorting)       | prestack          |
+-------+--------------------------------+-------------------+
| 2     | CDP ensemble                   | prestack          |
+-------+--------------------------------+-------------------+
| 3     | Single fold continuous profile | \-                |
+-------+--------------------------------+-------------------+
| 4     | Horizontally stacked           | poststack         |
+-------+--------------------------------+-------------------+
| 5     | Common source point            | prestack          |
+-------+--------------------------------+-------------------+
| 6     | Common receiver point          | prestack          |
+-------+--------------------------------+-------------------+
| 7     | Common offset point            | prestack          |
+-------+--------------------------------+-------------------+
| 8     | Common mid-point               | prestack          |
+-------+--------------------------------+-------------------+
| 9     | Common conversion point        | prestack          |
+-------+--------------------------------+-------------------+

--------------------

.. [1] `From SEG wiki: SEG-Y (https://wiki.seg.org/wiki/SEG-Y) <https://wiki.seg.org/wiki/SEG-Y>`_

