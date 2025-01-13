.. figure:: https://github.com/JintaoLee-Roger/images/raw/main/cigsegy/assets/logo.svg
    :alt: logo


**A SEG-Y tool developed by** `Computational Interpretation Group (CIG) <https://cig.ustc.edu.cn/main.htm>`_

**cigsegy** is a tool for exchanging data between **SEG-Y** format and 
**NumPy** array inside Python environment.

It can be used to read and convert a SEG-Y format data into a Numpy array, 
even if the SEG-Y format data is missing some traces or its inline/crossline 
step is not equal to 1.

It can also be used to create a SEG-Y format data from a Numpy array. In this
mode, users can use headers from a existed SEG-Y file or create new headers by 
setting some parameters.

Tutorial and reference documentation is provided 
at `cigsegy.readthedocs.io <https://cigsegy.readthedocs.io/en/latest/>`_.
And the source code is always available at `github.com/JintaoLee-Roger/cigsegy <https://github.com/JintaoLee-Roger/cigsegy>`_.



Core Features
=============

- Fast (Implemented in c++)
- python wraping and **numpy** array supports
- dealing with normal and **irregular** SEG-Y volume [1]_.
- creating a SEG-Y file using the **existed header** of a SEG-Y


Quick Start
===========

1. Install cigsegy via PyPi

.. code-block:: bash

    pip install cigsegy


2. Print the 3200 bytes textual header of a SEG-Y file

.. code-block:: python

    >>> import cigsegy
    >>> cigsegy.textual_header('rogan.sgy')
    # C01 CLIENT: STUART PETROLEUM LTD    AREA:COOPER BASIN   SOUTH AUSTRALIA
    # ...
    # C06 INLINE RANGE: 360 - 1684(2)  CROSSLINE RANGE 1764 - 2532(1)
    # C07 -------------PROCESSING FLOW---------------
    # ...
    # C35  DESC                   BYTE LOCATION       FORMAT 
    # C36  3D INLINE NUMBER        9- 12             32 BIT INTEGER 
    # C37  3D CROSSLINE  NUMBER   21- 24             32 BIT INTEGER 
    # C38  CDP_X                  73- 76             32 BIT INTEGER 
    # C39  CDP_Y                  77- 80             32 BIT INTEGER 
    # C40  

You can get some key information to read the SEG-Y file, such as inline location 
is 9 (C36), crossline location is 21 (C37), X location is 73 (C38), Y location 
is 77 (C39), inline step is 2 (C06), crossline step is 1 (C06).

3. Scan the SEG-Y file and get some meta information

.. code-block:: python

    >>> cigsegy.metaInfo('rogan.sgy', iline=9, xline=21, istep=2, xstep=1, xloc=73, yloc=77)
    # In python, the shape is (n-inline, n-crossline, n-time) = (663, 769, 1001).

    # shape: (n-time, n-crossline, n-inline) = (1001, 769, 663)
    # sample interval: 4000, data format code: 4-bytes IBM floating-point
    # inline range: 360 - 1684, crossline range: 1764 - 2532
    # interval of inline: 35.0, interval of crossline: 17.5, time start: 0
    # inline field: 9, crossline field: 21
    # inline step: 2, crossline step: 1
    # Is regular file (no missing traces): false

You will get some information about this SEG-Y file, such as, the data shape, 
intervals, data format ...

.. Note::

    If you are unsure about the values of some parameters, 
    you can ignore them and cigsegy will try to guess them automatically.
 
    .. code-block:: python

        >>> cigsegy.metaInfo('fx.segy', iline=9, xline=21) # ignore istep, xstep, ...


4. Read the SEG-Y

Please note that the shape is like (n-inlines, n-crosslines, n-time_samples)

.. code-block:: python

    >>> d = cigsegy.fromfile('rogan.sgy', iline=9, xline=21, istep=2, xstep=1)
    >>> d.shape
    # (663, 769, 1001)


If you need a binary file without any headers, i.e., save the numpy array

.. code-block:: python

    >>> cigsegy.tofile('rogan.sgy', 'out.dat', iline=9, xline=21, istep=2, xstep=1)

.. Note::
    When using ``cigsegy.tofile()``, you **don't** have to worry about 
    running out of memory. Therefore, this function is very useful when 
    dealing with **huge** files.


5. Create a SEG-Y using a numpy array and headers from another SEG-Y file

There is often such a workflow:
    a. Display SEG-Y format data ``orig.segy`` in specialized software, such as Petrel.
    b. Use Python code to process this data and obtain new data ``afterprocess``, which is in NumPy array format
    c. To display this processed data in specialized software, it needs to be converted back to SEG-Y format and use the headers from the original data, i.e., using the NumPy array ``afterprocess`` and the header of ``orig.segy`` to create a new SEG-Y file ``out.segy``.

.. code-block:: python

    # assume the iline/xline/istep/xstep of **orig.segy** are 9/21/1/1
    >>> cigsegy.create_by_sharing_header('out.segy', 'orig.segy', afterprocess, \
        keylocs=[9, 21])

6. Create a SEG-Y using a numpy array and some parameters

.. code-block:: python

    # d is a numpy array, d.shape == (n-inlines, n-crosslines, n-time)
    >>> cigsegy.create('out.segy', d, format=5, start_time=0, iline_interval=15, ...)


7. Access the SEG-Y file as a 3D numpy array, without reading the whole file into memory

.. code-block:: python

    >>> from cigsegy import SegyNP
    >>> d = SegyNP('rogan.sgy', keylocs=[9, 21])
    >>> d.shape # (ni, nx, nt), use as a numpy array, 3D geometry
    >>> sx = d[100] # the 100-th inline profile
    >>> sx = d[100:200] # return a 3D array with shape (100, nx, nt)
    >>> sx = d[:, 200, :] # the 200-th crossline profile
    >>> sx = d[:, :, 100] # the 100-th time slice, note, it may be slow if the file is large
    >>> sx.min(), sx.max() 
    # get the min and max value, but they are evaluated from a part of data, 
    # so they may not be the real min and max value
    >>> sx.ntrace # get the number of traces for the file



License
=======

cigsegy is provided under a MIT license that can be found in the `LICENSE <https://github.com/JintaoLee-Roger/cigsegy/blob/main/LICENSE>`_ file. By using, distributing, or contributing to this project, you agree to the terms and conditions of this license.


TODO
====

- Add convenient functions to support **unsorted** prestack gathers.


Citations
===========
If you find this work useful in your research and want to cite it, please consider use this:

Plain Text

.. code-block:: python

    Li, Jintao. "CIGSEGY: A tool for exchanging data between SEG-Y format and NumPy array inside Python environment". URL: https://github.com/JintaoLee-Roger/cigsegy


BibTex

.. code-block:: latex
    
    @misc{cigsegy,
    author = {Li, Jintao},
    title = {{CIGSEGY}: A tool for exchanging data between SEG-Y format and NumPy array inside Python environment},
    howpublished = {\url{https://github.com/JintaoLee-Roger/cigsegy}},
    }



=========

.. [1] Here **irregular** SEG-Y volume means the area covered by a SEG-Y file is not a rectangle but a polygon (meaning that some lines are missing some traces), or its inline/crossline intervals are not 1. 