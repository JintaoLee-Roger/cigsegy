
Changelog
#########



v1.1.9
--------
- fixed a bug with out-of-int32 range offsets (occurred on Windows)
- changed code for creating large files to speed things up
- optimize ``SegyNP``



v1.1.8
---------
- fixed a bug when create file on windows
- deal with interreput signal
- optimize ``SegyNP``



v1.1.7
--------

- Refine the ``scan`` function to support more situations.
- Add supports for dealing with more data sample formats, such as 4-byte, two's complement integer.
- Add a new class ``SegyNP`` to simulate the segy file accessed as a numpy array.
- Add functions: ``scan_unsorted3D`` and ``load_unsorted3D`` to support 3D unsorted data.
- Remove the comparison part of the documents, as ``segysak`` has a large update.
- ``use_guess`` in the functions like ``metaInfo`` has been deseperated.
- Added more atomic operations, enabling finer control of SEG-Y files
- And more ...


v1.1.6
-------

- Fixed a bug of ``create_by_sharing_header`` when the file is huge.
- Add functions: ``scan_prestack`` and ``load_prestack3D`` to support 3D prestack data (4D array).


v1.1.5
------

- Disable the progress bar in jupyter notbook.
- Allow ``numpy slice`` in ``create_by_sharing_header`` function. Now you **don't** need ``data.clone()`` to create a memory continuous sub-array to the function.
- Add a ``tools.read_header`` to read the full binary/trace header. This function is based on 3 functions of ``Pysegy``, i.e., ``Pysegy.get_binary_header()``, ``Pysegy.get_trace_header()``, ``Pysegy.get_trace()``.
- Add a function ``tools.get_metaInfo()`` obtain the meta information in dict.
- Add documents and post on `read the docs <https://cigsegy.readthedocs.io/>`_ website.

v1.1.4
------

- Support create a segy from a sub-volume array and a full volume segy header, using ``cigsegy.create_by_sharing_header('sub.segy', 'header-full.segy', dsub, iline=, xline=, istep=, xstep=, offset=)``
- Add a function ``Pysegy.cut(outsegy, ...)`` to cut a sub-volume from the full segy file and save to a segy file
- For all functions to create a segy file, you can use ``custom_info`` to custom the first 12 lines of textual header
- Support the segy files with inline first order, (i.e., the order is time->inline->crossline, the default is time->crossline->inline)
- Fixed a bug of ``scan()`` function
- Add a function ``tools.plot_region(...)`` to plot the region of a segy file

v1.1.2
------

- Add some useful functions, such as ``trace_count()``, ``get_traceInfo()``
- ``fromfile``, ``tofile``, ``create_by_sharing_header`` functions support ``istep`` and ``xstep`` parameters to deal with the segy whose steps of inline/crossline is not 1
- Add ``tools.guess`` function to guess the locations and steps of inline/crossline
- Add ``tools.fromfile_by_guess``, ``tools.tofile_by_guess``, ``create_by_sharing_header_guess`` to support unknown segy files
- Support windows/MSVC