cigsegy CLI
###########

``cigsegy`` installs a Python command-line interface together with the Python
package.  The command is a thin wrapper around the public Python API, so the
behavior stays aligned with ``textual_header``, ``metaInfo``, ``fromfile``,
``collect``, ``tofile``, and ``SegyWriter``.

Run:

.. code-block:: bash

   cigsegy --help

Every subcommand includes common examples and shows default values in its help
message.


Common Workflow
===============

Print the textual header:

.. code-block:: bash

   cigsegy textual input.sgy
   cigsegy textual input.sgy --coding e --output textual.txt

Scan geometry metadata:

.. code-block:: bash

   cigsegy meta input.sgy
   cigsegy meta gather.sgy

If the inferred geometry is wrong or ambiguous, pass the known header fields:

.. code-block:: bash

   cigsegy meta input.sgy --iline 189 --xline 193
   cigsegy meta gather.sgy --iline 189 --xline 193 --offset 37 --is4d

Inspect binary or trace headers:

.. code-block:: bash

   cigsegy header input.sgy --binary
   cigsegy header input.sgy --trace 100
   cigsegy header input.sgy --trace 100 --output trace100.json

Read a regular SEG-Y volume into NumPy:

.. code-block:: bash

   cigsegy fromfile input.sgy volume.npy
   cigsegy fromfile gather.sgy gathers.npy

``fromfile`` materializes the volume in memory before saving it.  For data that
is too large for RAM, stream directly to ``.npy``:

.. code-block:: bash

   cigsegy to-npy input.sgy volume.npy
   cigsegy to-npy gather.sgy gathers.npy

Explicit geometry is still available when needed:

.. code-block:: bash

   cigsegy fromfile input.sgy volume.npy --iline 189 --xline 193
   cigsegy fromfile gather.sgy gathers.npy --iline 189 --xline 193 --offset 37 --is4d
   cigsegy to-npy input.sgy volume.npy --iline 189 --xline 193

Collect traces in file order:

.. code-block:: bash

   cigsegy collect input.sgy traces.npy
   cigsegy collect input.sgy traces.npy --beg 0 --end 1000 --tbeg 100 --tend 800
   cigsegy collect input.sgy traces.npy --indices 0,10,20,30

Export raw IEEE float32 samples:

.. code-block:: bash

   cigsegy tofile input.sgy samples.dat
   cigsegy tofile input.sgy samples.dat --as2d

Raw output has no shape metadata.  Use ``to-npy`` when you want a NumPy file
that can be reopened with ``numpy.load(..., mmap_mode="r")``.

Write SEG-Y from a template:

.. code-block:: bash

   cigsegy create template.sgy processed.npy processed.sgy --overwrite
   cigsegy create template.sgy processed.dat processed.sgy --shape 100,200,300

Pass explicit key locations when the template geometry cannot be inferred:

.. code-block:: bash

   cigsegy create template.sgy processed.npy processed.sgy --iline 189 --xline 193 --overwrite

Change the time sample interval when writing back:

.. code-block:: bash

   cigsegy create input_2ms.sgy output_1ms.npy output_1ms.sgy \
       --sample-interval-us 1000 --non-strict

Thin a regular template by logical axes:

.. code-block:: bash

   cigsegy create input.sgy thin.npy thin.sgy \
       --iline-slice ::2 --xline-slice ::2 --sample-slice ::2


Native C++ Tool
===============

The repository still keeps ``tools/CIGSEGY.cpp`` as an optional standalone C++
program for project-specific or deployment-specific use.  It is not built
automatically during Python package installation.  Prefer the Python
``cigsegy`` command for normal package users; build the native tool manually
only when a specific project needs a small external executable.


Build Manually
--------------

Linux or macOS:

.. code-block:: bash

   clang++ -std=c++17 -O3 \
       -o CIGSEGY \
       tools/CIGSEGY.cpp cigsegy/cpp/segyrw.cpp \
       -Itools -Icigsegy/cpp

Windows, from a 64-bit MSVC developer prompt:

.. code-block:: bat

   cl /EHsc /std:c++17 /utf-8 /O2 ^
      /I"tools" /I"cigsegy\cpp" ^
      tools\CIGSEGY.cpp cigsegy\cpp\segyrw.cpp ^
      /Fe:CIGSEGY.exe


Native Examples
---------------

Print the textual header:

.. code-block:: bash

   CIGSEGY -p input.sgy

Print metadata:

.. code-block:: bash

   CIGSEGY -m input.sgy

Convert samples to a raw binary file:

.. code-block:: bash

   CIGSEGY -o output.dat input.sgy

Specify header locations:

.. code-block:: bash

   CIGSEGY -o output.dat -z 189 -c 193 --istep 1 --xstep 1 input.sgy

Inspect headers:

.. code-block:: bash

   CIGSEGY -b input.sgy
   CIGSEGY -t 100 input.sgy

Create a SEG-Y from a new binary file while sharing an existing header:

.. code-block:: bash

   CIGSEGY -i template.sgy -n new.dat -o new.sgy
