Installing
##########

There are several ways to install ``cigsegy``, whose source code 
lives at `JintaoLee-Roger/cigsegy on Github <https://github.com/JintaoLee-Roger/cigsegy>`_.
I recommand using PyPi for installing it.

Install via PyPi
-----------------

You can install cigsegy via PyPi, just use:

.. code-block:: bash

    pip install cigsegy


Install Locally
---------------

If you want to install cigsegy locally, you can use the following commands:

.. code-block:: bash

    pip install -U pip
    pip install -e . --config-settings editable_mode=strict


If you need the executable files, please see ``tools/run.sh``:

.. code-block:: bash

    # Install on Linux/MacOS
    clang++ -std=c++17 -o CIGSEGY CIGSEGY.cpp ../cigsegy/cpp/segyrw.cpp -I../cigsegy/cpp

    # Install on Windows
    # please use "x64" Native Tools Command Prompt for VS 2019
    # if use "x86" toolchain, you may fail to deal with SEG-Y files larger than 2GB.
    cl /EHsc /std:c++17 /utf-8 /O2 /I"tools" /I"cigsegy\cpp" tools\CIGSEGY.cpp cigsegy\cpp\segyrw.cpp /Fe:CIGSEGY.exe