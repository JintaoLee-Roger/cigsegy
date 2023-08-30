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


First, you need to install ``fmt`` and ``pybind11``.

.. code-block:: bash

    # linux
    sudo apt-get install python3-pybind11 libfmt-dev

    # mac
    brew install pybind11 fmt


You can also install ``fmt`` and ``pybind11`` manually.

.. code-block:: bash

    # install fmt
    mkdir thridPart && cd thridPart/
    git clone https://github.com/fmtlib/fmt.git
    # Now fmt is installed into /xxx/cigsegy/thridPart/fmt

    # install pybind11 using pypi
    pip install pybind11
    # Now pybind11 is installed into /xxx/lib/python3.8/site-packages/pybind11/


If you need the python package only, run:

.. code-block:: bash

    pip install -U pip
    pip install .

    # if fmt is not in the env path
    pip install . --install-option="--fmt_root=/xxx/fmt"

    # if you need to build a wheel
    # pip wheel . --build-option="--fmt_root=/xxx/fmt"


If you need the two executables files:

.. code-block:: bash

    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=your-install-dir
    make
    make install

    # if install fmt manually and pybind11 via pypi, using this
    cmake .. -DCMAKE_INSTALL_PREFIX=your-install-dir \
    -Dpybind11_ROOT=/xxx/lib/python3.8/site-packages/pybind11/ \
    -Dfmt_ROOT=/xxx/cigsegy/thridPart/fmt/build/fmt/

    # selecting other python version, add -DPYTHON_EXECUTABL
    cmake .. -DCMAKE_INSTALL_PREFIX=your-install-dir -DPYTHON_EXECUTABLE=/xxx/bin/python -DPYTHON_LIBRARIES=/xxx/lib/