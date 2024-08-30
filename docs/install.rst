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
    pip install . --config-settings editable_mode=strict


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