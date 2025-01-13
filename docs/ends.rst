Third part dependencies
#######################

**cigsegy** relies on the following third-party libraries.

- ``src/include/mio.hpp`` is from `mandreyel/mio <https://github.com/mandreyel/mio>`_: Cross-platform C++11 header-only library for memory mapped file IO

- ``tools/cxxopts.hpp`` is from `jarro2783/cxxopts <https://github.com/jarro2783/cxxopts>`_: Lightweight C++ command line option parser

- `pybind/pybind11 <https://github.com/pybind/pybind11>`_: Seamless operability between C++11 and Python


Limitations
###########

- Don't support extended textual/trace file headers
- Need to refine the unsorted file reading
- Need more high performance functions to read/write prestack data


Acknowledge
###########


Thanks to `jiarunyangustc <https://github.com/jiarunyangustc>`_ and `shenghanlin <https://github.com/shenghanlin>`_ for submitting important bugs and requests.