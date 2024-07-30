
C++ API
#######

Core functions
==============

.. doxygenfunction:: segy::read

.. doxygenfunction:: segy::tofile

.. doxygenfunction:: segy::fromfile_without_scan

.. doxygenfunction:: segy::create_by_sharing_header(const std::string &segy_name, const std::string &header_segy, const float *src, int sizeX, int sizeY, int sizeZ, int iline = 189, int xline = 193, int istep = 1, int xstep = 1, int offsetX = -1, int offsetY = -1, int offsetZ = -1, const std::vector<std::string> &custom_info = std::vector<std::string>())

.. doxygenfunction:: segy::create_by_sharing_header(const std::string &segy_name, const std::string &header_segy, const std::string &src_file, int sizeX, int sizeY, int sizeZ, int iline = 189, int xline = 193, int istep = 1, int xstep = 1, int offsetX = -1, int offsetY = -1, int offsetZ = -1, const std::vector<std::string> &custom_info = std::vector<std::string>())

.. doxygenfunction:: modify_trace_key

.. doxygenfunction:: modify_bin_key

.. doxygenfunction:: segy::disable_progressbar


SegyIO class
============

.. doxygenclass:: segy::SegyIO
    :members:
    :undoc-members:
    :private-members:



