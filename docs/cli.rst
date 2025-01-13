cigsegy - CLI
#############

CIGSEGY offers a executable to inspect SEG-Y files on the command line.

You can build locally,
or can download them in `github releases <https://github.com/JintaoLee-Roger/cigsegy/releases>`_.

- ``CIGSEGY``

.. code-block:: bash

    $ SEGYRead
    # D:\cigsegy\CIGSEGY.exe - a tool for segy file access and conversion.
    # Usage:
    # D:\cigsegy\CIGSEGY.exe [OPTION...] positional parameters

    # -o, --out arg               out binary file name
    # -n, --new_binary arg        new binary file to create new segy file
    # -f, --fills arg             the number to fill the miss trace, can be any
    #                             float or nan, or NAN
    # -z, --inline-loc arg        inline field in trace header, default is 189
    # -c, --crossline-loc arg     crossline field in trace header, default is
    #                             193
    #     --istep arg             inline step
    #     --xstep arg             crossline step
    #     --xloc arg              X field in trace header, default is 73
    #     --yloc arg              Y field in trace header, default is 77
    # -p, --print_textual_header  print 3200 bytes textual header
    # -m, --meta_info             print meta info
    #     --ignore-header         reading segy by ignoring header and specify
    #                             shape
    # -b, --bheader               show the 400 bytes binary header
    # -t, --theader arg           show the header of the i-th trace

    # Examples:
    # D:\cigsegy\CIGSEGY.exe -p f3.segy             : show textual header
    # D:\cigsegy\CIGSEGY.exe -m f3.segy             : show meta information
    # D:\cigsegy\CIGSEGY.exe -o f3.dat f3.segy      : convert
    # D:\cigsegy\CIGSEGY.exe -i f3.segy -o f3.dat   : convert
    # D:\cigsegy\CIGSEGY.exe -o f3.dat -z 5 f3.segy : convert by specify inline field
    # D:\cigsegy\CIGSEGY.exe -o f3.dat -z 5 --istep 2 f3.segy : convert by specify inline field and step
    # D:\cigsegy\CIGSEGY.exe -o f3.dat -f nan f3.segy : convert and fill with nan
    # D:\cigsegy\CIGSEGY.exe -o f3.dat --ignore-header f3.segy : ignore header and specify shape
    # D:\cigsegy\CIGSEGY.exe -i f3.segy -n new.dat -o new.segy : create new segy file from new binary
    # D:\cigsegy\CIGSEGY.exe -i f3.segy -b : show binary header
    # D:\cigsegy\CIGSEGY.exe -i f3.segy -t 100 : show the header of the 100-th trace
