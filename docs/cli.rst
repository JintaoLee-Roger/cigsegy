cigsegy - CLI
#############

CIGSEGY offers two executables to inspect SEG-Y files on the command line.

You can build locally via cmake (can be found in ``your-install-dir/bin/``),
or can download them in `github releases <https://github.com/JintaoLee-Roger/cigsegy/releases>`_.

- ``SEGYRead``

.. code-block:: bash

    $ SEGYRead
    # SEGYRead - a tool for segy file reading to binary file
    # Usage:
    # ./SEGYRead [OPTION...] positional parameters

    # -o, --out arg               out binary file name
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
    # -d, --dimensions arg        the dimensions (x, y, z) or (nt, ncrossline, 
    #                             ninline), use as '-d 128,128,256' (Required)

    # Examples:
    # ./SEGYRead -p f3.segy             : show textual header
    # ./SEGYRead -m f3.segy             : show meta information
    # ./SEGYRead -o f3.dat f3.segy      : convert
    # ./SEGYRead -i f3.segy -o f3.dat   : convert
    # ./SEGYRead -o f3.dat -z 5 f3.segy : convert by specify inline field
    # ./SEGYRead -o f3.dat -z 5 --istep 2 f3.segy : convert by specify inline field and step
    # ./SEGYRead -o f3.dat -f nan f3.segy : convert and fill with nan
    # ./SEGYRead -o f3.dat --ignore-header -d 236,789,890 f3.segy : ignore header and specify shape



- ``SEGYCeate``

.. code-block:: bash

    $ SEGYCreate
    # SEGYCreate - a tool for creating a segy file from a binary file
    # Usage:
    # ./SEGYCreate [OPTION...] positional parameters

    # -o, --out arg            out segy file name (Required)
    # -d, --dimensions arg     the dimensions (x, y, z) or (nt, ncrossline, 
    #                         ninline), use as '-d 128,128,256' (Required)
    # -z, --inline-loc arg     set inline field in trace header, default is 189
    # -c, --crossline-loc arg  set crossline field in trace header, default is 
    #                         193
    #     --dt arg             set sample interval, default is 4000
    # -f, --data_format arg    data format code, 4 bytes (1 for IBM, 5 for 
    #                         IEEE) floating point, defualt is 5
    #     --dx arg             set X interval
    #     --dy arg             set Y interval
    #     --min-inline arg     set start inline number
    #     --min-crossline arg  set start crossline number
    #     --start-time arg     set start time for each trace
    # -s, --share              create a segy file using a exist segy file 
    #                         header
    #     --header             the shared header segy file
    #     --istep arg          inline step for the header segy file
    #     --xstep arg          crossline step for the header segy file

    # Examples:
    # ./SEGYCreate -i test.dat -o test.segy -d 128,128,256 : convert
    # ./SEGYCreate -o test.segy -d 128,128,256 test.dat : convert
    # ./SEGYCreate -o test.segy -d 128,128,256 -f 5 test.dat : specify data format
    # ./SEGYCreate -o test.segy -d 128,128,256 --dt 2000 test.dat : specify time interval
    # ./SEGYCreate -s -i test.dat -o test.segy --header header.segy -d 128,128,256 -z 189 -c 193 : using sharing header segy file