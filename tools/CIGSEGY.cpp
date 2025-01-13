/*********************************************************************
** Copyright (c) 2023 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
**
** @File: SEGYRead.cpp
** @Description :
*********************************************************************/

#include "cxxopts.hpp"
#include "segyrw.h"
#include "utils.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <map>


std::map<int, std::pair<std::string, int>> kBinaryHeaderHelp = {
    {1, {"Job ID", 4}},
    {5, {"Line number", 4}},
    {9, {"Reel Number", 4}},
    {13, {"N traces per ensemble", 2}},
    {15, {"N auxiliary traces per ensemble", 2}},
    {17, {"Sample interval(dt)", 2}},
    {19, {"dt of original", 2}},
    {21, {"N samples per traces(ns)", 2}},
    {23, {"ns of original", 2}},
    {25, {"Data sample format code (1-IBM, 5-IEEE)", 2}},
    {27, {"Ensemble fold", 2}},
    {29, {"Trace sorting code", 2}},
    {31, {"vertical sum code", 2}},
    {33, {"Sweep freq at start(Hz)", 2}},
    {35, {"Sweep freq at end(Hz)", 2}},
    {37, {"Sweep length(ms)", 2}},
    {39, {"Sweep type code", 2}},
    {41, {"Trace number of sweep channel", 2}},
    {43, {"Sweep trace taper length in ms at start", 2}},
    {45, {"Sweep trace taper length in ms at end", 2}},
    {47, {"Taper type", 2}},
    {49, {"Correlated data traces", 2}},
    {51, {"Binary gain recovered", 2}},
    {53, {"Amplitude recovery method", 2}},
    {55, {"Measurement system (units)", 2}},
    {57, {"Impulse signal polarity", 2}},
    {59, {"Vibratory polarity code", 2}},
    {61, {"Extended number of data traces per ensemble", 4}},
    {65, {"Extended number of auxiliary traces per ensemble", 4}},
    {69, {"Extended number of samples per data trace", 4}},
    {73, {"Extended sample interval, IEEE double precision (64-bit)", 8}},
    {81, {"Extended sample interval of original field recording, IEEE double precision (64-bit)", 8}},
    {89, {"Extended number of samples per data trace in original recording", 4}},
    {93, {"Extended ensemble fold", 4}},
    {97, {"The integer constant 16909060_10 (01020304_16)", 4}},
    {101, {"Unassigned", 200}},
    {301, {"Major SEG-Y Format Revision Number. This is an 8-bit unsigned value", 1}},
    {302, {"Minor SEG-Y Format Revision Number. This is an 8-bit unsigned value with a radix point between the first and second bytes.", 1}},
    {303, {"Fixed length trace flag", 2}},
    {305, {"Number of 3200-byte, Extended Textual File Header records following the Binary Header", 2}},
    {307, {"Max number of additional 240-byte trace header", 4}},
    {311, {"Time basis code", 2}},
    {313, {"number of trace header in this file, 64-bit unsigned integer value", 8}},
    {321, {"Byte offset of first trace relative to start of file or stream if known, otherwise zero. (64-bit unsigned integer value)", 8}},
    {329, {"Number of 3200-byte data trailer stanza records following the last trace (4 byte signed integer).", 4}},
    {333, {"Unassigned", 68}},
};

// 定义 kTraceHeaderHelp
std::map<int, std::pair<std::string, int>> kTraceHeaderHelp = {
    {1, {"Trace sequence number within line", 4}},
    {5, {"Trace sequence number within SEG-Y file", 4}},
    {9, {"Original field record number", 4}},
    {13, {"Trace number within the original field record", 4}},
    {17, {"Energy source point number", 4}},
    {21, {"Ensemble number", 4}},
    {25, {"Trace number within the ensemble", 4}},
    {29, {"Trace identification code", 2}},
    {31, {"Number of vertically summed traces yielding this trace", 2}},
    {33, {"Number of horizontally stacked traces yielding this trace", 2}},
    {35, {"Data use for", 2}},
    {37, {"Distance from center of the source point to the center of the receiver group", 4}},
    {41, {"Elevation of receiver group", 4}},
    {45, {"Surface elevation at source location", 4}},
    {49, {"Source depth below surface", 4}},
    {53, {"Seismic Datum elevation at receiver group", 4}},
    {57, {"Seismic Datum elevation at source", 4}},
    {61, {"Water column height at source location", 4}},
    {65, {"Water column height at receiver group location", 4}},
    {69, {"Scalar to be applied to all elevations and depths", 2}},
    {71, {"Scalar to be applied to all coordinates", 2}},
    {73, {"Source coordinate - X", 4}},
    {77, {"Source coordinate - Y", 4}},
    {81, {"Group coordinate - X", 4}},
    {85, {"Group coordinate - Y", 4}},
    {89, {"Coordinate units", 2}},
    {91, {"Weathering velocity", 2}},
    {93, {"Subweathering velocity", 2}},
    {95, {"Uphole time at source in ms", 2}},
    {97, {"Uphole time at group in ms", 2}},
    {99, {"Source static correction in ms", 2}},
    {101, {"Group static correction in ms", 2}},
    {103, {"Total static applied in ms", 2}},
    {105, {"Lag time A", 2}},
    {107, {"Lag Time B", 2}},
    {109, {"Delay recording time", 2}},
    {111, {"Mute time - Start time in ms", 2}},
    {113, {"Mute time - End time in ms", 2}},
    {115, {"Number of samples in this trace", 2}},
    {117, {"Sample interval for this trace", 2}},
    {119, {"Gain type of field instruments", 2}},
    {121, {"Instrument gain constant", 2}},
    {123, {"Instrument gain constant", 2}},
    {125, {"Correlated", 2}},
    {127, {"Sweep frequency at start", 2}},
    {129, {"Sweep frequency at end", 2}},
    {131, {"Sweep length in ms", 2}},
    {133, {"Sweep type", 2}},
    {135, {"Sweep trace taper length at start in ms", 2}},
    {137, {"Sweep trace taper length at end in ms", 2}},
    {139, {"Taper type", 2}},
    {141, {"Alias filter frequency (Hz)", 2}},
    {143, {"Alias filter slope (dB/octave)", 2}},
    {145, {"Notch filter frequency (Hz)", 2}},
    {147, {"Notch filter slope (dB/octave)", 2}},
    {149, {"Low-cut frequency (Hz)", 2}},
    {151, {"High-cut frequency (Hz)", 2}},
    {153, {"Low-cut slope (dB/octave)", 2}},
    {155, {"High-cut slope (dB/octave)", 2}},
    {157, {"Year data recorded", 2}},
    {159, {"Day of year", 2}},
    {161, {"Hour of day", 2}},
    {163, {"Minute of hour", 2}},
    {165, {"Second of minute", 2}},
    {167, {"Time basis code", 2}},
    {169, {"Trace weighting factor", 2}},
    {171, {"Geophone group number of roll switch position one", 2}},
    {173, {"Geophone group number of trace number one within original field record", 2}},
    {175, {"Geophone group number of last trace within original field record", 2}},
    {177, {"Gap size (total number of groups dropped)", 2}},
    {179, {"Over travel associated with taper at beginning or end of line", 2}},
    {181, {"X coordinate", 4}},
    {185, {"Y coordinate", 4}},
    {189, {"The in-line number", 4}},
    {193, {"The cross-line number", 4}},
    {197, {"Shotpoint number", 4}},
    {201, {"Scalar to be applied to the shotpoint number", 2}},
    {203, {"Trace value measurement unit", 2}},
    {205, {"Transduction Constant", 6}},
    {211, {"Transduction Units", 2}},
    {213, {"Device/Trace Identifier", 2}},
    {215, {"Scalar to be applied to bytes 95-114", 2}},
    {217, {"Source Type/Orientation", 2}},
    {219, {"Source Energy Direction with respect to the source orientation", 6}},
    {225, {"Source Measurement - Describes the source effort used to generate the trace", 6}},
    {231, {"Source Measurement Unit", 2}},
    {233, {"Either binary zeros or chars SEG00000", 8}},
};


// 辅助函数：将字节数组转换为整数
template <typename T>
T bytes_to_int(const unsigned char* data, bool big_endian = true) {
    T value = 0;
    std::memcpy(&value, data, sizeof(T));
    if (big_endian) {
        value = segy::swap_endian(value);
    }
    return value;
}

// 辅助函数：将字节数组转换为 double（IEEE 754 64-bit）
double bytes_to_double(const unsigned char* data, bool big_endian = true) {
    double value;
    std::memcpy(&value, data, sizeof(double));
    if (big_endian) {
        value = segy::swap_endian(value);
    }
    return value;
}

// 通用解析函数
void parse_header(const unsigned char* header, const std::map<int, std::pair<std::string, int>>& help_dict, const std::string& header_type) {
    std::cout << "Parsing " << header_type << " Header:\n";
    std::cout << "----------------------------------------\n";

    std::ostringstream output;
    for (const auto& [key, value] : help_dict) {
        const std::string& disc = value.first;
        int ksize = value.second;
        int64_t field_value = 0;

        if (ksize == 1) {
            field_value = header[key - 1];
        } else if (ksize == 2) {
            field_value = bytes_to_int<int16_t>(header + key - 1);
        } else if (ksize == 4) {
            field_value = bytes_to_int<int32_t>(header + key - 1);
        } else if (ksize == 8) {
            field_value = bytes_to_int<int64_t>(header + key - 1);
        } else {
            field_value = 0; // 未处理的字段
        }

        // 按照指定格式输出
        output << std::setw(3) << std::right << key << " - "
               << std::setw(3) << std::right << (key + ksize - 1) << ": "
               << std::setw(8) << std::left << field_value << " - "
               << disc << "\n";
    }
    std::cout << output.str();
}



int main(int argc, char *argv[]) {
    cxxopts::Options options(
        argv[0],
        std::string(argv[0]) + " - a tool for segy file access and conversion."
    );

    options.add_options()
        ("i,input", "input segy file: (Required)", cxxopts::value<std::string>())
        ("o,out", "out binary file name", cxxopts::value<std::string>())
        ("n,new_binary", "new binary file to create new segy file", cxxopts::value<std::string>())
        ("f,fills", "the number to fill the miss trace, can be any float or nan, or NAN", cxxopts::value<std::string>())
        ("z,inline-loc", "inline field in trace header, default is 189", cxxopts::value<int>())
        ("c,crossline-loc", "crossline field in trace header, default is 193", cxxopts::value<int>())
        ("istep", "inline step", cxxopts::value<int>())
        ("xstep", "crossline step", cxxopts::value<int>())
        ("xloc", "X field in trace header, default is 73", cxxopts::value<int>())
        ("yloc", "Y field in trace header, default is 77", cxxopts::value<int>())
        ("p,print_textual_header", "print 3200 bytes textual header")
        ("m,meta_info", "print meta info")
        ("ignore-header", "reading segy by ignoring header and specify shape")
        ("b,bheader", "show the 400 bytes binary header", cxxopts::value<bool>())
        ("t,theader", "show the header of the i-th trace", cxxopts::value<int>());

    options.parse_positional({"input"});
    options.add_example(std::string(argv[0]) + " -p f3.segy             : show textual header");
    options.add_example(std::string(argv[0]) + " -m f3.segy             : show meta information");
    options.add_example(std::string(argv[0]) + " -o f3.dat f3.segy      : convert");
    options.add_example(std::string(argv[0]) + " -i f3.segy -o f3.dat   : convert");
    options.add_example(std::string(argv[0]) + " -o f3.dat -z 5 f3.segy : convert by specify inline field");
    options.add_example(std::string(argv[0]) + " -o f3.dat -z 5 --istep 2 f3.segy : convert by specify inline field and step");
    options.add_example(std::string(argv[0]) + " -o f3.dat -f nan f3.segy : convert and fill with nan");
    options.add_example(std::string(argv[0]) + " -o f3.dat --ignore-header f3.segy : ignore header and specify shape");
    // create
    options.add_example(std::string(argv[0]) + " -i f3.segy -n new.dat -o new.segy : create new segy file from new binary");
    // show binary header
    options.add_example(std::string(argv[0]) + " -i f3.segy -b : show binary header");
    // show trace header
    options.add_example(std::string(argv[0]) + " -i f3.segy -t 100 : show the header of the 100-th trace");

    auto args = options.parse(argc, argv);

    if (argc == 1) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (!args.count("i")) {
        throw std::runtime_error("Missing input segy file");
    }

    std::string segy_name = args["i"].as<std::string>();
    std::cout << "Read segy file from: " << segy_name << std::endl;

    try {
        segy::SegyRW segyio(segy_name);

        if (args.count("p")) {
            std::cout << "Textual header:\n" << segyio.textual_header() << "\n";
        }

        if (args.count("z")) {
            segyio.setInlineLocation(args["z"].as<int>());
        }

        if (args.count("c")) {
            segyio.setCrosslineLocation(args["c"].as<int>());
        }

        if (args.count("istep")) {
            segyio.setInlineStep(args["istep"].as<int>());
        }

        if (args.count("xstep")) {
            segyio.setCrosslineStep(args["xstep"].as<int>());
        }

        if (args.count("xloc")) {
            segyio.setXLocation(args["xloc"].as<int>());
        }

        if (args.count("yloc")) {
            segyio.setYLocation(args["yloc"].as<int>());
        }

        if (args.count("f")) {
            float fills = 0.0f;
            std::string fill_str = args["f"].as<std::string>();
            if (fill_str == "nan" || fill_str == "NAN") {
                fills = NAN;
            } else {
                try {
                    fills = std::stof(fill_str);
                } catch (const std::invalid_argument&) {
                    throw std::runtime_error("Invalid fill value provided");
                }
            }
            segyio.setFill(fills);
        }

        if (args.count("b")) {
            unsigned char bheader[400];
            segyio.get_binary_header(bheader);
            parse_header(bheader, kBinaryHeaderHelp, "Binary");
        }

        if (args.count("t")) {
            std::cout << "yes\n";
            int t_n = args["t"].as<int>();
            if (t_n < 0 || t_n >= segyio.m_meta.ntrace) {
                throw std::runtime_error("Trace number out of bound");
            }
            unsigned char theader[240];
            segyio.get_trace_header(theader, t_n);
            parse_header(theader, kTraceHeaderHelp, "Trace : " + std::to_string(t_n) + " ");
        }

        bool is2d = false;
        if (args.count("ignore-header")) {
            segyio.set_segy_type(2);
            is2d = true;
        } else {
            segyio.set_segy_type(3);
        }

        segyio.scan();

        if (args.count("m")) {
            auto keys = segyio.m_keys;
            auto meta = segyio.m_meta;
            std::cout << "Meta information:\n";
            std::cout << "N traces: " << meta.ntrace << "\n";
            if (args.count("ignore-header")) {
                std::cout << "Shape (n-trace, n-time) = (" << meta.ntrace << ", " << meta.nt << ")\n";
                std::cout << "dt: " << meta.dt << ", dformat: " << meta.dformat << "\n";
            } else {
                std::cout << "Shape (n-inline, n-xline, n-time) = (" << meta.ni << ", " << meta.nx << ", " << meta.nt << ")\n";
                std::cout << "Geometry (inline, xline, time) / (start, interval, length)\n";
                std::cout << meta.start_iline << ", " << keys.istep << ", " << meta.ni << "\n";
                std::cout << meta.start_xline << ", " << keys.xstep << ", " << meta.nx << "\n";
                std::cout << meta.start_time << ", " << meta.dt / 1000 << ", " << meta.nt << "\n";
            }
        }

        if (args.count("o")) {
            if (args.count("n")) {
                std::string new_name = args["n"].as<std::string>();
                std::string out_name = args["o"].as<std::string>();
                std::cout << "Create segy file by sharing header, write segy to: " << out_name << "\n";
                std::vector<size_t> shape = segyio.shape();
                std::vector<size_t> start = {0, 0, 0};
                segyio.create_by_sharing_header(out_name, new_name, shape, start);
                segyio.close_file();
            } else {
                std::string out_name = args["o"].as<std::string>();
                std::cout << "Write binary file to: " << out_name << "\n";
                segyio.tofile(out_name, is2d);
                segyio.close_file();
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}