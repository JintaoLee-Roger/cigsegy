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
#include <string>

int main(int argc, char *argv[]) {
    cxxopts::Options options(
        argv[0],
        std::string(argv[0]) + " - a tool for segy file reading to binary file"
    );

    options.add_options()
        ("i,input", "input segy file: (Required)", cxxopts::value<std::string>())
        ("o,out", "out binary file name", cxxopts::value<std::string>())
        ("f,fills", "the number to fill the miss trace, can be any float or nan, or NAN", cxxopts::value<std::string>())
        ("z,inline-loc", "inline field in trace header, default is 189", cxxopts::value<int>())
        ("c,crossline-loc", "crossline field in trace header, default is 193", cxxopts::value<int>())
        ("istep", "inline step", cxxopts::value<int>())
        ("xstep", "crossline step", cxxopts::value<int>())
        ("xloc", "X field in trace header, default is 73", cxxopts::value<int>())
        ("yloc", "Y field in trace header, default is 77", cxxopts::value<int>())
        ("p,print_textual_header", "print 3200 bytes textual header")
        ("m,meta_info", "print meta info")
        ("ignore-header", "reading segy by ignoring header and specify shape");

    options.parse_positional({"input"});
    options.add_example(std::string(argv[0]) + " -p f3.segy             : show textual header");
    options.add_example(std::string(argv[0]) + " -m f3.segy             : show meta information");
    options.add_example(std::string(argv[0]) + " -o f3.dat f3.segy      : convert");
    options.add_example(std::string(argv[0]) + " -i f3.segy -o f3.dat   : convert");
    options.add_example(std::string(argv[0]) + " -o f3.dat -z 5 f3.segy : convert by specify inline field");
    options.add_example(std::string(argv[0]) + " -o f3.dat -z 5 --istep 2 f3.segy : convert by specify inline field and step");
    options.add_example(std::string(argv[0]) + " -o f3.dat -f nan f3.segy : convert and fill with nan");
    options.add_example(std::string(argv[0]) + " -o f3.dat --ignore-header f3.segy : ignore header and specify shape");

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

    segy::SegyRW segyio(segy_name);

    if (args.count("ignore-header")) {
        if (!args.count("d")) {
            throw std::runtime_error("When using '--ignore-header', must specify shape '-d'");
        }
    }

    if (args.count("p")) {
        if (args.count("ignore-header")) {
            throw std::runtime_error("You have ignored header (--ignore-header).");
        }
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
        segyio.setYLocation(args["yloc"].as<int>()); // 确保 segy::SegyRW 类有 setYLocation 方法
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

    bool is2d = false;
    if (args.count("ignore-header")) {
        segyio.set_segy_type(2);
        is2d = true;
    } else {
        segyio.set_segy_type(3);
    }

    segyio.scan();

    if (args.count("m")) {
        if (args.count("ignore-header")) {
            throw std::runtime_error("You have ignored header (--ignore-header).");
        }
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
            std::cout << meta.start_time << ", " << meta.dt << ", " << meta.nt << "\n";
        }
    }

    if (args.count("o")) {
        std::string out_name = args["o"].as<std::string>();
        std::cout << "Write binary file to: " << out_name << "\n";
        segyio.tofile(out_name, is2d);
        segyio.close_file();
    }

    return 0;
}