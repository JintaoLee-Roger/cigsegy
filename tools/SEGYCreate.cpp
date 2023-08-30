/*********************************************************************
** Copyright (c) 2023 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
**
** @File: SEGYCreate.cpp
** @Description :
*********************************************************************/

#include "cxxopts.hpp"
#include "segy.h"
#include <fmt/format.h>
#include <stdexcept>
#include <vector>

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      argv[0],
      fmt::format("{} - a tool for creating a segy file from a binary file",
                  argv[0]));
  options.add_options()("i,input", "input binary file (Required)",
                        cxxopts::value<std::string>())(
      "o,out", "out segy file name (Required)", cxxopts::value<std::string>())(
      "d,dimensions",
      "the dimensions (x, y, z) or (nt, ncrossline, ninline), use as '-d "
      "128,128,256' (Required)",
      cxxopts::value<std::vector<int>>())(
      "z,inline-loc", "set inline field in trace header, default is 189",
      cxxopts::value<int>())(
      "c,crossline-loc", "set crossline field in trace header, default is 193",
      cxxopts::value<int>())("dt", "set sample interval, default is 4000",
                             cxxopts::value<int>())(
      "f,data_format",
      "data format code, 4 bytes (1 for IBM, 5 for IEEE) floating point, "
      "defualt is 5",
      cxxopts::value<int>())("dx", "set X interval", cxxopts::value<float>())(
      "dy", "set Y interval", cxxopts::value<float>())(
      "min-inline", "set start inline number", cxxopts::value<int>())(
      "min-crossline", "set start crossline number", cxxopts::value<int>())(
      "start-time", "set start time for each trace", cxxopts::value<int>())(
      "s,share", "create a segy file using a exist segy file header",
      cxxopts::value<bool>())("header", "the shared header segy file")(
      "istep", "inline step for the header segy file",
      cxxopts::value<int>())("xstep", "crossline step for the header segy file",
                             cxxopts::value<int>());

  options.parse_positional("input");
  options.add_example(fmt::format(
      "{} -i test.dat -o test.segy -d 128,128,256 : convert", argv[0]));
  options.add_example(fmt::format(
      "{} -o test.segy -d 128,128,256 test.dat : convert", argv[0]));
  options.add_example(fmt::format(
      "{} -o test.segy -d 128,128,256 -f 5 test.dat : specify data format",
      argv[0]));
  options.add_example(fmt::format("{} -o test.segy -d 128,128,256 --dt 2000 "
                                  "test.dat : specify time interval",
                                  argv[0]));
  options.add_example(fmt::format("{} -s -i test.dat -o test.segy --header "
                                  "header.segy -d 128,128,256 -z 189 -c 193 : "
                                  "using sharing header segy file",
                                  argv[0]));

  auto args = options.parse(argc, argv);

  if (argc == 1) {
    fmt::print("{}", options.help());
    exit(0);
  }
  if (!args.count("i")) {
    throw std::runtime_error("Missing input binary file");
  }
  if (!args.count("o")) {
    throw std::runtime_error("Missing out segy file");
  }
  if (!args.count("d")) {
    throw std::runtime_error(
        "Must specify the dimensions, e.g., use '-d 128,128,256'");
  }

  std::string binary_name = args["i"].as<std::string>();
  std::string segy_name = args["o"].as<std::string>();
  fmt::print("Read binary file from: {}\n", binary_name);
  fmt::print("Create segy file to: {}\n", segy_name);

  std::vector<int> dims = args["d"].as<std::vector<int>>();
  if (dims.size() != 3) {
    throw std::runtime_error(
        fmt::format("Can only create 3D data, now dimensions are: {}",
                    fmt::join(dims, ", ")));
  }

  if (args.count("s")) {
    if (!args.count("header")) {
      throw std::runtime_error(
          "--header option is necessary when -s/--share option exists");
    }
    int istep, xstep, iline, xline;
    istep = args.count("istep") ? args["istep"].as<int>() : 1;
    xstep = args.count("xstep") ? args["xstep"].as<int>() : 1;
    iline = args.count("z") ? args["z"].as<int>() : 189;
    xline = args.count("c") ? args["c"].as<int>() : 193;

    std::string header_name = args["header"].as<std::string>();
    segy::create_by_sharing_header(segy_name, header_name, binary_name, dims[0],
                                   dims[1], dims[2], iline, xline, istep,
                                   xstep);
    exit(0);
  }

  segy::SegyIO segy_create(binary_name, dims[0], dims[1], dims[2]);

  if (args.count("z")) {
    segy_create.setInlineLocation(args["z"].as<int>());
  }
  if (args.count("c")) {
    segy_create.setCrosslineLocation(args["c"].as<int>());
  }
  if (args.count("dt")) {
    segy_create.setSampleInterval(args["dt"].as<int>());
  }
  if (args.count("f")) {
    segy_create.setDataFormatCode(args["f"].as<int>());
  }
  if (args.count("dx")) {
    segy_create.setInlineInterval(args["dx"].as<float>());
  }
  if (args.count("dy")) {
    segy_create.setCrosslineInterval(args["dy"].as<float>());
  }
  if (args.count("min-inline")) {
    segy_create.setMinInline(args["min-inline"].as<int>());
  }
  if (args.count("min-crossline")) {
    segy_create.setMinCrossline(args["min-crossline"].as<int>());
  }
  if (args.count("start-time")) {
    segy_create.setStartTime(args["start-time"].as<int>());
  }

  segy_create.create(segy_name);
}