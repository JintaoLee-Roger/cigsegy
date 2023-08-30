/*********************************************************************
** Copyright (c) 2023 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
**
** @File: SEGYRead.cpp
** @Description :
*********************************************************************/

#include "cxxopts.hpp"
#include "segy.h"
#include <cmath>
#include <fmt/format.h>
#include <stdexcept>

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      argv[0],
      fmt::format("{} - a tool for segy file reading to binary file", argv[0]));
  options.add_options()("i,input", "input segy file: (Required)",
                        cxxopts::value<std::string>())(
      "o,out", "out binary file name", cxxopts::value<std::string>())(
      "f,fills",
      "the number to fill the miss trace, can be any float or nan, or NAN",
      cxxopts::value<std::string>())(
      "z,inline-loc", "inline field in trace header, default is 189",
      cxxopts::value<int>())(
      "c,crossline-loc", "crossline field in trace header, default is 193",
      cxxopts::value<int>())("istep", "inline step", cxxopts::value<int>())(
      "xstep", "crossline step", cxxopts::value<int>())(
      "xloc", "X field in trace header, default is 73", cxxopts::value<int>())(
      "yloc", "Y field in trace header, default is 77", cxxopts::value<int>())(
      "p,print_textual_header",
      "print 3200 bytes textual header")("m,meta_info", "print meta info")(
      "ignore-header", "reading segy by ignoring header and specify shape")(
      "d,dimensions",
      "the dimensions (x, y, z) or (nt, ncrossline, ninline), use as '-d "
      "128,128,256' (Required)",
      cxxopts::value<std::vector<int>>());

  options.parse_positional("input");
  options.add_example(
      fmt::format("{} -p f3.segy             : show textual header", argv[0]));
  options.add_example(fmt::format(
      "{} -m f3.segy             : show meta information", argv[0]));
  options.add_example(
      fmt::format("{} -o f3.dat f3.segy      : convert", argv[0]));
  options.add_example(
      fmt::format("{} -i f3.segy -o f3.dat   : convert", argv[0]));
  options.add_example(fmt::format(
      "{} -o f3.dat -z 5 f3.segy : convert by specify inline field", argv[0]));
  options.add_example(fmt::format("{} -o f3.dat -z 5 --istep 2 f3.segy : "
                                  "convert by specify inline field and step",
                                  argv[0]));
  options.add_example(fmt::format(
      "{} -o f3.dat -f nan f3.segy : convert and fill with nan", argv[0]));
  options.add_example(fmt::format("{} -o f3.dat --ignore-header -d 236,789,890 "
                                  "f3.segy : ignore header and specify shape",
                                  argv[0]));

  auto args = options.parse(argc, argv);
  if (argc == 1) {
    fmt::print("{}", options.help());
    exit(0);
  }
  if (!args.count("i")) {
    throw std::runtime_error("Missing input segy file");
  }

  std::string segy_name = args["i"].as<std::string>();
  fmt::print("Read segy file from: {}\n", segy_name);

  segy::SegyIO segyio(segy_name);

  if (args.count("ignore-header")) {
    if (!args.count("d")) {
      throw std::runtime_error(
          "when using '--ignore-header', must specify shape '-d'");
    }
  }

  if (args.count("p")) {
    if (args.count("ignore-header")) {
      throw std::runtime_error("You have ignored header (--ignore-header).");
    }
    fmt::print("textual header:\n{}\n", segyio.textual_header());
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
    segyio.setXLocation(args["yloc"].as<int>());
  }

  if (args.count("f")) {
    float fills = 0;
    if (args["f"].as<std::string>() == "nan" ||
        args["f"].as<std::string>() == "NAN") {
      fills = NAN;
    }
    segyio.setFillNoValue(fills);
  }

  if (args.count("m")) {
    if (args.count("ignore-header")) {
      throw std::runtime_error("You have ignored header (--ignore-header).");
    }
    fmt::print("meta information: \n{}\n", segyio.metaInfo());
  }

  if (args.count("o")) {
    std::string out_name = args["o"].as<std::string>();
    fmt::print("Write binary file to: {}\n", out_name);
    if (args.count("ignore-header")) {
      std::vector<int> dims = args["d"].as<std::vector<int>>();
      if (dims.size() != 3) {
        throw std::runtime_error(
            fmt::format("Can only create 3D data, now dimensions are: {}",
                        fmt::join(dims, ", ")));
      }
      segyio.set_size(dims[0], dims[1], dims[2]);
    }
    segyio.tofile(out_name);
    segyio.close_file();
  }
}