/*********************************************************************
** Copyright (c) 2023 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
**
** @File: PySegy.cpp
** @Description :
*********************************************************************/

#include "segy.h"
#include "utils.h"
#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>

namespace py = pybind11;
using npfloat = py::array_t<float, py::array::c_style | py::array::forcecast>;

class Pysegy : public segy::SegyIO {
public:
  using segy::SegyIO::create;
  using segy::SegyIO::cut;
  using segy::SegyIO::read;
  using segy::SegyIO::read_cross_slice;
  using segy::SegyIO::read_inline_slice;
  using segy::SegyIO::read_time_slice;
  using segy::SegyIO::read_trace;
  using segy::SegyIO::SegyIO;

  py::array_t<int> get_lineInfo();

  py::array_t<uchar> get_binary_header(bool raw = false);
  py::array_t<uchar> get_trace_header(int n, bool raw = false);
  py::array_t<uchar> get_trace(int n, bool raw = false);

  py::array_t<float> read(int startZ, int endZ, int startY, int endY,
                          int startX, int endX);
  py::array_t<float> read();
  py::array_t<float> read_inline_slice(int iZ);
  py::array_t<float> read_cross_slice(int iY);
  py::array_t<float> read_time_slice(int iX);
  py::array_t<float> read_trace(int iZ, int iY);

  void cut(const std::string &outname, int startZ, int endZ, int startY,
           int endY, int startX, int endX,
           const py::list &custom_info = py::list());
  void cut(const std::string &outname, int startZ, int endZ, int startY,
           int endY, const py::list &custom_info = py::list());
  void cut(const std::string &outname, int startX, int endX,
           const py::list &custom_info = py::list());

  void create(const std::string &segy_out_name, const npfloat &src,
              const py::list &custom_info = py::list());
  void create(const std::string &segy_out_name,
              const py::list &custom_info = py::list());
};

py::array_t<int> Pysegy::get_lineInfo() {
  // [inline, crossline_start, crossline_end,
  //  trace_strat, trace_end, count]
  auto info = line_info();
  int n = info.size();
  py::array_t<int> out({n, 6});
  int *outptr = static_cast<int *>(out.request().ptr);
  memcpy(outptr, info.data(), sizeof(int) * 6 * n);
  if (!is_crossline_fast_order()) {
    std::cout
        << "[Warining] As the fast order of your segy file "
           "is inline order (default is crossline order). So the lineInfo "
           "you obtained is seems as [crossline, inline_start, inline_end, "
           "trace_start, trace_end, count]\n"
        << std::endl;
  }
  return out;
}

py::array_t<uchar> Pysegy::get_binary_header(bool raw) {
  py::array_t<uchar> out(segy::kBinaryHeaderSize);
  uchar *outptr = static_cast<uchar *>(out.request().ptr);
  get_binary_header_full(outptr, raw);
  return out;
}

py::array_t<uchar> Pysegy::get_trace_header(int n, bool raw) {
  py::array_t<uchar> out(segy::kTraceHeaderSize);
  uchar *outptr = static_cast<uchar *>(out.request().ptr);
  get_trace_header_full(n, outptr, raw);
  return out;
}

py::array_t<uchar> Pysegy::get_trace(int n, bool raw) {
  int sizeX = shape(0);
  py::array_t<uchar> out(segy::kTraceHeaderSize + sizeX * sizeof(float));
  uchar *outptr = static_cast<uchar *>(out.request().ptr);
  get_trace_full(n, outptr, raw);
  return out;
}

// Be careful the order of the dimensions
// In segy, X (time) is the first, but when you read into python,
// Z (inline) is the first. If need change it to X first,
// use data.transpose() in python
py::array_t<float> Pysegy::read(int startZ, int endZ, int startY, int endY,
                                int startX, int endX) {

  if (startX >= endX || startY >= endY || startZ >= endZ) {
    throw std::runtime_error("Index 'end' must large than 'start'");
  }
  if (startX < 0 || endX > shape(0) || startY < 0 || endY > shape(1) ||
      startZ < 0 || endZ > shape(2)) {
    throw std::runtime_error("Index out of range");
  }

  int sizeX = endX - startX;
  int sizeY = endY - startY;
  int sizeZ = endZ - startZ;
  py::array_t<float> out({sizeZ, sizeY, sizeX});
  auto buff = out.request();
  float *ptr = static_cast<float *>(buff.ptr);
  read(ptr, startX, endX, startY, endY, startZ, endZ);
  return out;
}

py::array_t<float> Pysegy::read() {
  py::array_t<float> out({shape(2), shape(1), shape(0)});
  auto buff = out.request();
  float *ptr = static_cast<float *>(buff.ptr);
  read(ptr);
  return out;
}
py::array_t<float> Pysegy::read_inline_slice(int iZ) {
  py::array_t<float> out({shape(1), shape(0)});
  auto buff = out.request();
  float *ptr = static_cast<float *>(buff.ptr);
  read_inline_slice(ptr, iZ);
  return out;
}
py::array_t<float> Pysegy::read_cross_slice(int iY) {
  py::array_t<float> out({shape(2), shape(0)});
  auto buff = out.request();
  float *ptr = static_cast<float *>(buff.ptr);
  read_cross_slice(ptr, iY);
  return out;
}
py::array_t<float> Pysegy::read_time_slice(int iX) {
  py::array_t<float> out({shape(2), shape(1)});
  auto buff = out.request();
  float *ptr = static_cast<float *>(buff.ptr);
  read_time_slice(ptr, iX);
  return out;
}

py::array_t<float> Pysegy::read_trace(int iZ, int iY) {
  py::array_t<float> out(shape(0));
  auto buff = out.request();
  float *ptr = static_cast<float *>(buff.ptr);
  read_trace(ptr, iY, iZ);
  return out;
}

void Pysegy::cut(const std::string &outname, int startZ, int endZ, int startY,
                 int endY, int startX, int endX, const py::list &custom_info) {
  auto customvec = custom_info.cast<std::vector<std::string>>();
  cut(outname, startX, endX, startY, endY, startZ, endZ, customvec);
}

void Pysegy::cut(const std::string &outname, int startZ, int endZ, int startY,
                 int endY, const py::list &custom_info) {
  auto customvec = custom_info.cast<std::vector<std::string>>();
  cut(outname, startY, endY, startZ, endZ, customvec);
}

void Pysegy::cut(const std::string &outname, int startX, int endX,
                 const py::list &custom_info) {
  auto customvec = custom_info.cast<std::vector<std::string>>();
  cut(outname, startX, endX, customvec);
}

void Pysegy::create(const std::string &segy_out_name, const npfloat &src,
                    const py::list &custom_info) {
  auto buff = src.request();
  if (buff.ndim != 3) {
    throw std::runtime_error("Input data must be a 3D data.");
  }
  auto r = src.shape();
  set_size(r[2], r[1], r[0]);
  auto customvec = custom_info.cast<std::vector<std::string>>();

  float *ptr = static_cast<float *>(buff.ptr);
  create(segy_out_name, ptr, customvec);
}

void Pysegy::create(const std::string &segy_out_name,
                    const py::list &custom_info) {
  auto customvec = custom_info.cast<std::vector<std::string>>();
  create(segy_out_name, customvec);
}

void create_by_sharing_header(const std::string &segy_name,
                              const std::string &header_segy,
                              const npfloat &src, int iline = 189,
                              int xline = 193, int istep = 1, int xstep = 1,
                              const py::object &offset = py::none(),
                              const py::list &custom_info = py::list()) {
  auto buff = src.request();
  if (buff.ndim != 3) {
    throw std::runtime_error("Input data must be a 3D data.");
  }

  auto r = src.shape();
  float *ptr = static_cast<float *>(buff.ptr);

  if (offset.is_none()) {
    segy::create_by_sharing_header(segy_name, header_segy, ptr, r[2], r[1],
                                   r[0], iline, xline, istep, xstep);
  } else {
    int offsetX, offsetY, offsetZ;
    if (py::isinstance<py::dict>(offset)) {
      py::dict off = offset.cast<py::dict>();
      offsetZ = off["iline"].cast<int>();
      offsetY = off["xline"].cast<int>();
      offsetX = off["time"].cast<int>();
    } else if (py::isinstance<py::sequence>(offset)) {
      py::sequence off = offset.cast<py::sequence>();
      offsetX = py::cast<int>(off[2]);
      offsetY = py::cast<int>(off[1]);
      offsetZ = py::cast<int>(off[0]);
    } else {
      throw std::runtime_error("Unkown type of offset");
    }
    auto customvec = custom_info.cast<std::vector<std::string>>();

    segy::create_by_sharing_header(segy_name, header_segy, ptr, r[2], r[1],
                                   r[0], iline, xline, istep, xstep, offsetX,
                                   offsetY, offsetZ, customvec);
  }
}

void create_by_sharing_header(const std::string &segy_name,
                              const std::string &header_segy,
                              const std::string &src_file,
                              const py::sequence &shape, int iline = 189,
                              int xline = 193, int istep = 1, int xstep = 1,
                              const py::object &offset = py::none(),
                              const py::list &custom_info = py::list()) {
  if (shape.size() != 3) {
    throw std::runtime_error("dimensions must be 3");
  }

  if (offset.is_none()) {
    segy::create_by_sharing_header(
        segy_name, header_segy, src_file, py::cast<int>(shape[2]),
        py::cast<int>(shape[1]), py::cast<int>(shape[0]), iline, xline, istep,
        xstep);
  } else {
    int offsetX, offsetY, offsetZ;
    if (py::isinstance<py::dict>(offset)) {
      py::dict off = offset.cast<py::dict>();
      offsetZ = off["iline"].cast<int>();
      offsetY = off["xline"].cast<int>();
      offsetX = off["time"].cast<int>();
    } else if (py::isinstance<py::sequence>(offset)) {
      py::sequence off = offset.cast<py::sequence>();
      offsetX = py::cast<int>(off[2]);
      offsetY = py::cast<int>(off[1]);
      offsetZ = py::cast<int>(off[0]);
    } else {
      throw std::runtime_error("Unkown type of offset");
    }
    auto customvec = custom_info.cast<std::vector<std::string>>();

    segy::create_by_sharing_header(
        segy_name, header_segy, src_file, py::cast<int>(shape[2]),
        py::cast<int>(shape[1]), py::cast<int>(shape[0]), iline, xline, istep,
        xstep, offsetX, offsetY, offsetZ, customvec);
  }
}

py::array_t<float> fromfile_ignore_header(const std::string &segy_name,
                                          int sizeZ, int sizeY, int sizeX,
                                          int format = 5) {
  py::array_t<float> out({sizeZ, sizeY, sizeX});
  auto buff = out.request();
  float *ptr = static_cast<float *>(buff.ptr);
  segy::read_ignore_header(segy_name, ptr, sizeX, sizeY, sizeZ, format);
  return out;
}

py::array_t<float> fromfile(const std::string &segy_name, int iline = 189,
                            int xline = 193, int istep = 1, int xstep = 1) {
  Pysegy segy_data(segy_name);
  segy_data.setInlineLocation(iline);
  segy_data.setCrosslineLocation(xline);
  segy_data.setSteps(istep, xstep);
  segy_data.scan();
  py::array_t<float> out = segy_data.read();

  return out;
}

py::array_t<float> collect(const std::string &segy_in, int beg = -1,
                           int end = 0) {
  Pysegy segy(segy_in);

  int need = 0;
  if (beg < 0) {
    need = segy.trace_count();
  } else {
    if (end == 0) {
      need = 1;
    } else if (end < 0) {
      need = segy.trace_count() - beg;
    } else {
      need = end - beg;
    }
  }

  auto data = py::array_t<float>({need, segy.shape(0)});
  auto buff = data.request();
  float *ptr = static_cast<float *>(buff.ptr);

  segy.collect(ptr, beg, end);

  return data;
}

py::array_t<float> _load_prestack3D(const std::string &segy_name,
                                    const py::sequence &shape, int min_iline,
                                    int max_iline, int min_xline, int max_xline,
                                    int min_offset, int max_offset, int istep,
                                    int xstep, int ostep, int iline, int xline,
                                    int offset = 37, float fill = 0) {
  auto out =
      py::array_t<float>({py::cast<int>(shape[0]), py::cast<int>(shape[1]),
                          py::cast<int>(shape[2]), py::cast<int>(shape[3])});
  float *ptr = static_cast<float *>(out.request().ptr);
  segy::load_prestack3D(ptr, segy_name, py::cast<int>(shape[3]), min_iline,
                        max_iline, min_xline, max_xline, min_offset, max_offset,
                        istep, xstep, ostep, iline, xline, offset, fill);

  return out;
}

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(cigsegy, m) {
  py::class_<Pysegy>(m, "Pysegy")
      .def(py::init<std::string>())
      .def(py::init<int, int, int>())
      .def(py::init<std::string, int, int, int>())
      .def("get_lineInfo", &Pysegy::get_lineInfo)
      .def("get_metaInfo", &Pysegy::get_metaInfo)
      .def("get_binary_header", &Pysegy::get_binary_header,
           py::arg("raw") = false)
      .def("get_trace_header", &Pysegy::get_trace_header, py::arg("n"),
           py::arg("raw") = false)
      .def("get_trace", &Pysegy::get_trace, py::arg("n"),
           py::arg("raw") = false)
      .def("setInlineLocation", &Pysegy::setInlineLocation, py::arg("iline"))
      .def("setCrosslineLocation", &Pysegy::setCrosslineLocation,
           py::arg("xline"))
      .def("setXLocation", &Pysegy::setXLocation, py::arg("xfield"))
      .def("setYLocation", &Pysegy::setYLocation, py::arg("yfield"))
      .def("setInlineStep", &Pysegy::setInlineStep, py::arg("step"))
      .def("setCrosslineStep", &Pysegy::setCrosslineStep, py::arg("step"))
      .def("setSteps", &Pysegy::setSteps, py::arg("istep"), py::arg("xstep"))
      .def("setFillNoValue", &Pysegy::setFillNoValue, py::arg("fills"))
      .def("scan", &Pysegy::scan)
      .def("tofile", &Pysegy::tofile, py::arg("binary_out_name"))
      .def("read", overload_cast_<>()(&Pysegy::read), "read hole volume")
      .def("read",
           overload_cast_<int, int, int, int, int, int>()(&Pysegy::read),
           "read with index", py::arg("startZ"), py::arg("endZ"),
           py::arg("startY"), py::arg("endY"), py::arg("startX"),
           py::arg("endX"))
      .def("read_inline_slice",
           overload_cast_<int>()(&Pysegy::read_inline_slice),
           "read inline slice", py::arg("iZ"))
      .def("read_cross_slice", overload_cast_<int>()(&Pysegy::read_cross_slice),
           "read crossline slice", py::arg("iY"))
      .def("read_time_slice", overload_cast_<int>()(&Pysegy::read_time_slice),
           "read time slice", py::arg("iX"))
      .def("read_trace", overload_cast_<int, int>()(&Pysegy::read_trace),
           "read trace", py::arg("iZ"), py::arg("iY"))
      .def("cut",
           overload_cast_<const std::string &, int, int, int, int, int, int,
                          const py::list &>()(&Pysegy::cut),
           "cut a sub volume to a segy", py::arg("outname"), py::arg("startZ"),
           py::arg("endZ"), py::arg("startY"), py::arg("endY"),
           py::arg("startX"), py::arg("endX"),
           py::arg("custom_info") = py::list())
      .def("cut",
           overload_cast_<const std::string &, int, int, int, int,
                          const py::list &>()(&Pysegy::cut),
           "cut a sub volume to a segy", py::arg("outname"), py::arg("startZ"),
           py::arg("endZ"), py::arg("startY"), py::arg("endY"),
           py::arg("custom_info") = py::list())
      .def("cut",
           overload_cast_<const std::string &, int, int, const py::list &>()(
               &Pysegy::cut),
           "cut a sub volume to a segy", py::arg("outname"), py::arg("startX"),
           py::arg("endX"), py::arg("custom_info") = py::list())
      .def("setSampleInterval", &Pysegy::setSampleInterval, py::arg("dt"))
      .def("setDataFormatCode", &Pysegy::setDataFormatCode, py::arg("format"))
      .def("setStartTime", &Pysegy::setStartTime, py::arg("start_time"))
      .def("setInlineInterval", &Pysegy::setInlineInterval, py::arg("dx"))
      .def("setCrosslineInterval", &Pysegy::setCrosslineInterval, py::arg("dy"))
      .def("setMinInline", &Pysegy::setMinInline, py::arg("minInline"))
      .def("setMinCrossline", &Pysegy::setMinCrossline, py::arg("minXline"))
      .def("textual_header", &Pysegy::textual_header, py::arg("coding") = 'u')
      .def("metaInfo", &Pysegy::metaInfo)
      //     .def("create", overload_cast_<const std::string
      //     &>()(&Pysegy::create),
      //          "create a segy", py::arg("segy_out_name"))
      //     .def("create",
      //          overload_cast_<const std::string &,
      //                         const pybind11::array_t<float>
      //                         &>()(&Pysegy::create),
      //          "create a segy from memory", py::arg("segy_out_name"),
      //          py::arg("src"))
      .def("create",
           static_cast<void (Pysegy::*)(const std::string &, const py::list &)>(
               &Pysegy::create),
           "create a segy", py::arg("segy_out_name"),
           py::arg("custom_info") = py::list())
      .def("create",
           static_cast<void (Pysegy::*)(const std::string &, const npfloat &,
                                        const py::list &)>(&Pysegy::create),
           "create a segy from memory", py::arg("segy_out_name"),
           py::arg("src"), py::arg("custom_info") = py::list())
      .def("set_size", &Pysegy::set_size, py::arg("sizeX"), py::arg("sizeY"),
           py::arg("sizeZ"))
      .def("close_file", &Pysegy::close_file)
      .def_property_readonly("is_crossline_fast_order",
                             &Pysegy::is_crossline_fast_order)
      .def_property_readonly("trace_count", &Pysegy::trace_count);

  m.def("fromfile_ignore_header", &fromfile_ignore_header,
        "read by ignoring header and specify shape", py::arg("segy_name"),
        py::arg("sizeZ"), py::arg("sizeY"), py::arg("sizeX"),
        py::arg("format") = 5);
  m.def("fromfile", &fromfile, "read from a file", py::arg("segy_name"),
        py::arg("iline") = 189, py::arg("xline") = 193, py::arg("istep") = 1,
        py::arg("xstep") = 1);
  m.def("tofile_ignore_header", &segy::tofile_ignore_header,
        "convert to binary file by ignoring header and specify shape",
        py::arg("segy_name"), py::arg("out_name"), py::arg("sizeX"),
        py::arg("sizeY"), py::arg("sizeZ"), py::arg("format") = 5);
  m.def("tofile", &segy::tofile, "convert to binary file", py::arg("segy_name"),
        py::arg("out_name"), py::arg("iline") = 189, py::arg("xline") = 193,
        py::arg("istep") = 1, py::arg("xstep") = 1);
  m.def("collect", &collect, "colloct all trace (as 2D data)",
        py::arg("segy_in"), py::arg("beg") = -1, py::arg("end") = 0);
  m.def("create_by_sharing_header",
        overload_cast_<const std::string &, const std::string &,
                       const npfloat &, int, int, int, int, const py::object &,
                       const py::list &>()(&create_by_sharing_header),
        "create a segy file using a existed segy header", py::arg("segy_name"),
        py::arg("header_segy"), py::arg("src"), py::arg("iline") = 189,
        py::arg("xline") = 193, py::arg("istep") = 1, py::arg("xstep") = 1,
        py::arg("offset") = py::none(), py::arg("custom_info") = py::list());
  m.def("create_by_sharing_header",
        overload_cast_<const std::string &, const std::string &,
                       const std::string &, const py::sequence &, int, int, int,
                       int, const py::object &, const py::list &>()(
            &create_by_sharing_header),
        "create a segy file using a existed segy header", py::arg("segy_name"),
        py::arg("header_segy"), py::arg("src_file"), py::arg("shape"),
        py::arg("iline") = 189, py::arg("xline") = 193, py::arg("istep") = 1,
        py::arg("xstep") = 1, py::arg("offset") = py::none(),
        py::arg("custom_info") = py::list());
  m.def("disable_progressbar", &segy::disable_progressbar,
        "disable progress bar");
  m.def("_load_prestack3D", &_load_prestack3D, "load prestack 3D SEG-Y",
        py::arg("segy_name"), py::arg("shape"), py::arg("min_iline"),
        py::arg("max_iline"), py::arg("min_xline"), py::arg("max_xline"),
        py::arg("min_offset"), py::arg("max_offset"), py::arg("istep"),
        py::arg("xstep"), py::arg("ostep"), py::arg("iline"), py::arg("xline"),
        py::arg("offset") = 37, py::arg("fill") = 0);

  m.attr("kBinaryHeaderHelp") =
      py::cast(&segy::kBinaryHeaderHelp, py::return_value_policy::reference);
  m.attr("kTraceHeaderHelp") =
      py::cast(&segy::kTraceHeaderHelp, py::return_value_policy::reference);

  py::class_<segy::MetaInfo>(m, "MetaInfo")
      .def(py::init<>())
      .def_readwrite("sizeX", &segy::MetaInfo::sizeX)
      .def_readwrite("sizeY", &segy::MetaInfo::sizeY)
      .def_readwrite("sizeZ", &segy::MetaInfo::sizeZ)
      .def_readwrite("trace_count", &segy::MetaInfo::trace_count)
      .def_readwrite("sample_interval", &segy::MetaInfo::sample_interval)
      .def_readwrite("data_format", &segy::MetaInfo::data_format)
      .def_readwrite("Y_interval", &segy::MetaInfo::Y_interval)
      .def_readwrite("Z_interval", &segy::MetaInfo::Z_interval)
      .def_readwrite("start_time", &segy::MetaInfo::start_time)
      .def_readwrite("scalar", &segy::MetaInfo::scalar)
      .def_readwrite("min_inline", &segy::MetaInfo::min_inline)
      .def_readwrite("max_inline", &segy::MetaInfo::max_inline)
      .def_readwrite("min_crossline", &segy::MetaInfo::min_crossline)
      .def_readwrite("max_crossline", &segy::MetaInfo::max_crossline)
      .def_readwrite("isNormalSegy", &segy::MetaInfo::isNormalSegy)
      .def_readwrite("fillNoValue", &segy::MetaInfo::fillNoValue)
      .def_readwrite("inline_field", &segy::MetaInfo::inline_field)
      .def_readwrite("crossline_field", &segy::MetaInfo::crossline_field)
      .def_readwrite("X_field", &segy::MetaInfo::X_field)
      .def_readwrite("Y_field", &segy::MetaInfo::Y_field)
      .def_readwrite("inline_step", &segy::MetaInfo::inline_step)
      .def_readwrite("crossline_step", &segy::MetaInfo::crossline_step)
      .def_readwrite("trace_sorting_code", &segy::MetaInfo::trace_sorting_code);
}
