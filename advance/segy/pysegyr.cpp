#include "segyc.hpp"
#include "segyrw.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using npfloat = py::array_t<float, py::array::c_style | py::array::forcecast>;
using npint32 = py::array_t<int32_t, py::array::c_style | py::array::forcecast>;

class SegyRWpy : public segy::SegyRW {
public:
  using segy::SegyRW::create_by_sharing_header;
  using segy::SegyRW::SegyRW;

  npfloat read4d(int is, int ie, int xs, int xe, int os, int oe, int ts,
                 int te) {
    // TODO: check bound
    int nt = te - ts;
    int ni = ie - is;
    int nx = xe - xs;
    int no = oe - os;
    auto data = py::array_t<float>({ni, nx, no, nt});
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::read4d(ptr, is, ie, xs, xe, os, oe, ts, te);
    return data;
  }

  npfloat read3d(int is, int ie, int xs, int xe, int ts, int te) {
    int nt = te - ts;
    int ni = ie - is;
    int nx = xe - xs;
    auto data = py::array_t<float>({ni, nx, nt});
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::read3d(ptr, is, ie, xs, xe, ts, te);
    return data;
  }

  npfloat read() {
    py::array_t<float> out(shape());
    auto buff = out.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::read(ptr);
    return out;
  }

  void create_by_sharing_header(const std::string &segy_name,
                                const npfloat &src, const py::list &ranges,
                                bool is2d = false,
                                const std::string &textual = "") {
    auto buff = src.request();
    float *ptr = static_cast<float *>(buff.ptr);
    auto rge = ranges.cast<std::vector<int>>();
    create_by_sharing_header(segy_name, ptr, rge, is2d, textual);
  }

  npfloat itrace(int n) {
    py::array_t<float> out(m_meta.nt);
    auto buff = out.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::itrace(ptr, n);
    return out;
  }

  npfloat collect(int beg, int end, int tbeg, int tend) {
    // TODO: check bound
    auto data = py::array_t<float>({end - beg, tend - tbeg});
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::collect(ptr, beg, end, tbeg, tend);
    return data;
  }

  npfloat collect(const npint32 &index, int tbeg, int tend) {
    auto buff = index.request();
    if (buff.ndim != 1) {
      throw std::runtime_error("Input index must be a 1D data.");
    }
    int N = index.shape()[0];
    int32_t *idx = static_cast<int32_t *>(buff.ptr);

    auto data = py::array_t<float>({N, tend - tbeg});
    auto buff2 = data.request();
    float *ptr = static_cast<float *>(buff2.ptr);
    SegyRW::collect(ptr, idx, N, tbeg, tend);
    return data;
  }

  npint32 get_trace_keys(const py::list &keys, const py::list &length, int beg,
                         int end) {
    auto keysvec = keys.cast<std::vector<int>>();
    auto lengthvec = length.cast<std::vector<int>>();
    int n1 = end - beg;
    int n2 = keysvec.size();
    py::array_t<int> out({n1, n2});
    auto buff = out.request();
    int *ptr = static_cast<int *>(buff.ptr);
    SegyRW::get_trace_keys(ptr, keysvec, lengthvec, beg, end);
    return out;
  }

  // for write
  void write_itrace(const npfloat &data, int n) {
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::write_itrace(ptr, n);
  }

  void write_traces(const npfloat &data, int beg, int end, int tbeg, int tend) {
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::write_traces(ptr, beg, end, tbeg, tend);
  }

  void write_traces(const npfloat &data, const npint32 &index, int tbeg,
                    int tend) {
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    auto buff2 = index.request();
    int32_t *idx = static_cast<int32_t *>(buff2.ptr);
    SegyRW::write_traces(ptr, idx, index.shape()[0], tbeg, tend);
  }

  void write(const npfloat &data) {
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::write(ptr);
  }

  void write3d(const npfloat &data, int is, int ie, int xs, int xe, int ts,
               int te) {
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::write3d(ptr, is, ie, xs, xe, ts, te);
  }

  void write4d(const npfloat &data, int is, int ie, int xs, int xe, int os,
               int oe, int ts, int te) {
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyRW::write4d(ptr, is, ie, xs, xe, os, oe, ts, te);
  }
};

class segyc : public segy::SegyC {
public:
  using segy::SegyC::SegyC;

  void add_trace(const npfloat &src, int idx, int iline, int xline,
                 int offset = 0) {
    // TODO: check src' length == nt
    auto buff = src.request();
    float *ptr = static_cast<float *>(buff.ptr);
    segy::SegyC::add_trace(ptr, idx, iline, xline, offset);
  }

  void create(const npfloat &src, const npint32 &xyico) {
    // TODO: check bound
    auto buff = src.request();
    float *ptr = static_cast<float *>(buff.ptr);
    auto buff2 = xyico.request();
    int32_t *ptr2 = static_cast<int32_t *>(buff2.ptr);
    segy::SegyC::create(ptr, ptr2);
  }
};

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(_CXX_SEGYR, m) {
  py::class_<SegyRWpy>(m, "SegyRWpy")
      .def(py::init<std::string>())

      // location
      .def("setLocations", &SegyRWpy::setLocations, py::arg("iline"),
           py::arg("xline"), py::arg("offset") = 37)
      .def("setInlineLocation", &SegyRWpy::setInlineLocation, py::arg("iline"))
      .def("setCrosslineLocation", &SegyRWpy::setCrosslineLocation,
           py::arg("xline"))
      .def("setOffsetLocation", &SegyRWpy::setOffsetLocation, py::arg("offset"))
      .def("setXLocation", &SegyRWpy::setXLocation, py::arg("xloc"))
      .def("setYLocation", &SegyRWpy::setYLocation, py::arg("yloc"))
      .def("setXYLocation", &SegyRWpy::setXYLocations, py::arg("xloc"),
           py::arg("yloc"))
      .def("setInlineStep", &SegyRWpy::setInlineStep, py::arg("istep"))
      .def("setCrosslineStep", &SegyRWpy::setCrosslineStep, py::arg("xstep"))
      .def("setOffsetStep", &SegyRWpy::setOffsetStep, py::arg("ostep"))
      .def("setSteps", &SegyRWpy::setSteps, py::arg("istep"), py::arg("xstep"),
           py::arg("ostep") = 1)
      .def("setFill", &SegyRWpy::setFill, py::arg("fill"))

      // read func base
      .def("textual_header", &SegyRWpy::textual_header, py::arg("coding") = 'u')
      .def("get_trace_keys", &SegyRWpy::get_trace_keys, py::arg("leys"),
           py::arg("length"), py::arg("beg"), py::arg("end"))
      .def("itrace", &SegyRWpy::itrace, py::arg("n"))
      .def("collect",
           overload_cast_<const npint32 &, int, int>()(&SegyRWpy::collect),
           py::arg("index"), py::arg("tbeg"), py::arg("tend"))
      .def("collect", overload_cast_<int, int, int, int>()(&SegyRWpy::collect),
           py::arg("beg"), py::arg("end"), py::arg("tbeg"), py::arg("tend"))

      // read
      .def("set_segy_type", &SegyRWpy::set_segy_type, py::arg("ndim"))
      .def("scan", &SegyRWpy::scan)
      .def("read4d", &SegyRWpy::read4d, py::arg("ib"), py::arg("ie"),
           py::arg("xb"), py::arg("xe"), py::arg("ob"), py::arg("oe"),
           py::arg("tb"), py::arg("te"))
      .def("read3d", &SegyRWpy::read3d, py::arg("ib"), py::arg("ie"),
           py::arg("xb"), py::arg("xe"), py::arg("tb"), py::arg("te"))
      .def("read", &SegyRWpy::read)
      .def("tofile", &SegyRWpy::tofile, py::arg("binary_out_name"),
           py::arg("as_2d") = false)
      .def("cut", &SegyRWpy::cut, py::arg("outname"), py::arg("ranges"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           overload_cast_<const std::string &, const npfloat &,
                          const py::list &, bool, const std::string &>()(
               &SegyRWpy::create_by_sharing_header),
           py::arg("segy_name"), py::arg("src"), py::arg("ranges"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           overload_cast_<const std::string &, const std::string &,
                          const std::vector<int> &, const std::vector<int> &,
                          bool, const std::string &>()(
               &SegyRWpy::create_by_sharing_header))

      // for write
      .def("set_bkeyi2", &SegyRWpy::set_bkeyi2, py::arg("loc"), py::arg("val"))
      .def("set_bkeyi4", &SegyRWpy::set_bkeyi4, py::arg("loc"), py::arg("val"))
      //  .def("set_bkeyi8", &SegyRWpy::set_bkeyi8, py::arg("loc"),
      //  py::arg("val"))
      .def("set_keyi2", &SegyRWpy::set_keyi2, py::arg("n"), py::arg("loc"),
           py::arg("val"))
      .def("set_keyi4", &SegyRWpy::set_keyi4, py::arg("n"), py::arg("loc"),
           py::arg("val"))
      //  .def("set_keyi8", &SegyRWpy::set_keyi8, py::arg("n"), py::arg("loc"),
      //       py::arg("val"))
      .def("set_iline", &SegyRWpy::set_iline, py::arg("n"), py::arg("val"))
      .def("set_xline", &SegyRWpy::set_xline, py::arg("n"), py::arg("val"))
      .def("set_offset", &SegyRWpy::set_offset, py::arg("n"), py::arg("val"))
      .def("set_coordx", &SegyRWpy::set_coordx, py::arg("n"), py::arg("val"))
      .def("set_coordy", &SegyRWpy::set_coordy, py::arg("n"), py::arg("val"))

      .def("write_itrace", &SegyRWpy::write_itrace, py::arg("data"),
           py::arg("n"))
      .def("write_traces",
           overload_cast_<const npfloat &, int, int, int, int>()(
               &SegyRWpy::write_traces),
           py::arg("data"), py::arg("beg"), py::arg("end"), py::arg("tbeg"),
           py::arg("tend"))
      .def("write_traces",
           overload_cast_<const npfloat &, const npint32 &, int, int>()(
               &SegyRWpy::write_traces),
           py::arg("data"), py::arg("index"), py::arg("tbeg"), py::arg("tend"))
      .def("write", &SegyRWpy::write, py::arg("data"))
      .def("write3d", &SegyRWpy::write3d, py::arg("data"), py::arg("is"),
           py::arg("ie"), py::arg("xs"), py::arg("xe"), py::arg("ts"),
           py::arg("te"))
      .def("write4d", &SegyRWpy::write4d, py::arg("data"), py::arg("is"),
           py::arg("ie"), py::arg("xs"), py::arg("xe"), py::arg("os"),
           py::arg("oe"), py::arg("ts"), py::arg("te"))

      ;

  py::class_<segyc>(m, "segyc")
      .def(py::init<std::string, int, int>())
      .def(py::init<std::string, int, int, int>())
      .def(py::init<std::string, int, int, int, int>())

      // location
      .def("setLocations", &segyc::setLocations, py::arg("iline"),
           py::arg("xline"), py::arg("offset") = 37)
      .def("setInlineLocation", &segyc::setInlineLocation, py::arg("iline"))
      .def("setCrosslineLocation", &segyc::setCrosslineLocation,
           py::arg("xline"))
      .def("setOffsetLocation", &segyc::setOffsetLocation, py::arg("offset"))
      .def("setXLocation", &segyc::setXLocation, py::arg("xloc"))
      .def("setYLocation", &segyc::setYLocation, py::arg("yloc"))
      .def("setXYLocation", &segyc::setXYLocations, py::arg("xloc"),
           py::arg("yloc"))
      .def("setInlineStep", &segyc::setInlineStep, py::arg("istep"))
      .def("setCrosslineStep", &segyc::setCrosslineStep, py::arg("xstep"))
      .def("setOffsetStep", &segyc::setOffsetStep, py::arg("ostep"))
      .def("setSteps", &segyc::setSteps, py::arg("istep"), py::arg("xstep"),
           py::arg("ostep") = 1)

      .def("setStartTime", &segyc::setStartTime, py::arg("start_time"))
      .def("setInlineInterval", &segyc::setInlineInterval, py::arg("di"))
      .def("setCrosslineInterval", &segyc::setCrosslineInterval, py::arg("dx"))
      .def("setStartInline", &segyc::setStartInline, py::arg("il"))
      .def("setStartCrossline", &segyc::setStartCrossline, py::arg("xl"))
      .def("setStartOffset", &segyc::setStartOffset, py::arg("of"))

      .def("copy_textual_from", &segyc::copy_textual_from, py::arg("segyname"))
      .def("copy_bheader_from", &segyc::copy_bheader_from, py::arg("segyname"))
      .def("copy_theader_from", &segyc::copy_theader_from, py::arg("segyname"))

      .def("set_bkeyi2", &segyc::set_bkeyi2, py::arg("loc"), py::arg("val"))
      .def("set_bkeyi4", &segyc::set_bkeyi4, py::arg("loc"), py::arg("val"))
      .def("set_bkeyi8", &segyc::set_bkeyi8, py::arg("loc"), py::arg("val"))
      .def("set_keyi2", &segyc::set_keyi2, py::arg("n"), py::arg("loc"),
           py::arg("val"))
      .def("set_keyi4", &segyc::set_keyi4, py::arg("n"), py::arg("loc"),
           py::arg("val"))
      .def("set_keyi8", &segyc::set_keyi8, py::arg("n"), py::arg("loc"),
           py::arg("val"))
      .def("set_iline", &segyc::set_iline, py::arg("n"), py::arg("val"))
      .def("set_xline", &segyc::set_xline, py::arg("n"), py::arg("val"))
      .def("set_offset", &segyc::set_offset, py::arg("n"), py::arg("val"))
      .def("set_coordx", &segyc::set_coordx, py::arg("n"), py::arg("val"))
      .def("set_coordy", &segyc::set_coordy, py::arg("n"), py::arg("val"))
      .def("add_trace", &segyc::add_trace, py::arg("src"), py::arg("idx"),
           py::arg("iline"), py::arg("xline"), py::arg("offset") = 0)
      .def("create", &segyc::create, py::arg("src"), py::arg("xyico"));
  ;
}
