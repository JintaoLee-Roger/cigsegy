/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

// #include "segyc.hpp"
#include "segyrw.h"
#include "sutils.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using npfloat = py::array_t<float, py::array::c_style | py::array::forcecast>;
using npint32 = py::array_t<int32_t, py::array::c_style | py::array::forcecast>;
using npuchar = py::array_t<uchar, py::array::c_style | py::array::forcecast>;

namespace segy {

class Pysegy : public SegyRW {
public:
  using segy::SegyRW::create_by_sharing_header;
  using SegyRW::SegyRW;

  npfloat read4d(int is, int ie, int xs, int xe, int os, int oe, int ts,
                 int te) {
    // TODO: check bound
    int nt = te - ts;
    int ni = ie - is;
    int nx = xe - xs;
    int no = oe - os;
    auto data = py::array_t<float>({ni, nx, no, nt});
    float *ptr = data.mutable_data();
    SegyRW::read4d(ptr, is, ie, xs, xe, os, oe, ts, te);
    return data;
  }

  npfloat read3d(int is, int ie, int xs, int xe, int ts, int te) {
    int nt = te - ts;
    int ni = ie - is;
    int nx = xe - xs;
    auto data = py::array_t<float>({ni, nx, nt});
    float *ptr = data.mutable_data();
    SegyRW::read3d(ptr, is, ie, xs, xe, ts, te);
    return data;
  }

  npfloat read() {
    py::array_t<float> out(shape());
    float *ptr = out.mutable_data();
    SegyRW::read(ptr);
    return out;
  }

  void create_by_sharing_header(const std::string &segy_name,
                                const npfloat &src, const py::list &ranges,
                                bool is2d = false,
                                const std::string &textual = "") {
    const float *ptr = src.data();
    auto rge = ranges.cast<std::vector<int>>();
    create_by_sharing_header(segy_name, ptr, rge, is2d, textual);
  }

  npfloat itrace(int n) {
    py::array_t<float> out(m_meta.nt);
    float *ptr = out.mutable_data();
    SegyRW::itrace(ptr, n);
    return out;
  }

  npfloat collect(int beg, int end, int tbeg, int tend) {
    // TODO: check bound
    auto data = py::array_t<float>({end - beg, tend - tbeg});
    float *ptr = data.mutable_data();
    SegyRW::collect(ptr, beg, end, tbeg, tend);
    return data;
  }

  npfloat collect(const npint32 &index, int tbeg, int tend) {
    if (index.ndim() != 1) {
      throw std::runtime_error("Input index must be a 1D data.");
    }
    int N = index.shape()[0];
    const int32_t *idx = index.data();

    auto data = py::array_t<float>({N, tend - tbeg});
    float *ptr = data.mutable_data();
    SegyRW::collect(ptr, idx, N, tbeg, tend);
    return data;
  }

  npuchar get_binary_header() {
    py::array_t<uchar> out(segy::kBinaryHeaderSize);
    uchar *ptr = out.mutable_data();
    SegyRW::get_binary_header(ptr);
    return out;
  }
  npuchar get_trace_header(int n) {
    py::array_t<uchar> out(segy::kTraceHeaderSize);
    uchar *ptr = out.mutable_data();
    SegyRW::get_trace_header(ptr, n);
    return out;
  }

  npint32 get_trace_keys(const py::list &keys, const py::list &length, int beg,
                         int end) {
    auto keysvec = keys.cast<std::vector<int>>();
    auto lengthvec = length.cast<std::vector<int>>();
    int n1 = end - beg;
    int n2 = keysvec.size();
    py::array_t<int> out({n1, n2});
    int *ptr = out.mutable_data();
    SegyRW::get_trace_keys(ptr, keysvec, lengthvec, beg, end);
    return out;
  }

  //  for write
  void write_itrace(const npfloat &data, int n) {
    const float *ptr = data.data();
    SegyRW::write_itrace(ptr, n);
  }

  void write_traces(const npfloat &data, int beg, int end, int tbeg, int tend) {
    const float *ptr = data.data();
    SegyRW::write_traces(ptr, beg, end, tbeg, tend);
  }

  void write_traces(const npfloat &data, const npint32 &index, int tbeg,
                    int tend) {
    const float *ptr = data.data();
    const int32_t *idx = index.data();
    SegyRW::write_traces(ptr, idx, index.shape()[0], tbeg, tend);
  }

  void write(const npfloat &data) {
    const float *ptr = data.data();
    SegyRW::write(ptr);
  }

  void write3d(const npfloat &data, int is, int ie, int xs, int xe, int ts,
               int te) {
    const float *ptr = data.data();
    SegyRW::write3d(ptr, is, ie, xs, xe, ts, te);
  }

  void write4d(const npfloat &data, int is, int ie, int xs, int xe, int os,
               int oe, int ts, int te) {
    const float *ptr = data.data();
    SegyRW::write4d(ptr, is, ie, xs, xe, os, oe, ts, te);
  }

  std::map<std::string, int> get_keylocs() const {
    return {{"iline", m_keys.iline},   {"xline", m_keys.xline},
            {"offset", m_keys.offset}, {"xloc", m_keys.xloc},
            {"yloc", m_keys.yloc},     {"istep", m_keys.istep},
            {"xstep", m_keys.xstep},   {"ostep", m_keys.ostep}};
  }

  std::map<std::string, pybind11::object> get_metainfo() const {
    return {{"nt", pybind11::int_(m_meta.nt)},
            {"no", pybind11::int_(m_meta.no)},
            {"nx", pybind11::int_(m_meta.nx)},
            {"ni", pybind11::int_(m_meta.ni)},
            {"ntrace", pybind11::int_(m_meta.ntrace)},
            {"tracesize", pybind11::int_(m_meta.tracesize)},
            {"dt", pybind11::int_(m_meta.dt)},
            {"dx", pybind11::float_(m_meta.dx)},
            {"di", pybind11::float_(m_meta.di)},
            {"scalar", pybind11::int_(m_meta.scalar)},
            {"dformat", pybind11::int_(m_meta.dformat)},
            {"start_time", pybind11::int_(m_meta.start_time)},
            {"start_iline", pybind11::int_(m_meta.start_iline)},
            {"end_iline", pybind11::int_(m_meta.end_iline)},
            {"start_xline", pybind11::int_(m_meta.start_xline)},
            {"end_xline", pybind11::int_(m_meta.end_xline)},
            {"start_offset", pybind11::int_(m_meta.start_offset)},
            {"end_offset", pybind11::int_(m_meta.end_offset)},
            {"trace_sorting_code", pybind11::int_(m_meta.trace_sorting_code)},
            {"esize", pybind11::int_(m_meta.esize)},
            {"fillNoValue", pybind11::float_(m_meta.fillNoValue)},
            {"ndim", pybind11::int_(ndim())}};
  }

  npint32 get_lineInfo() {
    if (ndim() == 2) {
      throw std::runtime_error("get_lineInfo function valid when ndim > 2");
    }
    py::array_t<int> out;
    if (ndim() == 3) {
      out = py::array_t<int>({m_meta.ni, 3});
    } else {
      out = py::array_t<int>({m_meta.ni, m_meta.nx, 4});
    }
    int *ptr = out.mutable_data();

    if (ndim() == 3) {
      for (auto linfo : m_iinfos) {
        ptr[0] = linfo.line;
        if (linfo.count == 0) {
          ptr[1] = -1;
          ptr[2] = -1;
        } else {
          ptr[1] = linfo.itstart;
          ptr[2] = linfo.itend;
        }
        ptr += 3;
      }
    } else {
      for (auto linfo : m_iinfos) {
        // TODO: fill
        for (auto xinfo : linfo.xinfos) {
          ptr[0] = linfo.line;
          ptr[1] = xinfo.line;
          if (xinfo.count == 0) {
            ptr[2] = -1;
            ptr[3] = -1;
          } else {
            ptr[2] = xinfo.itstart;
            ptr[3] = xinfo.itend;
          }
          ptr += 4;
        }
      }
    }

    return out;
  }
};

// class SegyCpy : public segy::SegyC {
// public:
//   using segy::SegyC::SegyC;

//   void add_trace(const npfloat &src, int idx, int iline, int xline,
//                  int offset = 0) {
//     // TODO: check src' length == nt
//     auto buff = src.request();
//     float *ptr = static_cast<float *>(buff.ptr);
//     segy::SegyC::add_trace(ptr, idx, iline, xline, offset);
//   }

//   void create(const npfloat &src, const npint32 &xyico) {
//     // TODO: check bound
//     auto buff = src.request();
//     float *ptr = static_cast<float *>(buff.ptr);
//     auto buff2 = xyico.request();
//     int32_t *ptr2 = static_cast<int32_t *>(buff2.ptr);
//     segy::SegyC::create(ptr, ptr2);
//   }
// };

npfloat ieees_to_ibms(const npfloat &ieee_arr, bool is_little_endian) {
  auto shape = ieee_arr.shape();
  std::vector<ssize_t> shape_vec(shape, shape + ieee_arr.ndim());
  auto size = ieee_arr.size();

  py::array_t<float> ibm_arr(shape_vec);

  const float *ieee_ptr = ieee_arr.data();
  float *ibm_ptr = ibm_arr.mutable_data();

  for (size_t i = 0; i < size; ++i) {
    ibm_ptr[i] = segy::ieee_to_ibm(ieee_ptr[i], is_little_endian);
  }

  return ibm_arr;
}

npfloat ibms_to_ieees(const npfloat &ibm_arr, bool is_big_endian) {
  auto shape = ibm_arr.shape();
  std::vector<ssize_t> shape_vec(shape, shape + ibm_arr.ndim());
  auto size = ibm_arr.size();

  py::array_t<float> ieee_arr(shape_vec);

  const float *ibm_ptr = ibm_arr.data();
  float *ieee_ptr = ieee_arr.mutable_data();

  for (size_t i = 0; i < size; ++i) {
    ieee_ptr[i] = segy::ibm_to_ieee(ibm_ptr[i], is_big_endian);
  }

  return ibm_arr;
}

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(_CXX_SEGY, m) {
  py::class_<Pysegy>(m, "Pysegy")
      .def(py::init<const std::string &>())
      .def("close", &Pysegy::close_file)

      // location
      .def("setLocations", &Pysegy::setLocations, py::arg("iline"),
           py::arg("xline"), py::arg("offset") = 37)
      .def("setInlineLocation", &Pysegy::setInlineLocation, py::arg("iline"))
      .def("setCrosslineLocation", &Pysegy::setCrosslineLocation,
           py::arg("xline"))
      .def("setOffsetLocation", &Pysegy::setOffsetLocation, py::arg("offset"))
      .def("setXLocation", &Pysegy::setXLocation, py::arg("xloc"))
      .def("setYLocation", &Pysegy::setYLocation, py::arg("yloc"))
      .def("setXYLocations", &Pysegy::setXYLocations, py::arg("xloc"),
           py::arg("yloc"))
      .def("setInlineStep", &Pysegy::setInlineStep, py::arg("istep"))
      .def("setCrosslineStep", &Pysegy::setCrosslineStep, py::arg("xstep"))
      .def("setOffsetStep", &Pysegy::setOffsetStep, py::arg("ostep"))
      .def("setSteps", &Pysegy::setSteps, py::arg("istep"), py::arg("xstep"),
           py::arg("ostep") = 1)
      .def("setFill", &Pysegy::setFill, py::arg("fill"))

      // read func base
      .def("textual_header", &Pysegy::textual_header, py::arg("coding") = 'u')
      .def("bkeyi2", &Pysegy::bkeyi2, py::arg("loc"))
      .def("bkeyi4", &Pysegy::bkeyi4, py::arg("loc"))
      .def("keyi2", &Pysegy::keyi2, py::arg("n"), py::arg("loc"))
      .def("keyi4", &Pysegy::keyi4, py::arg("n"), py::arg("loc"))
      .def("iline", &Pysegy::iline, py::arg("n"))
      .def("xline", &Pysegy::xline, py::arg("n"))
      .def("offset", &Pysegy::offset, py::arg("n"))
      .def("coordx", &Pysegy::coordx, py::arg("n"))
      .def("coordy", &Pysegy::coordy, py::arg("n"))
      .def("get_binary_header", &Pysegy::get_binary_header)
      .def("get_trace_header", &Pysegy::get_trace_header, py::arg("n"))
      .def("get_trace_keys", &Pysegy::get_trace_keys, py::arg("leys"),
           py::arg("length"), py::arg("beg"), py::arg("end"))
      .def("itrace", &Pysegy::itrace, py::arg("n"))
      .def("collect",
           overload_cast_<const npint32 &, int, int>()(&Pysegy::collect),
           py::arg("index"), py::arg("tbeg"), py::arg("tend"))
      .def("collect", overload_cast_<int, int, int, int>()(&Pysegy::collect),
           py::arg("beg"), py::arg("end"), py::arg("tbeg"), py::arg("tend"))

      // meta info
      .def("get_keylocs", &Pysegy::get_keylocs)
      .def("get_metainfo", &Pysegy::get_metainfo)
      .def("get_lineInfo", &Pysegy::get_lineInfo)

      // read
      .def("set_segy_type", &Pysegy::set_segy_type, py::arg("ndim"))
      .def("scan", &Pysegy::scan)
      .def("read4d", &Pysegy::read4d, py::arg("ib"), py::arg("ie"),
           py::arg("xb"), py::arg("xe"), py::arg("ob"), py::arg("oe"),
           py::arg("tb"), py::arg("te"))
      .def("read3d", &Pysegy::read3d, py::arg("ib"), py::arg("ie"),
           py::arg("xb"), py::arg("xe"), py::arg("tb"), py::arg("te"))
      .def("read", &Pysegy::read)
      .def("tofile", &Pysegy::tofile, py::arg("binary_out_name"),
           py::arg("as_2d") = false)
      .def("cut", &Pysegy::cut, py::arg("outname"), py::arg("ranges"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           overload_cast_<const std::string &, const npfloat &,
                          const py::list &, bool, const std::string &>()(
               &Pysegy::create_by_sharing_header),
           py::arg("segy_name"), py::arg("src"), py::arg("ranges"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           overload_cast_<const std::string &, const std::string &,
                          const std::vector<int> &, const std::vector<int> &,
                          bool, const std::string &>()(
               &Pysegy::create_by_sharing_header))

      // for write
      .def("set_bkeyi2", &Pysegy::set_bkeyi2, py::arg("loc"), py::arg("val"))
      .def("set_bkeyi4", &Pysegy::set_bkeyi4, py::arg("loc"), py::arg("val"))
      .def("set_keyi2", &Pysegy::set_keyi2, py::arg("n"), py::arg("loc"),
           py::arg("val"))
      .def("set_keyi4", &Pysegy::set_keyi4, py::arg("n"), py::arg("loc"),
           py::arg("val"))
      .def("set_iline", &Pysegy::set_iline, py::arg("n"), py::arg("val"))
      .def("set_xline", &Pysegy::set_xline, py::arg("n"), py::arg("val"))
      .def("set_offset", &Pysegy::set_offset, py::arg("n"), py::arg("val"))
      .def("set_coordx", &Pysegy::set_coordx, py::arg("n"), py::arg("val"))
      .def("set_coordy", &Pysegy::set_coordy, py::arg("n"), py::arg("val"))

      .def("write_itrace", &Pysegy::write_itrace, py::arg("data"), py::arg("n"))
      .def("write_traces",
           overload_cast_<const npfloat &, int, int, int, int>()(
               &Pysegy::write_traces),
           py::arg("data"), py::arg("beg"), py::arg("end"), py::arg("tbeg"),
           py::arg("tend"))
      .def("write_traces",
           overload_cast_<const npfloat &, const npint32 &, int, int>()(
               &Pysegy::write_traces),
           py::arg("data"), py::arg("index"), py::arg("tbeg"), py::arg("tend"))
      // .def("write", &Pysegy::write, py::arg("data"))
      // .def("write3d", &Pysegy::write3d, py::arg("data"), py::arg("ib"),
      //      py::arg("ie"), py::arg("xb"), py::arg("xe"), py::arg("tb"),
      //      py::arg("te"))
      // .def("write4d", &Pysegy::write4d, py::arg("data"), py::arg("ib"),
      //      py::arg("ie"), py::arg("xb"), py::arg("xe"), py::arg("ob"),
      //      py::arg("oe"), py::arg("tb"), py::arg("te"))

      .def_property_readonly("trace_count", &Pysegy::ntrace)
      .def_property_readonly("nt", &Pysegy::nt)
      .def_property_readonly("ndim", &Pysegy::ndim)

      ;

  //   py::class_<SegyCpy>(m, "SegyCpy")
  //       .def(py::init<std::string, int, int>())
  //       .def(py::init<std::string, int, int, int>())
  //       .def(py::init<std::string, int, int, int, int>())

  //       // location
  //       .def("setLocations", &SegyCpy::setLocations, py::arg("iline"),
  //            py::arg("xline"), py::arg("offset") = 37)
  //       .def("setInlineLocation", &SegyCpy::setInlineLocation,
  //       py::arg("iline")) .def("setCrosslineLocation",
  //       &SegyCpy::setCrosslineLocation,
  //            py::arg("xline"))
  //       .def("setOffsetLocation", &SegyCpy::setOffsetLocation,
  //       py::arg("offset")) .def("setXLocation", &SegyCpy::setXLocation,
  //       py::arg("xloc")) .def("setYLocation", &SegyCpy::setYLocation,
  //       py::arg("yloc")) .def("setXYLocation", &SegyCpy::setXYLocations,
  //       py::arg("xloc"),
  //            py::arg("yloc"))
  //       .def("setInlineStep", &SegyCpy::setInlineStep, py::arg("istep"))
  //       .def("setCrosslineStep", &SegyCpy::setCrosslineStep,
  //       py::arg("xstep")) .def("setOffsetStep", &SegyCpy::setOffsetStep,
  //       py::arg("ostep")) .def("setSteps", &SegyCpy::setSteps,
  //       py::arg("istep"), py::arg("xstep"),
  //            py::arg("ostep") = 1)

  //       .def("setStartTime", &SegyCpy::setStartTime, py::arg("start_time"))
  //       .def("setInlineInterval", &SegyCpy::setInlineInterval, py::arg("di"))
  //       .def("setCrosslineInterval", &SegyCpy::setCrosslineInterval,
  //            py::arg("dx"))
  //       .def("setStartInline", &SegyCpy::setStartInline, py::arg("il"))
  //       .def("setStartCrossline", &SegyCpy::setStartCrossline, py::arg("xl"))
  //       .def("setStartOffset", &SegyCpy::setStartOffset, py::arg("of"))

  //       .def("copy_textual_from", &SegyCpy::copy_textual_from,
  //            py::arg("segyname"))
  //       .def("copy_bheader_from", &SegyCpy::copy_bheader_from,
  //            py::arg("segyname"))
  //       .def("copy_theader_from", &SegyCpy::copy_theader_from,
  //            py::arg("segyname"))

  //       .def("set_bkeyi2", &SegyCpy::set_bkeyi2, py::arg("loc"),
  //       py::arg("val")) .def("set_bkeyi4", &SegyCpy::set_bkeyi4,
  //       py::arg("loc"), py::arg("val")) .def("set_keyi2",
  //       &SegyCpy::set_keyi2, py::arg("n"), py::arg("loc"),
  //            py::arg("val"))
  //       .def("set_keyi4", &SegyCpy::set_keyi4, py::arg("n"), py::arg("loc"),
  //            py::arg("val"))
  //       .def("set_iline", &SegyCpy::set_iline, py::arg("n"), py::arg("val"))
  //       .def("set_xline", &SegyCpy::set_xline, py::arg("n"), py::arg("val"))
  //       .def("set_offset", &SegyCpy::set_offset, py::arg("n"),
  //       py::arg("val")) .def("set_coordx", &SegyCpy::set_coordx,
  //       py::arg("n"), py::arg("val")) .def("set_coordy",
  //       &SegyCpy::set_coordy, py::arg("n"), py::arg("val")) .def("add_trace",
  //       &SegyCpy::add_trace, py::arg("src"), py::arg("idx"),
  //            py::arg("iline"), py::arg("xline"), py::arg("offset") = 0)
  //       .def("create", &SegyCpy::create, py::arg("src"), py::arg("xyico"));
  ;

  m.def("ieee_to_ibm", &segy::ieee_to_ibm, py::arg("value"),
        py::arg("is_little_endian"));
  m.def("ibm_to_ieee", &segy::ibm_to_ieee, py::arg("value"),
        py::arg("is_big_endian"));
  m.def("ieees_to_ibms", &ieees_to_ibms, py::arg("ieee_arr"),
        py::arg("is_little_endian"));
  m.def("ibms_to_ieees", &ibms_to_ieees, py::arg("ibm_arr"),
        py::arg("is_big_endian"));
}
} // namespace segy
