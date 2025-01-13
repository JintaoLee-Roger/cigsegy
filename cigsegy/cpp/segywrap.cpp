/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

// #include "segyc.hpp"
#include "pybind11/cast.h"
#include "pybind11/detail/common.h"
#include "segyrw.h"
#include "utils.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using npfloat = py::array_t<float, py::array::c_style | py::array::forcecast>;
using npint32 = py::array_t<int32_t, py::array::c_style | py::array::forcecast>;
using npuchar = py::array_t<uchar, py::array::c_style | py::array::forcecast>;

namespace segy {

void checkSignals() {
    if (PyErr_CheckSignals() != 0) {
        throw py::error_already_set();
    }
}


class Pysegy : public SegyRW {
public:
  using SegyRW::SegyRW;

  npfloat read4d(size_t ib, size_t ie, size_t xb, size_t xe, size_t ob,
                 size_t oe, size_t tb, size_t te) {
    if (ndim() != 4) {
      throw std::runtime_error("read4d function valid when ndim == 4");
    }
    if (ib > ie || ie > m_meta.ni || xb > xe || xe > m_meta.nx || ob > oe ||
        oe > m_meta.no || tb > te || te > m_meta.nt) {
      throw std::out_of_range("Index out of bound.");
    }

    size_t nt = te - tb;
    size_t ni = ie - ib;
    size_t nx = xe - xb;
    size_t no = oe - ob;
    auto data = py::array_t<float>({ni, nx, no, nt});
    float *ptr = data.mutable_data();
    SegyRW::read4d(ptr, ib, ie, xb, xe, ob, oe, tb, te);
    return data;
  }

  npfloat read3d(size_t ib, size_t ie, size_t xb, size_t xe, size_t tb,
                 size_t te) {
    // check bound
    if (ndim() != 3) {
      throw std::runtime_error("read3d function valid when ndim == 3");
    }
    if (ib > ie || ie > m_meta.ni || xb > xe || xe > m_meta.nx || tb > te ||
        te > m_meta.nt) {
      throw std::out_of_range("Index out of bound.");
    }

    size_t nt = te - tb;
    size_t ni = ie - ib;
    size_t nx = xe - xb;
    auto data = py::array_t<float>({ni, nx, nt});
    float *ptr = data.mutable_data();
    SegyRW::read3d(ptr, ib, ie, xb, xe, tb, te);
    return data;
  }

  npfloat read() {
    py::array_t<float> out(shape());
    float *ptr = out.mutable_data();
    SegyRW::read(ptr);
    return out;
  }

  npfloat read_tslice(size_t t, size_t stepi = 1, size_t stepx = 1) {
    if (stepi < 1 || stepx < 1) {
      throw std::runtime_error("stepi and stepx must be greater than 0.");
    }
    size_t ni = (m_meta.ni + stepi - 1) / stepi;
    size_t nx = (m_meta.nx + stepx - 1) / stepx;
    auto data = py::array_t<float>({ni, nx});
    float *ptr = data.mutable_data();
    SegyRW::read_tslice(ptr, t, stepi, stepx);
    return data;
  }

  void create_by_sharing_header(const std::string &segy_name,
                                const npfloat &src, const py::list &start,
                                bool is2d = false,
                                const std::string &textual = "") {
    const float *ptr = src.data();
    std::vector<size_t> shape(src.shape(), src.shape() + src.ndim());
    auto sta = start.cast<std::vector<size_t>>();
    SegyRW::create_by_sharing_header(segy_name, ptr, shape, sta, is2d, textual);
  }

  void create_by_sharing_header(const std::string &segy_name,
                                const std::string &src_file,
                                const py::list &shape, const py::list &start,
                                bool is2d = false,
                                const std::string &textual = "") {
    auto sha = shape.cast<std::vector<size_t>>();
    auto sta = start.cast<std::vector<size_t>>();
    SegyRW::create_by_sharing_header(segy_name, src_file, sha, sta, is2d,
                                     textual);
  }

  npfloat itrace(size_t n) {
    if (n >= m_meta.ntrace) {
      throw std::out_of_range("Index out of bound. Index: " +
                              std::to_string(n));
    }

    py::array_t<float> out(m_meta.nt);
    float *ptr = out.mutable_data();
    SegyRW::itrace(ptr, n);
    return out;
  }

  npfloat collect(size_t beg, size_t end, size_t tbeg, size_t tend) {
    if (beg > end || end > m_meta.ntrace || tbeg > tend || tend > m_meta.nt) {
      throw std::out_of_range("Index out of bound.");
    }

    auto data = py::array_t<float>({end - beg, tend - tbeg});
    float *ptr = data.mutable_data();
    SegyRW::collect(ptr, beg, end, tbeg, tend);
    return data;
  }

  npfloat collect(const npint32 &index, size_t tbeg, size_t tend) {
    if (index.ndim() != 1 && index.size() == 0) {
      throw std::runtime_error("Input index must be a 1D data.");
    }
    if (tbeg > tend || tend > m_meta.nt) {
      throw std::out_of_range("`tbeg` or `tend` index out of bound.");
    }

    size_t N = index.shape()[0];
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

  npuchar get_trace_header(size_t n) {
    if (n > ntrace()) {
      throw std::out_of_range("Index out of bound." + std::to_string(n));
    }

    py::array_t<uchar> out(segy::kTraceHeaderSize);
    uchar *ptr = out.mutable_data();
    SegyRW::get_trace_header(ptr, n);
    return out;
  }

  npint32 get_trace_keys(const py::list &keys, const py::list &length,
                         size_t beg, size_t end) {
    if (beg > end || end > ntrace()) {
      throw std::out_of_range("`beg` or `end` Index out of bound.");
    }
    if (keys.size() != length.size()) {
      throw std::runtime_error("`keys` and `length` must have the same size.");
    }

    auto keysvec = keys.cast<std::vector<size_t>>();
    auto lengthvec = length.cast<std::vector<size_t>>();

    size_t n1 = end - beg;
    size_t n2 = keysvec.size();
    py::array_t<int> out({n1, n2});
    int *ptr = out.mutable_data();
    SegyRW::get_trace_keys(ptr, keysvec, lengthvec, beg, end);
    return out;
  }

  //  for write
  void write_itrace(const npfloat &data, size_t n) {
    if (n >= m_meta.ntrace) {
      throw std::out_of_range("Index out of bound: " + std::to_string(n));
    }
    if (data.ndim() != 1 || data.size() != m_meta.nt) {
      throw std::runtime_error("Input data shape not match.");
    }

    const float *ptr = data.data();
    SegyRW::write_itrace(ptr, n);
  }

  void write_traces(const npfloat &data, size_t beg, size_t end, size_t tbeg,
                    size_t tend) {
    if (beg > end || end > m_meta.ntrace || tbeg > tend || tend > m_meta.nt) {
      throw std::out_of_range("Index out of bound.");
    }
    if (data.size() != (uint64_t)(end - beg) * (tend - tbeg)) {
      throw std::runtime_error("Input data size not match.");
    }

    const float *ptr = data.data();
    SegyRW::write_traces(ptr, beg, end, tbeg, tend);
  }

  void write_traces(const npfloat &data, const npint32 &index, size_t tbeg,
                    size_t tend) {
    if (index.ndim() != 1) {
      throw std::runtime_error("Input index must be a 1D data.");
    }
    if (tbeg > tend || tend > m_meta.nt) {
      throw std::out_of_range("`tbeg` or `tend` index out of bound.");
    }
    if (data.size() != index.size() * (tend - tbeg)) {
      throw std::runtime_error("Input data size not match.");
    }

    const float *ptr = data.data();
    const int32_t *idx = index.data();
    SegyRW::write_traces(ptr, idx, index.shape()[0], tbeg, tend);
  }

  void write(const npfloat &data) {
    if (ndim() == 2 && data.size() != m_meta.ntrace * m_meta.nt) {
      throw std::runtime_error("Input data size not match.");
    }
    if (ndim() == 3 &&
        data.size() != (uint64_t)m_meta.ni * m_meta.nx * m_meta.nt) {
      throw std::runtime_error("Input data size not match.");
    }
    if (ndim() == 4 && data.size() != (uint64_t)m_meta.ni * m_meta.nx *
                                          m_meta.no * m_meta.nt) {
      throw std::runtime_error("Input data size not match.");
    }

    const float *ptr = data.data();
    SegyRW::write(ptr);
  }

  void write3d(const npfloat &data, size_t ib, size_t ie, size_t xb, size_t xe,
               size_t tb, size_t te) {
    if (ndim() != 3) {
      throw std::runtime_error("write3d function valid when ndim == 3");
    }
    if (ib > ie || ie > m_meta.ni || xb > xe || xe > m_meta.nx || tb > te ||
        te > m_meta.nt) {
      throw std::out_of_range("Index out of bound.");
    }
    if (data.size() != (uint64_t)(ie - ib) * (xe - xb) * (te - tb)) {
      throw std::runtime_error("Input data size not match.");
    }

    const float *ptr = data.data();
    SegyRW::write3d(ptr, ib, ie, xb, xe, tb, te);
  }

  void write4d(const npfloat &data, size_t ib, size_t ie, size_t xb, size_t xe,
               size_t ob, size_t oe, size_t tb, size_t te) {
    if (ndim() != 4) {
      throw std::runtime_error("write4d function valid when ndim == 4");
    }
    if (ib > ie || ie > m_meta.ni || xb > xe || xe > m_meta.nx || ob > oe ||
        oe > m_meta.no || tb > te || te > m_meta.nt) {
      throw std::out_of_range("Index out of bound.");
    }
    if (data.size() !=
        (uint64_t)(ie - ib) * (xe - xb) * (oe - ob) * (te - tb)) {
      throw std::runtime_error("Input data size not match.");
    }

    const float *ptr = data.data();
    SegyRW::write4d(ptr, ib, ie, xb, xe, ob, oe, tb, te);
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
      // iline, xline_start, xline_end, trace_start, trace_end
      out = py::array_t<int>({(int)m_meta.ni, 5});
      int *ptr = out.mutable_data();
      std::fill(ptr, ptr + out.size(), -1);
      for (auto& linfo : m_iinfos) {
        ptr[0] = linfo.line;
        if (!(linfo.count == 0 && linfo.idx.size() == 0)) {
          ptr[1] = linfo.lstart;
          ptr[2] = linfo.lend;
          ptr[3] = linfo.itstart;
          ptr[4] = linfo.itend;
        }
        ptr += 5;
      }
    } else {
      // iline, xline, offset_start, offset_end, trace_start, trace_end
      out = py::array_t<int>({(int)m_meta.ni, (int)m_meta.nx, 6});
      int *ptr = out.mutable_data();
      std::fill(ptr, ptr + out.size(), -1);

      for (auto& linfo : m_iinfos) {
        size_t xs = (linfo.lstart - m_meta.start_xline) / m_keys.xstep;
        size_t xe = (linfo.lend - m_meta.start_xline) / m_keys.xstep;

        if (xs > 0) {
          for (size_t ix = 0; ix < xs; ix++) {
            ptr[0] = linfo.line;
            ptr[1] = linfo.lstart + ix * m_keys.xstep;
            ptr += 6;
          }
        }

        for (auto& xinfo : linfo.xinfos) {
          ptr[0] = linfo.line;
          ptr[1] = xinfo.line;
          if (!(xinfo.count == 0 && xinfo.idx.size() == 0)) {
            ptr[2] = xinfo.lstart;
            ptr[3] = xinfo.lend;
            ptr[4] = xinfo.itstart;
            ptr[5] = xinfo.itend;
          }
          ptr += 6;
        }

        if (xe < m_meta.nx) {
          for (size_t ix = xe; ix < m_meta.nx; ix++) {
            ptr[0] = linfo.line;
            ptr[1] = linfo.lend + ix * m_keys.xstep;
            ptr += 6;
          }
        }
      }
      if ((ptr - out.mutable_data()) != out.size()) {
        std::runtime_error("Error when get lineinfo");
      }
    }
    return out;
  }
};

void create_segy(const std::string &segyname, const npfloat &src,
                 const npint32 &keys, const std::string &textual,
                 const npuchar &bheader, const npuchar &theader) {
  if (textual.size() != kTextualHeaderSize) {
    throw std::runtime_error("textual header size must be 3200, but got " +
                             std::to_string(textual.size()));
  }
  if (bheader.size() != kBinaryHeaderSize) {
    throw std::runtime_error("binary header size must be 400, but got " +
                             std::to_string(bheader.size()));
  }
  if (theader.size() != kTraceHeaderSize) {
    throw std::runtime_error("trace header size must be 240, but got " +
                             std::to_string(theader.size()));
  }
  if (keys.ndim() != 2) {
    throw std::runtime_error("keys must be a 2D array.");
  }
  size_t keysize = keys.shape()[1];
  std::vector<size_t> shape(src.shape(), src.shape() + src.ndim());
  uint64_t ntrace = 1;
  for (size_t i = 0; i < shape.size() - 1; i++) {
    ntrace *= shape[i];
  }
  if (ntrace != keys.shape()[0]) {
    throw std::runtime_error("Input data size not match with keys.");
  }
  if (keysize < 4 || keysize > 5) {
    throw std::runtime_error("keys size must be 4 or 5");
  }
  if (shape.size() == 4 && keysize != 5) {
    throw std::runtime_error("keys size must be 5 when ndim == 4.");
  }

  const float *ptr = src.data();
  const int32_t *kptr = keys.data();
  const uchar *bptr = bheader.data();
  const uchar *tptr = theader.data();
  create_segy(segyname, ptr, kptr, shape, textual, bptr, tptr, keysize);
}

npfloat ieees_to_ibms(const npfloat &ieee_arr, bool is_litte_endian_input, bool is_big_endian_output) {
  std::vector<size_t> shape_vec(ieee_arr.shape(),
                                ieee_arr.shape() + ieee_arr.ndim());
  size_t size = ieee_arr.size();

  py::array_t<float> ibm_arr(shape_vec);

  const float *ieee_ptr = ieee_arr.data();
  float *ibm_ptr = ibm_arr.mutable_data();

  for (size_t i = 0; i < size; ++i) {
    ibm_ptr[i] = segy::ieee_to_ibm(ieee_ptr[i], is_litte_endian_input, is_big_endian_output);
  }

  return ibm_arr;
}

npfloat ibms_to_ieees(const npfloat &ibm_arr, bool is_big_endian) {
  std::vector<int> shape_vec(ibm_arr.shape(), ibm_arr.shape() + ibm_arr.ndim());

  size_t size = ibm_arr.size();

  py::array_t<float> ieee_arr(shape_vec);

  const float *ibm_ptr = ibm_arr.data();
  float *ieee_ptr = ieee_arr.mutable_data();

  for (size_t i = 0; i < size; ++i) {
    ieee_ptr[i] = segy::ibm_to_ieee(ibm_ptr[i], is_big_endian);
  }

  return ibm_arr;
}

PYBIND11_MODULE(_CXX_SEGY, m) {

  g_check_signals_callback = checkSignals;


  // set global variable
  m.def("set_progress_callback", [](py::function func) {
      g_progress_callback = [func](int current, int total) {
          py::gil_scoped_acquire acquire;
          func(current, total);
      };
  });
  m.def("set_global_show_progress", [](bool show) {
      g_show_progress = show;
  });

  py::class_<Pysegy>(m, "Pysegy")
      .def(py::init<const std::string &, bool>(), py::arg("segyname"),
           py::arg("write") = false)
      .def("close", &Pysegy::close_file)
      .def("show_progress", &Pysegy::show_progress, py::arg("show"))

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
           py::overload_cast<const npint32 &, size_t, size_t>(&Pysegy::collect),
           py::arg("index"), py::arg("tbeg"), py::arg("tend"))
      .def("collect",
           py::overload_cast<size_t, size_t, size_t, size_t>(&Pysegy::collect),
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
      .def("read_tslice", &Pysegy::read_tslice, py::arg("t"),
           py::arg("stepi") = 1, py::arg("stepx") = 1)
      .def("tofile", &Pysegy::tofile, py::arg("binary_out_name"),
           py::arg("as_2d") = false)
      .def("cut", &Pysegy::cut, py::arg("outname"), py::arg("ranges"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           py::overload_cast<const std::string &, const npfloat &,
                             const py::list &, bool, const std::string &>(
               &Pysegy::create_by_sharing_header),
           py::arg("segy_name"), py::arg("src"), py::arg("start"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           py::overload_cast<const std::string &, const std::string &,
                             const py::list &, const py::list &, bool,
                             const std::string &>(
               &Pysegy::create_by_sharing_header),
           py::arg("segy_name"), py::arg("src_file"), py::arg("shape"),
           py::arg("start"), py::arg("is2d") = false, py::arg("textual") = "")

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
           py::overload_cast<const npfloat &, size_t, size_t, size_t, size_t>(
               &Pysegy::write_traces),
           py::arg("data"), py::arg("beg"), py::arg("end"), py::arg("tbeg"),
           py::arg("tend"))
      .def("write_traces",
           py::overload_cast<const npfloat &, const npint32 &, size_t, size_t>(
               &Pysegy::write_traces),
           py::arg("data"), py::arg("index"), py::arg("tbeg"), py::arg("tend"))
      .def("write", &Pysegy::write, py::arg("data"))
      .def("write3d", &Pysegy::write3d, py::arg("data"), py::arg("ib"),
           py::arg("ie"), py::arg("xb"), py::arg("xe"), py::arg("tb"),
           py::arg("te"))
      .def("write4d", &Pysegy::write4d, py::arg("data"), py::arg("ib"),
           py::arg("ie"), py::arg("xb"), py::arg("xe"), py::arg("ob"),
           py::arg("oe"), py::arg("tb"), py::arg("te"))

      .def_property_readonly("ntrace", &Pysegy::ntrace)
      .def_property_readonly("nt", &Pysegy::nt)
      .def_property_readonly("ndim", &Pysegy::ndim)
      .def_property_readonly("shape", &Pysegy::shape)
      .def_property_readonly("is_scanned", &Pysegy::is_scanned)

      ;

  m.def(
      "create_segy",
      py::overload_cast<const std::string &, const npfloat &, const npint32 &,
                        const std::string &, const npuchar &, const npuchar &>(
          &create_segy),
      py::arg("segyname"), py::arg("src"), py::arg("keys"), py::arg("textual"),
      py::arg("bheader"), py::arg("theader"));
  m.def("ieee_to_ibm", py::overload_cast<float, bool, bool>(&segy::ieee_to_ibm),
        py::arg("value"), py::arg("is_litte_endian_input"), py::arg("is_big_endian_output")=true);
  m.def("ibm_to_ieee", py::overload_cast<float, bool>(&segy::ibm_to_ieee),
        py::arg("value"), py::arg("is_big_endian"));
  m.def("ieees_to_ibms", &ieees_to_ibms, py::arg("ieee_arr"),
        py::arg("is_litte_endian_input"), py::arg("is_big_endian_output"));
  m.def("ibms_to_ieees", &ibms_to_ieees, py::arg("ibm_arr"),
        py::arg("is_big_endian"));
}
} // namespace segy
