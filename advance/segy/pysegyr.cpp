#include "segynd.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using npfloat = py::array_t<float, py::array::c_style | py::array::forcecast>;
using npint32 = py::array_t<int32_t, py::array::c_style | py::array::forcecast>;

class segyr : public segy::SegyND {
public:
  using segy::SegyND::create_by_sharing_header;
  using segy::SegyND::SegyND;

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
    SegyND::read4d(ptr, is, ie, xs, xe, os, oe, ts, te);
    return data;
  }

  npfloat read3d(int is, int ie, int xs, int xe, int ts, int te) {
    int nt = te - ts;
    int ni = ie - is;
    int nx = xe - xs;
    auto data = py::array_t<float>({ni, nx, nt});
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyND::read3d(ptr, is, ie, xs, xe, ts, te);
    return data;
  }

  npfloat read() {
    py::array_t<float> out(shape());
    auto buff = out.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyND::read(ptr);
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
    SegyND::itrace(ptr, n);
    return out;
  }

  npfloat collect(int beg, int end, int tbeg, int tend) {
    // TODO: check bound
    auto data = py::array_t<float>({end - beg, tend - tbeg});
    auto buff = data.request();
    float *ptr = static_cast<float *>(buff.ptr);
    SegyND::collect(ptr, beg, end, tbeg, tend);
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
    SegyND::collect(ptr, idx, N, tbeg, tend);
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
    SegyND::get_trace_keys(ptr, keysvec, lengthvec, beg, end);
    return out;
  }
};

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(_CXX_SEGYR, m) {
  py::class_<segyr>(m, "segyr")
      .def(py::init<std::string>())

      // location
      .def("setLocations", &segyr::setLocations, py::arg("iline"),
           py::arg("xline"), py::arg("offset") = 37)
      .def("setInlineLocation", &segyr::setInlineLocation, py::arg("iline"))
      .def("setCrosslineLocation", &segyr::setCrosslineLocation,
           py::arg("xline"))
      .def("setOffsetLocation", &segyr::setOffsetLocation, py::arg("offset"))
      .def("setXLocation", &segyr::setXLocation, py::arg("xloc"))
      .def("setYLocation", &segyr::setYLocation, py::arg("yloc"))
      .def("setXYLocation", &segyr::setXYLocations, py::arg("xloc"),
           py::arg("yloc"))
      .def("setInlineStep", &segyr::setInlineStep, py::arg("istep"))
      .def("setCrosslineStep", &segyr::setCrosslineStep, py::arg("xstep"))
      .def("setOffsetStep", &segyr::setOffsetStep, py::arg("ostep"))
      .def("setSteps", &segyr::setSteps, py::arg("istep"), py::arg("xstep"),
           py::arg("ostep") = 1)
      .def("setFill", &segyr::setFill, py::arg("fill"))

      // read func base
      .def("textual_header", &segyr::textual_header, py::arg("coding") = 'u')
      .def("get_trace_keys", &segyr::get_trace_keys, py::arg("leys"),
           py::arg("length"), py::arg("beg"), py::arg("end"))
      .def("itrace", &segyr::itrace, py::arg("n"))
      .def("collect",
           overload_cast_<const npint32 &, int, int>()(&segyr::collect),
           py::arg("index"), py::arg("tbeg"), py::arg("tend"))
      .def("collect", overload_cast_<int, int, int, int>()(&segyr::collect),
           py::arg("beg"), py::arg("end"), py::arg("tbeg"), py::arg("tend"))

      // read
      .def("set_segy_type", &segyr::set_segy_type, py::arg("ndim"))
      .def("scan", &segyr::scan)
      .def("read4d", &segyr::read4d, py::arg("ib"), py::arg("ie"),
           py::arg("xb"), py::arg("xe"), py::arg("ob"), py::arg("oe"),
           py::arg("tb"), py::arg("te"))
      .def("read3d", &segyr::read3d, py::arg("ib"), py::arg("ie"),
           py::arg("xb"), py::arg("xe"), py::arg("tb"), py::arg("te"))
      .def("read", &segyr::read)
      .def("tofile", &segyr::tofile, py::arg("binary_out_name"),
           py::arg("as_2d") = false)
      .def("cut", &segyr::cut, py::arg("outname"), py::arg("ranges"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           overload_cast_<const std::string &, const npfloat &,
                          const py::list &, bool, const std::string &>()(
               &segyr::create_by_sharing_header),
           py::arg("segy_name"), py::arg("src"), py::arg("ranges"),
           py::arg("is2d") = false, py::arg("textual") = "")
      .def("create_by_sharing_header",
           overload_cast_<const std::string &, const std::string &,
                          const std::vector<int> &, const std::vector<int> &,
                          bool, const std::string &>()(
               &segyr::create_by_sharing_header))

      ;
}
