#ifndef CIG_SEGY_nd_H
#define CIG_SEGY_nd_H
#define FMT_HEADER_ONLY

#include "mio.hpp"
#include "progressbar.hpp"
#include "segybase.hpp"
#include "sutils.h"

#include <fmt/format.h>
#include <fstream>
#include <stdexcept>
#include <vector>

#define MMAX(a, b) ((a) > (b) ? (a) : (b))
#define MMIN(a, b) ((a) < (b) ? (a) : (b))

namespace segy {

struct XLineInfo {
  int xline;
  int xtstart; // xline start trace index
  int xtend;   // xlind end trace index
  int count;   // the number of traces for this xline
  int ostart;  // offset start index
  int oend;    // offset end index

  // offset trace index of this xline,
  // only used when this xline is not continous, such as
  // 2, 3, 4, x, x, 7, 8, 9. oend - ostart + 1 = 8, but count=6
  std::vector<int> oidx;
};

struct LineInfo {
  int iline;
  int itstart; // inline trace start index
  int itend;   // inline trace end index
  int count;   // the number of traces for this line
  int xstart;  // xline start index
  int xend;    // xline end index
  std::vector<XLineInfo> xinfos;

  // xline trace index of this line,
  // only used when this line is not continous and is 3D, such as
  // 2, 3, 4, x, x, 7, 8, 9. xend - xstart + 1 = 8, but count=6
  std::vector<int> xidx;
};

class SegyND : public SegyBase {
public:
  explicit SegyND(const std::string &segyname);
  ~SegyND() override;

  void set_segy_type(int ndim);
  void scan();
  std::vector<int> shape();
  void read4d(float *dst, int is, int ie, int xs, int xe, int os, int oe,
              int ts, int te);
  void read3d(float *dst, int is, int ie, int xs, int xe, int ts, int te);
  void read(float *dst);
  void tofile(const std::string &binary_out_name, bool is2d = false);
  void cut(const std::string &outname, const std::vector<int> &ranges,
           bool is2d = false, const std::string &textual = "");
  void create_by_sharing_header(const std::string &segy_name, const float *src,
                                const std::vector<int> &ranges,
                                bool is2d = false,
                                const std::string &textual = "");
  void create_by_sharing_header(const std::string &segy_name,
                                const std::string &src_name,
                                const std::vector<int> &shape,
                                const std::vector<int> &ranges,
                                bool is2d = false,
                                const std::string &textual = "");

protected:
  mio::mmap_source m_src;
  std::vector<LineInfo> m_iinfos;
  int m_ndim = 2;
  void scanBinaryHeader();

private:
  bool isScan = false;
  uint64_t m_sizeL = 1;

  bool isPreStack();
  inline int il2ii(int il) { return (il - m_meta.start_iline) / m_keys.istep; }
  inline int ii2il(int ii) { return ii * m_keys.istep + m_meta.start_iline; }
  inline int xl2ix(int xl) { return (xl - m_meta.start_xline) / m_keys.xstep; }
  inline int ix2xl(int ix) { return ix * m_keys.xstep + m_meta.start_xline; }
  inline int of2io(int of) { return (of - m_meta.start_offset) / m_keys.ostep; }
  inline int io2of(int io) { return io * m_keys.ostep + m_meta.start_offset; }
  bool isCtnX(LineInfo &linfo);
  bool isCtnO(XLineInfo &xinfo);
  void find_idxX(std::array<int32_t, 4> &idx, LineInfo &linfo, int xs, int xe);
  void find_idxO(std::array<int32_t, 4> &idx, XLineInfo &xinfo, int os, int oe);

  void _read3d_x(float *dst, LineInfo &linfo, int xs, int xe, int ts, int te);
  void _read4d_o(float *dst, XLineInfo &xinfo, int os, int oe, int ts, int te);
  void _read4d_xo(float *dst, LineInfo &linfo, int xs, int xe, int os, int oe,
                  int ts, int te);
  uint64_t _copy3d_x(char *dst, const float *src, LineInfo &linfo, int xs,
                     int xe, int ts, int te, bool fromsrc);
  uint64_t _copy4d_o(char *dst, const float *src, XLineInfo &xinfo, int os,
                     int oe, int ts, int te, bool fromsrc);
  uint64_t _copy4d_xo(char *dst, const float *src, LineInfo &linfo, int xs,
                      int xe, int os, int oe, int ts, int te, bool fromsrc);
  void _create_from_segy(const std::string &outname, const float *src,
                         const std::vector<int> &ranges, bool is2d,
                         const std::string &textual, bool fromsrc);

  uint64_t _need_size_b(int ndim = -1);
  uint64_t _need_size_s(uint64_t ntrace);
};

inline SegyND::SegyND(const std::string &segyname) : SegyBase(segyname) {
  std::error_code error;
  this->m_src.map(segyname, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }

  m_data_ptr = m_src.data();
  this->scanBinaryHeader();
}

inline SegyND::~SegyND() {
  if (m_src.is_mapped()) {
    m_src.unmap();
  }
}

inline void SegyND::set_segy_type(int ndim) {
  if (ndim < 2 || ndim > 4) {
    std::runtime_error(
        "Error SEG-Y type, 2 for 2D, 3 for poststack, 4 for prestack");
  }
  m_ndim = ndim;
}

inline void SegyND::scan() {
  m_meta.start_time = keyi2(0, kTStartTimeField);
  m_meta.scalar = keyi2(0, kTScalarField);

  // steps
  int istep = m_keys.istep;
  int xstep = m_keys.xstep;
  int ostep = m_keys.ostep;

  // inline range
  int is = iline(0);
  int ie = iline(m_meta.ntrace - 1);
  int ni = (is - ie) / m_keys.istep + 1;
  m_iinfos.resize(ni);

  // global xline range and offset range
  int gxstart = xline(0);
  int gxend = xline(m_meta.ntrace - 1);
  int gostart = offset(0);
  int goend = offset(m_meta.ntrace - 1);

  // is a 4D pretack SEG-Y?
  bool is4D = true;
  if (m_ndim == 2) {
    is4D = isPreStack();
  } else if (m_ndim == 3) {
    is4D = false;
  }

  // use skipi and skipx to process when missing line and xline (in 4D)
  int skipi = 0;
  int skipx = 0;

  // trace index for iline
  int it = 0;
  int jumpl = 1;
  int jumpx = 1;
  for (size_t ii = 0; ii < ni; ++ii) {
    CHECK_SIGNALS();
    LineInfo &linfo = m_iinfos[ii];

    // when missing line
    if (skipi > 0) {
      linfo.iline = m_iinfos[ii - 1].iline + istep;
      linfo.count = 0;
      skipi--;
      continue;
    }

    // record the number of this line, trace start and xline start
    int iiline = iline(it);
    int itstart = it;
    int xs = xline(it);

    // jump
    it += jumpl;

    // jump too small
    if (iline(it) == iiline) {
      while (iline(it) == iiline && it < m_meta.ntrace) {
        jumpl++;
        it++;
      }
    }
    // jump too large
    else if (iline(it - 1) != iiline) {
      while (iline(it - 1) != iiline && it >= itstart) {
        it--;
        jumpl--;
      }
    }

    if (iline(it) == iiline || iline((it - 1) != iiline)) {
      throw std::runtime_error("error1");
    }

    // assign lineInfo
    int xend = xline(it - 1);
    linfo.iline = iiline;
    linfo.itstart = itstart;
    linfo.itend = it - 1;
    linfo.xstart = xs;
    linfo.xend = xend;
    linfo.count = it - itstart;

    // update global xline range
    gxstart = xstep > 0 ? MMIN(gxstart, xs) : MMAX(gxstart, xs);
    gxend = xstep > 0 ? MMAX(gxend, xend) : MMIN(gxend, xend);

    if (is4D) {
      // assert (m_iinfos[ii].xend - xs) % xstep == 0
      int nx = (xend - xs) / xstep + 1;
      linfo.xinfos.resize(nx);
      int xt = itstart;
      int xtmax = it;

      for (size_t ix = 0; ix < nx; ix++) {
        XLineInfo &xinfo = linfo.xinfos[ix];
        // when missing xline
        if (skipx > 0) {
          xinfo.xline = linfo.xinfos[ix - 1].xline + xstep;
          xinfo.count = 0;
          skipx--;
          continue;
        }

        int xxline = iline(xt);
        int xtstart = xt;
        int ostart = offset(xt);

        xt += jumpx;

        // jump too small
        if (xline(xt) == xxline) {
          while (xline(xt) == xxline && xt < xtmax) {
            jumpx++;
            xt++;
          }
        }
        // jump too large
        else if (xline(xt - 1) != xxline) {
          while (xline(xt - 1) != xxline && xt >= itstart) {
            xt--;
            jumpx--;
          }
        }

        if (xline(xt) == xxline || xline((xt - 1) != xxline)) {
          throw std::runtime_error("error2");
        }

        // assign xline info
        int oend = offset(xt - 1);
        xinfo.xline = xxline;
        xinfo.xtstart = xtstart;
        xinfo.xtend = xt - 1;
        xinfo.ostart = ostart;
        xinfo.oend = oend;
        xinfo.count = xt - xtstart;

        // update global offset range
        gostart = ostep > 0 ? MMIN(gostart, ostart) : MMAX(gostart, ostart);
        goend = ostep > 0 ? MMAX(goend, oend) : MMIN(goend, oend);

        // missing xline
        if (xline(xt) != (xxline + xstep)) {
          skipx = xline(xt) - (xxline + xstep);
          if (skipx % xstep != 0) {
            throw std::runtime_error(
                "Cannot analysis this segy file, may inline step != 1");
          }
          skipx /= xstep;
        }
      }
    }

    // missing line
    if (xline(it) != (iiline + istep)) {
      skipi = xline(it) - (iiline + istep);
      if (skipi % istep != 0) {
        throw std::runtime_error(
            "Cannot analysis this segy file, may inline step != 1");
      }
      skipi /= istep;
    }
  }

  // Post process, assign m_meta
  m_meta.start_iline = is;
  m_meta.end_iline = ie;
  m_meta.ni = ni;
  m_meta.start_xline = gxstart;
  m_meta.end_xline = gxend;
  m_meta.nx = (gxend - gxstart) / xstep + 1;

  if (is4D) {
    m_meta.start_offset = gostart;
    m_meta.end_offset = goend;
    m_meta.no = (goend - gostart) / ostep + 1;
  }

  m_sizeL = m_meta.nx * m_meta.no;

  /********** if line or xline is not continouse, we record their idx for fast
   * indexing ************/
  for (auto linfo : m_iinfos) {
    if (!is4D) {
      if (linfo.count == 0) {
        continue;
      }
      if (linfo.count == (linfo.xend - linfo.xstart) / m_keys.xstep + 1) {
        continue;
      }

      // we set count = -1 for not continous line
      linfo.xidx.resize(m_meta.nx, -1);
      for (size_t xt = linfo.itstart; xt < linfo.itend + 1; xt++) {
        linfo.xidx[xl2ix(xline(xt))] = xt;
      }
      linfo.count = -1;
    } else {
      for (auto xinfo : linfo.xinfos) {
        if (xinfo.count == 0) {
          continue;
        }
        if (xinfo.count == (xinfo.oend - xinfo.ostart) / m_keys.ostep + 1) {
          continue;
        }

        // we set count = -1 for not continous xline
        xinfo.oidx.resize(m_meta.no, -1);
        for (size_t ot = xinfo.xtstart; ot < xinfo.xtend + 1; ot++) {
          xinfo.oidx[of2io(offset(ot))] = ot;
        }
        xinfo.count = -1;
      }
    }
  }

  { /********** calculate interval **********/

    // cal x, y interval
    // find a line contained more than 5 traces
    int s = m_meta.ni / 8;
    while (m_iinfos[s].count < 6 && s < m_meta.ni) {
      s++;
    }

    // Y_interval = ((x2-x1)^2 + (y2-y1)^2)^0.5 / abs(xline2 - xline1) * xstep
    int te = m_iinfos[s].itend;
    int ts = m_iinfos[s].itstart;
    m_meta.dx = std::sqrt(std::pow(coordx(te) - coordx(ts), 2) +
                          std::pow(coordy(te) - coordy(ts), 2)) /
                std::abs(xline(te) - xline(ts)) * xstep;

    s = m_meta.ni - 1;
    te = m_iinfos[s].itstart;

    // S1 = ((x2-x1)^2 + (y2-y1)^2)^0.5
    // S2 = abs(xline2 - xline1) / xstep * Y_interval
    // (S1^2 - S2^2)^0.5 = S3
    // Z_interval = S3 / abs(iline2-iline1) / istep
    float S1_2 = std::pow(coordx(te) - coordx(ts), 2) +
                 std::pow(coordy(te) - coordy(ts), 2);
    float S2_2 = std::pow(float(xline(te) - xline(ts)) / xstep * m_meta.dx, 2);
    m_meta.di = std::sqrt(S1_2 - S2_2) /
                (std::abs(iline(te) - iline(ts)) / float(istep));
  }

  isScan = true;
  m_ndim = 3;
  if (is4D) {
    m_ndim = 4;
  }
}

inline std::vector<int> SegyND::shape() {
  if (m_ndim == 2) {
    return {(int)m_meta.ntrace, m_meta.nt};
  } else if (m_ndim == 3) {
    return {m_meta.ni, m_meta.nx, m_meta.nt};
  } else {
    return {m_meta.ni, m_meta.nx, m_meta.no, m_meta.nt};
  }
}

inline void SegyND::read4d(float *dst, int is, int ie, int xs, int xe, int os,
                           int oe, int ts, int te) {
  if (!isScan || m_ndim != 4) {
    throw std::runtime_error("Not scan or is not a 4d pretrack SEG-Y");
  }

  int ni = ie - is;
  int nx = xe - xs;
  int no = oe - os;
  int nt = te - ts;
  uint64_t sizeOT = no * nt;
  uint64_t sizeXOT = nx * sizeOT;

  for (size_t ii = is; ii < ie; ii++) {
    CHECK_SIGNALS();
    // if (showpbar) {
    //   updatebar(bar, step);
    // }

    LineInfo &linfo = m_iinfos[ii];
    float *dstiline = dst + (ii - is) * sizeXOT;
    _read4d_xo(dstiline, linfo, xs, xe, os, oe, ts, te);
  }
}

inline void SegyND::read3d(float *dst, int is, int ie, int xs, int xe, int ts,
                           int te) {
  if (!isScan || m_ndim != 3) {
    throw std::runtime_error("Not scan or is not a 3d SEG-Y");
  }

  if (is >= ie || xs >= xe || ts >= te) {
    throw std::runtime_error("Index 'end' must large than 'start'");
  }
  if (is < 0 || ie > m_meta.ni || xs < 0 || xe > m_meta.nx || ts < 0 ||
      te > m_meta.nt) {
    throw std::runtime_error("Index out of range");
  }

  int ni = ie - is;
  int nx = xe - xs;
  int nt = te - ts;
  uint64_t sizeXT = nx * nt;

  // int step = ni >= 100 ? 1 : 100 / ni + 1;
  // int nbar = ni >= 100 ? ni : step * ni;
  // progressbar bar(nbar);

  for (size_t ii = is; ii < ie; ii++) {
    CHECK_SIGNALS();
    // if (showpbar) {
    //   updatebar(bar, step);
    // }

    LineInfo &linfo = m_iinfos[ii];
    float *dstiline = dst + (ii - is) * sizeXT;
    _read3d_x(dstiline, linfo, xs, xe, ts, te);
  }
}

inline void SegyND::read(float *dst) {
  if (m_ndim == 2) {
    collect(dst, 0, m_meta.ntrace, 0, m_meta.nt);
  } else if (m_ndim == 3) {
    read3d(dst, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.nt);
  } else if (m_ndim == 4) {
    read4d(dst, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.no, 0, m_meta.nt);
  }
}

inline void SegyND::tofile(const std::string &binary_out_name, bool is2d) {
  if (!is2d && m_ndim == 2) {
    throw std::runtime_error("Not a 3D or 4D SEG-Y");
  }
  if (!isScan && !is2d) {
    scan();
  }

  uint64_t need_size = 0;
  if (!is2d) {
    need_size = _need_size_b(m_ndim);
  } else {
    need_size = _need_size_b(2);
  }

  create_file(binary_out_name, need_size);

  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(binary_out_name, error);
  if (error) {
    throw std::runtime_error("mmap fail when write data");
  }

  float *dst = reinterpret_cast<float *>(rw_mmap.data());

  // or need split into serveral chunks?
  if (!is2d) {
    read(dst);
  } else {
    collect(dst, 0, (int)m_meta.ntrace, 0, m_meta.nt);
  }

  rw_mmap.unmap();
}

/*   priviate function for reading */
inline void SegyND::_read3d_x(float *dst, LineInfo &linfo, int xs, int xe,
                              int ts, int te) {
  int nt = te - ts;
  int nx = xe - xs;
  uint64_t sizeXT = nt * nx;

  /* no data to copy, just fill */
  if (linfo.count == 0 || xs > xl2ix(linfo.xend) || xe <= xl2ix(linfo.xstart)) {
    std::fill(dst, dst + sizeXT, m_meta.fillNoValue);
    return;
  }

  /* line is not continuous */
  else if (linfo.count == -1) {
    for (size_t ix = xs; ix < xe; ix++) {
      int it = linfo.xidx[ix];
      if (it == -1) {
        std::fill(dst, dst + nt, m_meta.fillNoValue);
      } else {
        m_readfunc(dst, trDataStart(it) + ts, nt);
      }
      dst += nt;
    }
    return;
  }

  /* continuous line TODO: remove if confition?*/
  else if (linfo.count > 0 && isCtnX(linfo)) {
    std::array<int32_t, 4> idx;
    find_idxX(idx, linfo, xs, xe);

    if (idx[0] > 0) {
      std::fill(dst, dst + idx[0] * nt, m_meta.fillNoValue);
      dst += idx[0] * nt;
    }
    collect(dst, idx[2], idx[3], ts, te);
    dst += (idx[3] - idx[2]) * nt;

    if (idx[1] > 0) {
      std::fill(dst, dst + idx[1] * nt, m_meta.fillNoValue);
    }

    // assert
    return;
  } else {
    throw std::runtime_error("error count");
  }
}

inline void SegyND::_read4d_o(float *dst, XLineInfo &xinfo, int os, int oe,
                              int ts, int te) {
  int nt = te - ts;
  int no = oe - os;
  uint64_t sizeOT = nt * no;

  /* no data to copy, just fill */
  if (xinfo.count == 0 || os > of2io(xinfo.oend) || oe <= of2io(xinfo.ostart)) {
    std::fill(dst, dst + sizeOT, m_meta.fillNoValue);
    return;
  }

  /* xline is not continuous */
  else if (xinfo.count == -1) {
    for (size_t io = os; io < oe; io++) {
      int it = xinfo.oidx[io];
      if (it == -1) {
        std::fill(dst, dst + nt, m_meta.fillNoValue);
      } else {
        m_readfunc(dst, trDataStart(it) + ts, nt);
      }
      dst += nt;
    }
    return;
  }

  /* continuous xline */
  else if (xinfo.count > 0 && isCtnO(xinfo)) {

    std::array<int32_t, 4> idx;
    find_idxO(idx, xinfo, os, oe);

    if (idx[0] > 0) {
      std::fill(dst, dst + idx[0] * nt, m_meta.fillNoValue);
      dst += idx[0] * nt;
    }
    collect(dst, idx[2], idx[3], ts, te);
    dst += (idx[3] - idx[2]) * nt;

    if (idx[1] > 0) {
      std::fill(dst, dst + idx[1] * nt, m_meta.fillNoValue);
    }

    // assert
    return;
  }
}

inline void SegyND::_read4d_xo(float *dst, LineInfo &linfo, int xs, int xe,
                               int os, int oe, int ts, int te) {
  int nt = te - ts;
  int no = oe - os;
  int nx = xe - xs;
  uint64_t sizeOT = nt * no;
  uint64_t sizeXOT = nx * sizeOT;

  /* no data to copy, just fill */
  if (linfo.count == 0 || xs > xl2ix(linfo.xend) || xe <= xl2ix(linfo.xstart)) {
    std::fill(dst, dst + sizeXOT, m_meta.fillNoValue);
    return;
  }

  std::array<int32_t, 4> idx;
  find_idxX(idx, linfo, xs, xe);

  // fill the left
  if (idx[0] > 0) {
    std::fill(dst, dst + idx[0] * sizeOT, m_meta.fillNoValue);
    dst += idx[0] * sizeOT;
  }

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = idx[2]; ix < idx[3]; ix++) {
    _read4d_o(dst, linfo.xinfos[ix], os, oe, ts, te);
    dst += sizeOT;
  }

  // fill the right
  if (idx[1] > 0) {
    std::fill(dst, dst + idx[1] * sizeOT, m_meta.fillNoValue);
    dst += idx[1] * sizeOT;
  }
  // assert
}

inline uint64_t SegyND::_copy3d_x(char *dst, const float *src, LineInfo &linfo,
                                  int xs, int xe, int ts, int te,
                                  bool fromsrc) {
  char *odst = dst;
  int nt = te - ts;
  int nx = xe - xs;
  bool tchanged = nt == m_meta.nt ? false : true;
  uint64_t sizeXT = nt * nx;

  /* no data to copy, just fill */
  if (linfo.count == 0 || xs > xl2ix(linfo.xend) || xe <= xl2ix(linfo.xstart)) {
    return 0;
  }

  /* line is not continuous */
  else if (linfo.count == -1) {
    for (size_t ix = xs; ix < xe; ix++) {
      int it = linfo.xidx[ix];
      if (it >= 0) {
        // copy trace header
        memcpy(dst, trheader(it), kTraceHeaderSize);
        if (tchanged) {
          if (ts > 0) { // TODO: ts need *1000 ?
            set_keyi2(dst, kTStartTimeField, ts);
          }
          set_keyi2(dst, kTSampleCountField, nt);
        }
        dst += kTraceHeaderSize;

        // copy trace
        if (fromsrc) {
          // copy from src, i.e., cut
          memcpy(dst, trDataStart(it) + ts, nt * m_meta.esize);
        } else {
          // copy from data, i.e., create_by_sharing_header
          m_wfunc(dst, src + (ix - xs) * nt, nt);
        }
        dst += nt * m_meta.esize;
      }
    }
    return dst - odst;
  }

  /* continuous line TODO: remove if confition?*/
  else if (linfo.count > 0 && isCtnX(linfo)) {
    std::array<int32_t, 4> idx;
    find_idxX(idx, linfo, xs, xe);
    const float *srcf = src + idx[0] * nt;

    for (size_t tx = idx[2]; tx < idx[3]; tx++) {
      // copy trace header
      memcpy(dst, trheader(tx), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) { // TODO: ts need *1000 ?
          set_keyi2(dst, kTStartTimeField, ts);
        }
        set_keyi2(dst, kTSampleCountField, nt);
      }
      dst += kTraceHeaderSize;

      // copy data
      if (fromsrc) {
        // copy from src, i.e., cut
        memcpy(dst, trDataStart(tx) + ts, nt * m_meta.esize);
      } else {
        // copy from data, i.e., create_by_sharing_header
        m_wfunc(dst, srcf + (tx - idx[2]) * nt, nt);
      }
      dst += nt * m_meta.esize;
    }
    // assert
    return dst - odst;
  } else {
    throw std::runtime_error("error count");
  }
}

inline uint64_t SegyND::_copy4d_o(char *dst, const float *src, XLineInfo &xinfo,
                                  int os, int oe, int ts, int te,
                                  bool fromsrc) {
  char *odst = dst;
  int nt = te - ts;
  int no = oe - os;
  bool tchanged = nt == m_meta.nt ? false : true;
  uint64_t sizeOT = nt * no;

  /* no data to copy, just fill */
  if (xinfo.count == 0 || os > of2io(xinfo.oend) || oe <= of2io(xinfo.ostart)) {
    return 0;
  }

  /* xline is not continuous */
  else if (xinfo.count == -1) {
    for (size_t io = os; io < oe; io++) {
      int it = xinfo.oidx[io];
      if (it >= 0) {
        // copy trace header
        memcpy(dst, trheader(it), kTraceHeaderSize);
        if (tchanged) {
          if (ts > 0) { // TODO: ts need *1000 ?
            set_keyi2(dst, kTStartTimeField, ts);
          }
          set_keyi2(dst, kTSampleCountField, nt);
        }
        dst += kTraceHeaderSize;

        // copy trace
        if (fromsrc) {
          // copy from src, i.e., cut
          memcpy(dst, trDataStart(it) + ts, nt * m_meta.esize);
        } else {
          // copy from data, i.e., create_by_sharing_header
          m_wfunc(dst, src + (io - os) * nt, nt);
        }
        dst += nt * m_meta.esize;
      }
    }
    return dst - odst;
  }

  /* continuous xline */
  else if (xinfo.count > 0 && isCtnO(xinfo)) {

    std::array<int32_t, 4> idx;
    find_idxO(idx, xinfo, os, oe);
    const float *srcf = src + idx[0] * nt;

    for (size_t to = idx[2]; to < idx[3]; to++) {
      // copy trace header
      memcpy(dst, trheader(to), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) { // TODO: ts need *1000 ?
          set_keyi2(dst, kTStartTimeField, ts);
        }
        set_keyi2(dst, kTSampleCountField, nt);
      }
      dst += kTraceHeaderSize;

      // copy data
      if (fromsrc) {
        // copy from src, i.e., cut
        memcpy(dst, trDataStart(to) + ts, nt * m_meta.esize);
      } else {
        // copy from data, i.e., create_by_sharing_header
        m_wfunc(dst, srcf + (to - idx[2]) * nt, nt);
      }
      dst += nt * m_meta.esize;
    }

    // assert
    return dst - odst;
  } else {
    throw std::runtime_error("error");
  }
}

inline uint64_t SegyND::_copy4d_xo(char *dst, const float *src, LineInfo &linfo,
                                   int xs, int xe, int os, int oe, int ts,
                                   int te, bool fromsrc) {
  char *odst = dst;
  int nt = te - ts;
  int no = oe - os;
  int nx = xe - xs;
  bool tchanged = nt == m_meta.nt ? false : true;
  uint64_t sizeOT = nt * no;
  uint64_t sizeXOT = nx * sizeOT;

  /* no data to copy, just fill */
  if (linfo.count == 0 || xs > xl2ix(linfo.xend) || xe <= xl2ix(linfo.xstart)) {
    return 0;
  }

  uint64_t jump = 0;

  std::array<int32_t, 4> idx;
  find_idxX(idx, linfo, xs, xe);

  const float *srcf = src + idx[0] * sizeXOT;

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = idx[2]; ix < idx[3]; ix++) {
    jump = _copy4d_o(dst, srcf + (ix - idx[2]) * sizeXOT, linfo.xinfos[ix], os,
                     oe, ts, te, fromsrc);
    dst += jump;
  }

  return dst - odst;
  // assert
}

inline bool SegyND::isPreStack() {
  int o0 = offset(0);
  int n = m_meta.ntrace < 500 ? m_meta.ntrace : 500;
  for (size_t i = 1; i < n; i++) {
    if (o0 != offset(i)) {
      return true;
    }
  }
  return false;
}

inline bool SegyND::isCtnX(LineInfo &linfo) {
  return linfo.count == ((linfo.xend - linfo.xstart) / m_keys.xstep + 1);
}
inline bool SegyND::isCtnO(XLineInfo &xinfo) {
  return xinfo.count == ((xinfo.oend - xinfo.ostart) / m_keys.ostep + 1);
}

inline uint64_t SegyND::_need_size_b(int ndim) {
  uint64_t need_size = 0;
  if (ndim == -1) {
    ndim = m_ndim;
  }

  // TODO: check?

  if (ndim == 2) {
    need_size = m_meta.ntrace * m_meta.nt;
  } else if (ndim == 3) {
    need_size = static_cast<uint64_t>(m_meta.ni) * m_meta.nx * m_meta.nt;
  } else if (ndim == 4) {
    need_size =
        static_cast<uint64_t>(m_meta.ni) * m_meta.nx * m_meta.no * m_meta.nt;
  }

  return need_size * sizeof(float);
}

inline uint64_t SegyND::_need_size_s(uint64_t ntrace) {
  return kTraceHeaderStart + m_meta.tracesize * ntrace;
}

// Use only if xs/xe and xinfo have overlapped parts
inline void SegyND::find_idxX(std::array<int32_t, 4> &idx, LineInfo &linfo,
                              int xs, int xe) {
  memset(&idx, 0, 4 * 4);
  int start = xl2ix(linfo.xstart);
  int end = xl2ix(linfo.xend) + 1;
  if (xs < start) {
    idx[0] = start - xs;
    xs = start;
  }
  if (xe > end) {
    idx[1] = xe - end;
    xe = end;
  }
  idx[2] = linfo.itstart + (xs - start);
  idx[3] = idx[2] + (xe - xs);
}

// Use only if os/oe and xinfo have overlapped parts
inline void SegyND::find_idxO(std::array<int32_t, 4> &idx, XLineInfo &xinfo,
                              int os, int oe) {
  memset(&idx, 0, 4 * 4);
  int start = of2io(xinfo.ostart);
  int end = of2io(xinfo.oend);
  if (os < start) {
    idx[0] = start - os;
    os = start;
  }
  if (oe > end) {
    idx[1] = oe - end;
    oe = end;
  }
  idx[2] = xinfo.xtstart + (os - start);
  idx[3] = idx[2] + (oe - os);
}

inline void SegyND::scanBinaryHeader() {
  int dformat = bkeyi2(kBSampleFormatField);
  auto it = kElementSize.find(dformat);
  if (it != kElementSize.end()) {
    m_meta.esize = it->second;
  } else {
    throw std::runtime_error(fmt::format("Unknown data format {}.\n", dformat));
  }
  m_meta.dformat = dformat;

  m_meta.nt = bkeyi2(kBSampleCountField);
  m_meta.tracesize = kTraceHeaderSize + m_meta.nt * m_meta.esize;
  m_meta.dt = bkeyi2(kBSampleIntervalField);
  m_meta.ntrace = (m_src.size() - kTraceHeaderStart) / m_meta.tracesize;

  m_meta.trace_sorting_code = bkeyi2(kBTraceSortingCodeField);
  setRWFunc(m_meta.dformat);
}

inline void SegyND::_create_from_segy(const std::string &outname,
                                      const float *src,
                                      const std::vector<int> &ranges, bool is2d,
                                      const std::string &textual,
                                      bool fromsrc) {

  int tstart, tend, is, ie, xs, xe, os, oe, ts, te;
  uint64_t maxsize = 0, tracesize = 0;
  if (is2d || m_ndim == 2) {
    if (ranges.size() != 4) {
      throw std::runtime_error("error ranges");
    }
    tstart = ranges[0];
    tend = ranges[1];
    ts = ranges[2];
    te = ranges[3];
    if (tstart < 0 || tend > m_meta.ntrace) {
      throw std::runtime_error("error ranges");
    }
    tracesize = kTraceHeaderSize + (te - ts) * m_meta.esize;
    maxsize = kTraceHeaderStart + tracesize * (tend - tstart);

  } else if (m_ndim == 3) {
    if (ranges.size() != 6) {
      throw std::runtime_error("error ranges");
    }
    is = ranges[0];
    ie = ranges[1];
    xs = ranges[2];
    xe = ranges[3];
    ts = ranges[4];
    te = ranges[5];
    tracesize = kTraceHeaderSize + (te - ts) * m_meta.esize;
    maxsize = kTraceHeaderStart + tracesize * (ie - is) * (xe - xs);
  } else if (m_ndim == 4) {
    if (ranges.size() != 8) {
      throw std::runtime_error("error ranges");
    }
    is = ranges[0];
    ie = ranges[1];
    xs = ranges[2];
    xe = ranges[3];
    os = ranges[4];
    oe = ranges[5];
    ts = ranges[6];
    te = ranges[7];
    if (os < 0 || oe > m_meta.no) {
      throw std::runtime_error("error ranges");
    }
    tracesize = kTraceHeaderSize + (te - ts) * m_meta.esize;
    maxsize = kTraceHeaderStart + tracesize * (ie - is) * (xe - xs) * (oe - os);
  }

  if (ts < 0 || te > m_meta.nt) {
    throw std::runtime_error("error ranges");
  }
  if (!is2d and m_ndim > 2) {
    if (is < 0 || ie > m_meta.ni) {
      throw std::runtime_error("error ranges");
    }
    if (xs < 0 || xe > m_meta.nx) {
      throw std::runtime_error("error ranges");
    }
  }

  // set dimension
  int nt = te - ts;
  int no = oe - os;
  int nx = xe - xs;
  int ni = ie - is;

  bool tchanged = nt == m_meta.nt ? false : true;

  // We first create a large file (assume there is no missing)
  // to store all data, and then mmap it to
  create_file(outname, maxsize);

  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(outname, error);
  if (error) {
    throw std::runtime_error("mmap fail when write data");
  }

  char *outptr = rw_mmap.data();

  // copy textual header
  if (textual.size() == 0) {
    memcpy(outptr, m_data_ptr, kTextualHeaderSize);
  } else if (textual.size() == kTextualHeaderSize) {
    memcpy(outptr, textual.data(), kTextualHeaderSize);
  } else {
    throw std::runtime_error("textual header size not match");
  }
  outptr += kTextualHeaderSize;

  // copy binary header
  memcpy(outptr, brheader(), kBinaryHeaderSize);
  if (tchanged) {
    set_bkeyi2(outptr, kBSampleCountField, te - ts);
  }
  outptr += kBinaryHeaderSize;

  // copy trace
  if (is2d || m_ndim == 2) {
    for (size_t it = tstart; it < tend; it++) {
      CHECK_SIGNALS();
      memcpy(outptr, trheader(it), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) { // TODO: ts need *1000 ?
          set_keyi2(outptr, kTStartTimeField, ts);
        }
        set_keyi2(outptr, kTSampleCountField, nt);
      }
      outptr += kTraceHeaderSize;
      memcpy(outptr, trDataStart(it) + ts, nt * m_meta.esize);
      outptr += nt * m_meta.esize;
    }
  } else {
    uint64_t jump = 0;
    for (size_t ii = is; ii < ie; ii++) {
      CHECK_SIGNALS();
      LineInfo &linfo = m_iinfos[ii];
      if (linfo.count == 0) {
        continue;
      }
      if (m_ndim == 3) {
        jump = _copy3d_x(outptr, src, linfo, xs, xe, ts, te, fromsrc);
        outptr += jump;
      } else {
        jump = _copy4d_xo(outptr, src, linfo, xs, xe, os, oe, ts, te, fromsrc);
        outptr += jump;
      }
    }
  }

  uint64_t to_remove = maxsize - (outptr - rw_mmap.data());
  rw_mmap.unmap();
  if (to_remove > 0) {
    truncate_file(outname, to_remove);
  }
}

inline void SegyND::cut(const std::string &segy_name,
                        const std::vector<int> &ranges, bool is2d,
                        const std::string &textual) {
  _create_from_segy(segy_name, nullptr, ranges, is2d, textual, true);
}

inline void SegyND::create_by_sharing_header(const std::string &segy_name,
                                             const float *src,
                                             const std::vector<int> &ranges,
                                             bool is2d,
                                             const std::string &textual) {
  _create_from_segy(segy_name, src, ranges, is2d, textual, false);
}

inline void SegyND::create_by_sharing_header(const std::string &segy_name,
                                             const std::string &src_name,
                                             const std::vector<int> &shape,
                                             const std::vector<int> &ranges,
                                             bool is2d,
                                             const std::string &textual) {
  std::error_code error;
  mio::mmap_source float_file;
  float_file.map(src_name, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
  // TODO: check file size
  const float *src = reinterpret_cast<const float *>(float_file.data());
  create_by_sharing_header(segy_name, src, ranges, is2d, textual);
}

} // namespace segy
#endif