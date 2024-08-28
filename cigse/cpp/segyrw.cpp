/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

#include "segyrw.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace segy {

// modify header keys function, for cut, create_by_sharing_header
inline static void set_bkeyi2(char *bheader, int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(bheader + loc - 1) = swap_endian<int16_t>(val);
}
inline static void set_keyi2(char *theader, int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(theader + loc - 1) = swap_endian<int16_t>(val);
}
inline static void set_keyi4(char *theader, int loc, int32_t val) {
  *reinterpret_cast<int32_t *>(theader + loc - 1) = swap_endian<int32_t>(val);
}

void SegyRW::scan() {
  m_meta.start_time = keyi2(0, kTStartTimeField);
  m_meta.scalar = keyi2(0, kTScalarField);

  // steps
  int istep = m_keys.istep;
  int xstep = m_keys.xstep;
  int ostep = m_keys.ostep;

  // range
  int is = iline(0);
  int ie = iline(m_meta.ntrace - 1);
  int ni = (ie - is) / m_keys.istep + 1;
  m_iinfos.resize(ni, LineInfo(true));

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
      linfo.line = m_iinfos[ii - 1].line + istep;
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
    it = it < m_meta.ntrace ? it : m_meta.ntrace - 1;

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

    if (it < ntrace() &&
        (iline(it) != (iiline + istep) || iline(it - 1) != iiline)) {
      std::ostringstream oss;
      oss << "Error when scan this file. We except `iline(i) == (line+istep)` "
             "and `iline(i-1)==line` when the line breaking, however, we got "
             "line = "
          << iiline << ", istep = " << istep << ", iline(i) = " << iline(it)
          << ", iline(i-1) = " << iline(it - 1) << ", i = " << it
          << ". Maybe this file is unsorted." << ntrace();
      throw std::runtime_error(oss.str());
    }

    // assign lineInfo
    int xend = xline(it - 1);
    linfo.line = iiline;
    linfo.itstart = itstart;
    linfo.itend = it - 1;
    linfo.lstart = xs;
    linfo.lend = xend;
    linfo.count = it - itstart;

    // update global xline range
    gxstart = xstep > 0 ? MMIN(gxstart, xs) : MMAX(gxstart, xs);
    gxend = xstep > 0 ? MMAX(gxend, xend) : MMIN(gxend, xend);

    if (is4D) {
      assert((m_iinfos[ii].lend - xs) % xstep == 0);
      int nx = (xend - xs) / xstep + 1;
      linfo.set_xinfos(nx);
      int xt = itstart;
      int xtmax = it;

      for (size_t ix = 0; ix < nx; ix++) {
        LineInfo &xinfo = linfo.xinfos[ix];
        // when missing xline
        if (skipx > 0) {
          xinfo.line = linfo.xinfos[ix - 1].line + xstep;
          xinfo.count = 0;
          skipx--;
          continue;
        }

        int xxline = iline(xt);
        int xtstart = xt;
        int ostart = offset(xt);

        xt += jumpx;
        xt = xt < m_meta.ntrace ? xt : m_meta.ntrace - 1;

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

        if (xt < ntrace() &&
            (xline(xt) != (xxline + xstep) || xline(xt - 1) != xxline)) {
          std::ostringstream oss;
          oss << "Error when scan this file. We except `xline(i) == "
                 "(line+xstep)` and `xline(i-1)==line` when the xline "
                 "breaking, however, we got line = "
              << xxline << ", xstep = " << xstep << ", iline(i) = " << xline(xt)
              << ", xline(i-1) = " << xline(xt - 1)
              << ". Maybe this file is unsorted.";
          throw std::runtime_error(oss.str());
        }

        // assign xline info
        int oend = offset(xt - 1);
        xinfo.line = xxline;
        xinfo.itstart = xtstart;
        xinfo.itend = xt - 1;
        xinfo.lstart = ostart;
        xinfo.lend = oend;
        xinfo.count = xt - xtstart;

        // update global offset range
        gostart = ostep > 0 ? MMIN(gostart, ostart) : MMAX(gostart, ostart);
        goend = ostep > 0 ? MMAX(goend, oend) : MMIN(goend, oend);

        // missing xline
        if (xline(xt) != (xxline + xstep)) {
          skipx = xline(xt) - (xxline + xstep);
          if (skipx % xstep != 0) {
            std::ostringstream oss;
            oss << "Error when scan this file. In this xline, there are "
                   "missing some offset lines, i.e., xline(i) != (line+xstep). "
                   "But we got xline(i) = "
                << xline(xt) << ", line = " << xxline << ", xstep = " << xstep
                << ". xline(i) - (line+xstep) = " << skipx
                << ", and skipx \% xstep != 0. Maybe the xstep is wrong.";
            throw std::runtime_error(oss.str());
          }
          skipx /= xstep;
        }
      }
    }

    // missing line
    if (iline(it) != (iiline + istep)) {
      skipi = xline(it) - (iiline + istep);
      if (skipi % istep != 0) {
        std::ostringstream oss;
        oss << "Error when scan this file. In this line, there are "
               "missing some xlines, i.e., iline(i) != (line+istep). "
               "But we got iline(i) = "
            << iline(it) << ", line = " << iiline << ", istep = " << istep
            << ". iline(i) - (line+istep) = " << skipi
            << ", and skipi \% istep != 0. Maybe the istep is wrong.";
        throw std::runtime_error(oss.str());
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

  // if line or xline is not continouse, we record their idx for fast indexing
  for (auto linfo : m_iinfos) {
    if (!is4D) {
      if (linfo.count == 0) {
        continue;
      }
      // continuous line
      if (linfo.count == ((linfo.lend - linfo.lstart) / m_keys.xstep + 1)) {
        continue;
      }

      // we set count = -1 for not continous line
      linfo.idx.resize(m_meta.nx, -1);
      for (size_t xt = linfo.itstart; xt < linfo.itend + 1; xt++) {
        linfo.idx[xl2ix(xline(xt))] = xt;
      }
      linfo.count = -1;
    } else {
      for (auto xinfo : linfo.xinfos) {
        if (xinfo.count == 0) {
          continue;
        }
        // continuous xline
        if (xinfo.count == ((xinfo.lend - xinfo.lstart) / m_keys.ostep + 1)) {
          continue;
        }

        // we set count = -1 for not continous xline
        xinfo.idx.resize(m_meta.no, -1);
        for (size_t ot = xinfo.itstart; ot < xinfo.itend + 1; ot++) {
          xinfo.idx[of2io(offset(ot))] = ot;
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

void SegyRW::read4d(float *dst, int is, int ie, int xs, int xe, int os, int oe,
                    int ts, int te) {
  if (!isScan || m_ndim != 4) {
    throw std::runtime_error("Not scan or is not a 4d pretrack SEG-Y");
  }

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

void SegyRW::read3d(float *dst, int is, int ie, int xs, int xe, int ts,
                    int te) {
  if (!isScan || m_ndim != 3) {
    throw std::runtime_error("Not scan or is not a 3d SEG-Y");
  }

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
    _read_inner(dstiline, linfo, xs, xe, ts, te);
  }
}

void SegyRW::read(float *dst) {
  if (m_ndim == 2) {
    collect(dst, 0, m_meta.ntrace, 0, m_meta.nt);
  } else if (m_ndim == 3) {
    read3d(dst, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.nt);
  } else if (m_ndim == 4) {
    read4d(dst, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.no, 0, m_meta.nt);
  }
}

void SegyRW::write(const float *data) {
  if (m_ndim == 2) {
    write_traces(data, 0, m_meta.ntrace, 0, m_meta.nt);
  } else if (m_ndim == 3) {
    write3d(data, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.nt);
  } else if (m_ndim == 4) {
    write4d(data, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.no, 0, m_meta.nt);
  }
}

void SegyRW::write3d(const float *data, int is, int ie, int xs, int xe, int ts,
                     int te) {
  if (!isScan || m_ndim != 3) {
    throw std::runtime_error("Not scan or is not a 3d SEG-Y");
  }

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
    const float *srcline = data + (ii - is) * sizeXT;
    _write_inner(srcline, linfo, xs, xe, ts, te);
  }
}

void SegyRW::write4d(const float *data, int is, int ie, int xs, int xe, int os,
                     int oe, int ts, int te) {
  if (!isScan || m_ndim != 4) {
    throw std::runtime_error("Not scan or is not a 4d pretrack SEG-Y");
  }

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
    const float *srcline = data + (ii - is) * sizeXOT;
    _write4d_xo(srcline, linfo, xs, xe, os, oe, ts, te);
  }
}

void SegyRW::tofile(const std::string &binary_out_name, bool is2d) {
  if (!is2d && m_ndim == 2) {
    throw std::runtime_error("ndim == 2, maybe you need set is2d=true");
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
    throw std::runtime_error("mmap fail in 'rw' mode: " + binary_out_name);
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
void SegyRW::_read_inner(float *dst, LineInfo &linfo, int ks, int ke, int ts,
                         int te) {
  float *odst = dst;
  int nt = te - ts;
  int nk = ke - ks;
  uint64_t sizeKT = nt * nk;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, ks, ke)) {
    std::fill(dst, dst + sizeKT, m_meta.fillNoValue);
    return;
  }

  /* line is not continuous */
  else if (linfo.count == -1) {
    for (size_t ik = ks; ik < ke; ik++) {
      int it = linfo.idx[ik];
      if (it == -1) {
        std::fill(dst, dst + nt, m_meta.fillNoValue);
      } else {
        m_readfunc(dst, trDataStart(it, ts), nt);
      }
      dst += nt;
    }
    return;
  }

  else {
    std::array<int32_t, 4> idx;
    find_idx(idx, linfo, ks, ke);

    if (idx[0] > 0) {
      std::fill(dst, dst + idx[0] * nt, m_meta.fillNoValue);
      dst += idx[0] * nt;
    }
    collect(dst, idx[2], idx[3], ts, te);
    dst += (idx[3] - idx[2]) * nt;

    if (idx[1] > 0) {
      std::fill(dst, dst + idx[1] * nt, m_meta.fillNoValue);
      dst += idx[1] * nt; // TODO: remove
    }
    assert((dst - odst) == sizeKT);
  }
}

void SegyRW::_read4d_xo(float *dst, LineInfo &linfo, int xs, int xe, int os,
                        int oe, int ts, int te) {
  float *odst = dst;
  int nt = te - ts;
  int no = oe - os;
  int nx = xe - xs;
  uint64_t sizeOT = nt * no;
  uint64_t sizeXOT = nx * sizeOT;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, xs, xe)) {
    std::fill(dst, dst + sizeXOT, m_meta.fillNoValue);
    return;
  }

  std::array<int32_t, 4> idx;
  find_idx(idx, linfo, xs, xe);

  // fill the left
  if (idx[0] > 0) {
    std::fill(dst, dst + idx[0] * sizeOT, m_meta.fillNoValue);
    dst += idx[0] * sizeOT;
  }

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = idx[2]; ix < idx[3]; ix++) {
    _read_inner(dst, linfo.xinfos[ix], os, oe, ts, te);
    dst += sizeOT;
  }

  // fill the right
  if (idx[1] > 0) {
    std::fill(dst, dst + idx[1] * sizeOT, m_meta.fillNoValue);
    dst += idx[1] * sizeOT;
  }
  assert((dst - odst) == sizeXOT);
  // assert
}

void SegyRW::_write_inner(const float *src, LineInfo &linfo, int ks, int ke,
                          int ts, int te) {
  int nt = te - ts;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, ks, ke)) {
    return;
  }

  /* line is not continuous */
  else if (linfo.count == -1) {
    for (size_t ik = ks; ik < ke; ik++) {
      int it = linfo.idx[ik];
      if (it >= 0) {
        m_wfunc(twDataStart(it, ts), src + (ik - ks) * nt, nt);
      }
    }
    return;
  }

  else {
    std::array<int32_t, 4> idx;
    find_idx(idx, linfo, ks, ke);

    write_traces(src + idx[0] * nt, idx[2], idx[3], ts, te);
  }
}

void SegyRW::_write4d_xo(const float *src, LineInfo &linfo, int xs, int xe,
                         int os, int oe, int ts, int te) {
  int nt = te - ts;
  int no = oe - os;
  uint64_t sizeOT = nt * no;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, xs, xe)) {
    return;
  }

  std::array<int32_t, 4> idx;
  find_idx(idx, linfo, xs, xe);

  const float *srcf = src + idx[0] * sizeOT;

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = idx[2]; ix < idx[3]; ix++) {
    _write_inner(srcf + (ix - idx[2]) * sizeOT, linfo.xinfos[ix], os, oe, ts,
                 te);
  }
}

uint64_t SegyRW::_copy_inner(char *dst, const float *src, LineInfo &linfo,
                             int ks, int ke, int ts, int te, bool fromsrc) {
  char *odst = dst;
  int nt = te - ts;
  bool tchanged = nt == m_meta.nt ? false : true;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, ks, ke)) {
    return 0;
  }

  /* line is not continuous */
  else if (linfo.count == -1) {
    for (size_t ik = ks; ik < ke; ik++) {
      int it = linfo.idx[ik];
      if (it >= 0) {
        // copy trace header
        memcpy(dst, trheader(it), kTraceHeaderSize);
        if (tchanged) {
          if (ts > 0) { // TODO: ts need *1000 ?
            segy::set_keyi2(dst, kTStartTimeField, ts * m_meta.dt / 1000);
          }
          segy::set_keyi2(dst, kTSampleCountField, nt);
        }
        dst += kTraceHeaderSize;

        // copy trace
        if (fromsrc) {
          // copy from src, i.e., cut
          memcpy(dst, trDataStart(it, ts), nt * m_meta.esize);
        } else {
          // copy from data, i.e., create_by_sharing_header
          m_wfunc(dst, src + (ik - ks) * nt, nt);
        }
        dst += nt * m_meta.esize;
      }
    }
    return dst - odst;
  }

  else {
    std::array<int32_t, 4> idx;
    find_idx(idx, linfo, ks, ke);
    const float *srcf = src + idx[0] * nt;

    for (size_t tx = idx[2]; tx < idx[3]; tx++) {
      // copy trace header
      memcpy(dst, trheader(tx), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) { // TODO: ts need *1000 ?
          segy::set_keyi2(dst, kTStartTimeField, ts * m_meta.dt / 1000);
        }
        segy::set_keyi2(dst, kTSampleCountField, nt);
      }
      dst += kTraceHeaderSize;

      // copy data
      if (fromsrc) {
        // copy from src, i.e., cut
        memcpy(dst, trDataStart(tx, ts), nt * m_meta.esize);
      } else {
        // copy from data, i.e., create_by_sharing_header
        m_wfunc(dst, srcf + (tx - idx[2]) * nt, nt);
      }
      dst += nt * m_meta.esize;
    }
    // assert
    return dst - odst;
  }
}

uint64_t SegyRW::_copy4d_xo(char *dst, const float *src, LineInfo &linfo,
                            int xs, int xe, int os, int oe, int ts, int te,
                            bool fromsrc) {
  char *odst = dst;
  int nt = te - ts;
  int no = oe - os;
  int nx = xe - xs;
  uint64_t sizeOT = nt * no;
  uint64_t sizeXOT = nx * sizeOT;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, xs, xe)) {
    return 0;
  }

  uint64_t jump = 0;

  std::array<int32_t, 4> idx;
  find_idx(idx, linfo, xs, xe);

  const float *srcf = src + idx[0] * sizeXOT;

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = idx[2]; ix < idx[3]; ix++) {
    jump = _copy_inner(dst, srcf + (ix - idx[2]) * sizeXOT, linfo.xinfos[ix],
                       os, oe, ts, te, fromsrc);
    dst += jump;
  }

  return dst - odst;
  // assert
}

void SegyRW::scanBinaryHeader() {
  int dformat = bkeyi2(kBSampleFormatField);
  auto it = kElementSize.find(dformat);
  if (it != kElementSize.end()) {
    m_meta.esize = it->second;
  } else {
    throw std::runtime_error("Unknown data format: " + std::to_string(dformat));
  }
  m_meta.dformat = dformat;

  m_meta.nt = bkeyi2(kBSampleCountField);
  m_meta.tracesize = kTraceHeaderSize + m_meta.nt * m_meta.esize;
  m_meta.dt = bkeyi2(kBSampleIntervalField);
  m_meta.ntrace = (m_sink.size() - kTraceHeaderStart) / m_meta.tracesize;

  m_meta.trace_sorting_code = bkeyi2(kBTraceSortingCodeField);
  setRWFunc(m_meta.dformat);
}

void SegyRW::_create_from_segy(const std::string &outname, const float *src,
                               const std::vector<int> &ranges, bool is2d,
                               const std::string &textual, bool fromsrc,
                               uint64_t check_size) {
  int rsize = ranges.size();
  if (rsize < 4) {
    throw std::runtime_error(
        "Error ranges size, ranges.size() must >= 4, but got " +
        std::to_string(rsize));
  }
  int ts = ranges[rsize - 2], te = ranges[rsize - 1];
  if (ts < 0 || te > m_meta.nt) {
    throw std::runtime_error("ts and te index out of range");
  }

  int tstart, tend, is, ie, xs, xe, os, oe;
  uint64_t maxsize = 0, tracesize = 0;
  if (is2d || m_ndim == 2) {
    if (ranges.size() != 4) {
      throw std::runtime_error("Error ranges. you set (is2d || ndim==2), but "
                               "ranges.size() != 4, ranges.size() = " +
                               std::to_string(ranges.size()));
    }
    tstart = ranges[0];
    tend = ranges[1];
    if (check_size > 0 and
        check_size != (uint64_t)(tend - tstart) * (te - ts) * sizeof(float)) {
      throw std::runtime_error("file size don't match the ranges");
    }
    tracesize = kTraceHeaderSize + (uint64_t)(te - ts) * m_meta.esize;
    maxsize = kTraceHeaderStart + tracesize * (tend - tstart);

  } else {
    is = ranges[0];
    ie = ranges[1];
    xs = ranges[2];
    xe = ranges[3];

    if (m_ndim == 3) {
      if (ranges.size() != 6) {
        throw std::runtime_error("Error ranges. ndim == 3, but ranges.size() "
                                 "!= 6, ranges.size() = " +
                                 std::to_string(ranges.size()));
      }
      if (check_size > 0 and check_size != (uint64_t)(ie - is) * (xe - xs) *
                                               (te - ts) * sizeof(float)) {
        throw std::runtime_error("file size don't match the ranges");
      }
      tracesize = kTraceHeaderSize + (uint64_t)(te - ts) * m_meta.esize;
      maxsize = kTraceHeaderStart + tracesize * (ie - is) * (xe - xs);
    } else {
      if (ranges.size() != 8) {
        throw std::runtime_error("Error ranges. ndim == 4, but ranges.size() "
                                 "!= 8, ranges.size() = " +
                                 std::to_string(ranges.size()));
      }
      os = ranges[4];
      oe = ranges[5];
      if (check_size > 0 and check_size != (uint64_t)(ie - is) * (xe - xs) *
                                               (oe - os) * (te - ts) *
                                               sizeof(float)) {
        throw std::runtime_error("file size don't match the ranges");
      }
      tracesize = kTraceHeaderSize + (uint64_t)(te - ts) * m_meta.esize;
      maxsize =
          kTraceHeaderStart + tracesize * (ie - is) * (xe - xs) * (oe - os);
    }
  }

  // set dimension
  int nt = te - ts;

  bool tchanged = nt == m_meta.nt ? false : true;

  // We first create a large file (assume there is no missing)
  // to store all data, and then mmap it to
  create_file(outname, maxsize);

  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(outname, error);
  if (error) {
    throw std::runtime_error("mmap fail in 'rw' mode: " + outname);
  }

  char *outptr = rw_mmap.data();

  // copy textual header
  if (textual.size() == 0) {
    memcpy(outptr, m_data_ptr, kTextualHeaderSize);
  } else if (textual.size() == kTextualHeaderSize) {
    memcpy(outptr, textual.data(), kTextualHeaderSize);
  } else {
    throw std::runtime_error("textual header size not match, textual.size() "
                             "must be 0 or 3200, but got " +
                             std::to_string(textual.size()));
  }
  outptr += kTextualHeaderSize;

  // copy binary header
  memcpy(outptr, brheader(), kBinaryHeaderSize);
  if (tchanged) {
    segy::set_bkeyi2(outptr, kBSampleCountField, te - ts);
  }
  outptr += kBinaryHeaderSize;

  // copy trace
  if (is2d || m_ndim == 2) {
    for (size_t it = tstart; it < tend; it++) {
      CHECK_SIGNALS();
      memcpy(outptr, trheader(it), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) { // TODO: ts need *1000 ?
          segy::set_keyi2(outptr, kTStartTimeField, ts);
        }
        segy::set_keyi2(outptr, kTSampleCountField, nt);
      }
      outptr += kTraceHeaderSize;
      memcpy(outptr, trDataStart(it, ts), nt * m_meta.esize);
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
        jump = _copy_inner(outptr, src, linfo, xs, xe, ts, te, fromsrc);
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

void SegyRW::cut(const std::string &segy_name, const std::vector<int> &ranges,
                 bool is2d, const std::string &textual) {
  _create_from_segy(segy_name, nullptr, ranges, is2d, textual, true);
}

void SegyRW::create_by_sharing_header(const std::string &segy_name,
                                      const float *src,
                                      const std::vector<int> &shape,
                                      const std::vector<int> &start, bool is2d,
                                      const std::string &textual) {
  int nd = is2d ? 2 : m_ndim;

  if (shape.size() != start.size() || shape.size() != nd) {
    throw std::runtime_error("shape's size must be equal to the start's size, "
                             "and shape's size must be equal to ndim");
  }
  std::vector<int> ranges(nd * 2);
  for (int i = 0; i < nd; ++i) {
    ranges[i * 2] = start[i];
    ranges[i * 2 + 1] = start[i] + shape[i];
  }
  _create_from_segy(segy_name, src, ranges, is2d, textual, false);
}

void SegyRW::create_by_sharing_header(const std::string &segy_name,
                                      const std::string &src_name,
                                      const std::vector<int> &shape,
                                      const std::vector<int> &start, bool is2d,
                                      const std::string &textual) {
  std::error_code error;
  mio::mmap_source float_file;
  float_file.map(src_name, error);
  if (error) {
    throw std::runtime_error("Cannot mmap binary file in 'r' mode: " +
                             src_name);
  }

  int nd = is2d ? 2 : m_ndim;
  if (shape.size() != start.size() || shape.size() != nd) {
    throw std::runtime_error("shape's size must be equal to the start's size, "
                             "and shape's size must be equal to ndim");
  }
  std::vector<int> ranges(nd * 2);
  for (int i = 0; i < nd; ++i) {
    ranges[i * 2] = start[i];
    ranges[i * 2 + 1] = start[i] + shape[i];
  }

  const float *src = reinterpret_cast<const float *>(float_file.data());
  _create_from_segy(segy_name, src, ranges, is2d, textual, false,
                    float_file.size());
}

inline static void modify_keys(char *dst, const int32_t *keys, int keysize) {
  set_keyi4(dst, 5, keys[0]);
  set_keyi4(dst, 9, keys[0]);
  set_keyi4(dst, 189, keys[0]);
  set_keyi4(dst, 17, keys[1]);
  set_keyi4(dst, 21, keys[1]);
  set_keyi4(dst, 193, keys[1]);
  if (keysize == 4) {
    set_keyi4(dst, 73, keys[2]);
    set_keyi4(dst, 77, keys[3]);
    set_keyi4(dst, 181, keys[2]);
    set_keyi4(dst, 185, keys[3]);
  } else if (keysize == 5) {
    set_keyi4(dst, 37, keys[2]);
    set_keyi4(dst, 73, keys[3]);
    set_keyi4(dst, 77, keys[4]);
    set_keyi4(dst, 181, keys[3]);
    set_keyi4(dst, 185, keys[4]);
  }
}

void create_segy(const std::string &segyname, const float *src,
                 const int32_t *keys, const std::vector<int> &shape,
                 const std::string &textual, const uchar *bheader,
                 const uchar *theader, int keysize) {
  int ndim = shape.size();
  int ni, nx, no;
  uint64_t nt;
  uint64_t ntrace;
  if (ndim == 2) {
    ntrace = shape[0];
    nt = shape[1];
  } else if (ndim == 3) {
    ni = shape[0];
    nx = shape[1];
    nt = shape[2];
    ntrace = ni * nx;
  } else if (ndim == 4) {
    ni = shape[0];
    nx = shape[1];
    no = shape[2];
    nt = shape[3];
    ntrace = ni * nx * no;
  } else {
    throw std::runtime_error("shape's size must be 2, 3 or 4, but got " +
                             std::to_string(ndim));
  }
  uint64_t needsize = kTraceHeaderStart + ntrace * (nt * 4 + kTraceHeaderSize);
  create_file(segyname, needsize);

  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(segyname, error);
  if (error) {
    throw std::runtime_error("mmap fail in 'rw' mode: " + segyname);
  }
  char *dst = rw_mmap.data();
  int dformat = swap_endian<int16_t>(bheader + kBSampleFormatField);
  int esize = 0;

  auto it = kElementSize.find(dformat);
  if (it != kElementSize.end()) {
    esize = it->second;
  } else {
    throw std::runtime_error("Unknown data format: " + std::to_string(dformat));
  }

  WriteFunc wfunc;
  setWFunc(wfunc, dformat);

  memcpy(dst, textual.data(), kTextualHeaderSize);
  dst += kTextualHeaderSize;

  memcpy(dst, bheader, kBinaryHeaderSize);
  dst += kBinaryHeaderSize;

  for (size_t i = 0; i < ntrace; i++) {
    CHECK_SIGNALS();
    memcpy(dst, theader, kTraceHeaderSize);
    modify_keys(dst, keys + i * (uint64_t)keysize, keysize);
    dst += kTraceHeaderSize;
    wfunc(dst, src + i * nt, nt);
    dst += nt * esize;
  }

  assert(dst - rw_mmap.data() == needsize);
  rw_mmap.unmap();
}

} // namespace segy
