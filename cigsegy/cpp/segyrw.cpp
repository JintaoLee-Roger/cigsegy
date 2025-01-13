/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

#include "segyrw.h"
#include "utils.hpp"
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
// #include <chrono>

// static int tkscount = 0;
namespace segy {

bool g_show_progress = true;

std::function<void()> g_check_signals_callback = []() {};
std::function<void(int, int)> g_progress_callback = default_progress_callback;

// modify header keys function, for cut, create_by_sharing_header
inline static void set_bkeyi2(char *bheader, size_t loc, int16_t val) {
  *reinterpret_cast<int16_t *>(bheader + loc - 1) = swap_endian<int16_t>(val);
}
inline static void set_keyi2(char *theader, size_t loc, int16_t val) {
  *reinterpret_cast<int16_t *>(theader + loc - 1) = swap_endian<int16_t>(val);
}
inline static void set_keyi4(char *theader, size_t loc, int32_t val) {
  *reinterpret_cast<int32_t *>(theader + loc - 1) = swap_endian<int32_t>(val);
}

void SegyRW::scan() {
  int start_time = keyi2(0, kTStartTimeField);
  if (start_time < 0) {
    int start_time2 = keyi2(0, kTDelayTimeField);
    if (start_time2 < 0) {
      m_meta.start_time = 0;
    } else {
      m_meta.start_time = start_time2;
    }
  } else {
    m_meta.start_time = start_time;
  }
  m_meta.scalar = keyi2(0, kTScalarField);

  // steps
  int istep = m_keys.istep;
  int xstep = m_keys.xstep;
  int ostep = m_keys.ostep;

  // range
  int is = iline(0);
  int ie = iline(m_meta.ntrace - 1);
  int ni = (ie - is) / m_keys.istep + 1;
  if (ni < 0 || ni > kMaxSizeOneDimemsion) {
    std::ostringstream oss;
    oss << "Error when scan this file. The size of inline is error. ni = " << ni
        << ", max size = " << kMaxSizeOneDimemsion
        << ". If ni < 0, maybe istep is error. We got iline(0) = " << is
        << ", iline(ntrace-1) = " << ie << ", istep = " << istep << ".";
    throw std::runtime_error(oss.str());
  }

  m_iinfos.resize(ni, LineInfo(true));

  // global xline range and offset range
  int gxstart = xline(0);
  int gxend = xline(m_meta.ntrace - 1);
  int gostart = offset(0);
  int goend = offset(m_meta.ntrace - 1);

  // is a 4D pretack SEG-Y?
  bool is4D = true;
  if (m_ndim == 2) {
    is4D = isPreStack(); // NOTE: isPreStack() is very simple, maybe wrong
  } else if (m_ndim == 3) {
    is4D = false;
  }

  // use skipi and skipx to process when missing line and xline (in 4D)
  int skipi = 0;
  int skipx = 0;

  // trace index for iline
  size_t it = 0;
  size_t jumpl = m_meta.ntrace / ni;
  size_t jumpx = 1;
  // int step = cal_progress_steps(ni, show_progress_, 50);
  for (size_t ii = 0; ii < ni; ++ii) {
    g_check_signals_callback();
    // if (step && ii % step == 0) {
    //   g_progress_callback(ii, ni);
    // }
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
    size_t itstart = it;
    int xs = xline(it);

    // jump
    it += jumpl;
    it = it < m_meta.ntrace ? it : m_meta.ntrace - 1;

    // jump too small
    if (iline(it) == iiline) {
      while (it < m_meta.ntrace && iline(it) == iiline) {
        jumpl++;
        it++;
      }
    }
    // jump too large
    else if (iline(it - 1) != iiline) {
      while (it >= itstart && iline(it - 1) != iiline) {
        it--;
        jumpl--;
      }
    }

    // missing line
    if (it < ntrace() &&
        (iline(it) != (iiline + istep) || iline(it - 1) != iiline)) {
      if (iline(it - 1) == iiline) {
        skipi = iline(it) - (iiline + istep);
        if (skipi % istep != 0) {
          std::ostringstream oss;
          oss << "Error when scan this file. In this line, there are "
                 "missing some xlines, i.e., iline(i) != (line+istep). "
                 "But we got iline(i) = "
              << iline(it) << ", line = " << iiline << ", istep = " << istep
              << ". iline(i) - (line+istep) = " << skipi
              << ", and skipi % istep != 0. Maybe the istep is wrong.";
          throw std::runtime_error(oss.str());
        }
        skipi /= istep;
        if (skipi < 0) {
          std::ostringstream oss;
          oss << "Error when scan this file. skipi < 0. Maybe the istep is "
                 "wrong. "
              << "iline(i-1) = " << iline(it - 1)
              << ", iline(i) = " << iline(it) << ", line = " << iiline
              << ", istep = " << istep;
          throw std::runtime_error(oss.str());
        }
      } else {
        std::ostringstream oss;
        oss << "Error when scan this file. We except `iline(i) == "
               "(line+istep)` "
               "and `iline(i-1)==line` when the line breaking, however, we got "
               "line = "
            << iiline << ", istep = " << istep << ", iline(i) = " << iline(it)
            << ", iline(i-1) = " << iline(it - 1) << ", i = " << it
            << ". Maybe this file is unsorted." << ntrace();
        throw std::runtime_error(oss.str());
      }
    }

    // assign lineInfo
    int xend = xline(it - 1);
    linfo.line = iiline;
    linfo.itstart = itstart;
    linfo.itend = it - 1;
    linfo.lstart = xs;
    linfo.lend = xend;
    linfo.count = it - itstart;
    if (!is4D &&
        linfo.count != ((linfo.lend - linfo.lstart) / m_keys.xstep + 1)) {
      linfo.count = kInvalid;
    }

    // update global xline range
    gxstart = xstep > 0 ? MMIN(gxstart, xs) : MMAX(gxstart, xs);
    gxend = xstep > 0 ? MMAX(gxend, xend) : MMIN(gxend, xend);

    if (is4D) {
      assert((m_iinfos[ii].lend - xs) % xstep == 0);
      int nx = (xend - xs) / xstep + 1;
      if (nx < 0 || nx > kMaxSizeOneDimemsion) {
        std::ostringstream oss;
        oss << "Error when scan this file. The size of xline is error. nx = "
            << nx << ", max size = " << kMaxSizeOneDimemsion
            << ". If nx < 0, maybe xstep is error. xend = " << xend
            << ", xstart = " << xs << ", xstep = " << xstep;
        throw std::runtime_error(oss.str());
      }

      linfo.set_xinfos(nx);
      size_t xt = itstart;
      size_t xtmax = it;

      for (size_t ix = 0; ix < nx; ix++) {
        LineInfo &xinfo = linfo.xinfos[ix];
        // when missing xline
        if (skipx > 0) {
          xinfo.line = linfo.xinfos[ix - 1].line + xstep;
          xinfo.count = 0;
          skipx--;
          continue;
        }

        int xxline = xline(xt);
        size_t xtstart = xt;
        int ostart = offset(xt);

        xt += jumpx;
        xt = xt < xtmax ? xt : xtmax - 1;

        // jump too small
        if (xline(xt) == xxline) {
          while (xt < xtmax && xline(xt) == xxline) {
            jumpx++;
            xt++;
          }
        }
        // jump too large
        else if (xline(xt - 1) != xxline) {
          while (xt >= itstart && xline(xt - 1) != xxline) {
            xt--;
            jumpx--;
          }
        }

        if (xt < xtmax &&
            (xline(xt) != (xxline + xstep) || xline(xt - 1) != xxline)) {
          if (xline(xt - 1) == xxline) {
            skipx = xline(xt) - (xxline + xstep);
            if (skipx % xstep != 0) {
              std::ostringstream oss;
              oss << "Error when scan this file. In this xline, there are "
                     "missing some offset lines, i.e., xline(i) != "
                     "(line+xstep). "
                     "But we got xline(i) = "
                  << xline(xt) << ", line = " << xxline << ", xstep = " << xstep
                  << ". xline(i) - (line+xstep) = " << skipx
                  << ", and skipx % xstep != 0. Maybe the xstep is wrong.";
              throw std::runtime_error(oss.str());
            }
            skipx /= xstep;
            if (skipx < 0) {
              std::ostringstream oss;
              oss << "Error when scan this file. skipx < 0. Maybe the istep is "
                     "wrong. "
                  << "xline(i-1) = " << xline(xt - 1)
                  << ", xline(i) = " << xline(xt) << ", xxline = " << xxline
                  << ", xstep = " << xstep;
              throw std::runtime_error(oss.str());
            }
          } else {
            std::ostringstream oss;
            oss << "Error when scan this file. We except `xline(i) == "
                   "(line+xstep)` and `xline(i-1)==line` when the xline "
                   "breaking, however, we got line = "
                << xxline << ", xstep = " << xstep
                << ", iline(i) = " << xline(xt)
                << ", xline(i-1) = " << xline(xt - 1)
                << ". Maybe this file is unsorted.";
            throw std::runtime_error(oss.str());
          }
        }

        // assign xline info
        int oend = offset(xt - 1);
        xinfo.line = xxline;
        xinfo.itstart = xtstart;
        xinfo.itend = xt - 1;
        xinfo.lstart = ostart;
        xinfo.lend = oend;
        xinfo.count = xt - xtstart;
        if (xinfo.count != ((xinfo.lend - xinfo.lstart) / m_keys.ostep + 1)) {
          xinfo.count = kInvalid;
        }

        // update global offset range
        gostart = ostep > 0 ? MMIN(gostart, ostart) : MMAX(gostart, ostart);
        goend = ostep > 0 ? MMAX(goend, oend) : MMIN(goend, oend);
      }
    }
  }
  // if (step) {
  //   g_progress_callback(-1, ni);
  // }

  // Post process, assign m_meta
  m_meta.start_iline = is;
  m_meta.end_iline = ie;
  m_meta.ni = ni;
  m_meta.start_xline = gxstart;
  m_meta.end_xline = gxend;
  int nx = (gxend - gxstart) / xstep + 1;
  if (nx < 0 || nx > kMaxSizeOneDimemsion) {
    std::ostringstream oss;
    oss << "Error when scan this file. The size of xline is error. nx = " << nx
        << ", max size = " << kMaxSizeOneDimemsion
        << ". If nx < 0, maybe xstep is error. gxend = " << gxend
        << ", gxstart = " << gxstart << ", xstep = " << xstep;
    throw std::runtime_error(oss.str());
  }
  m_meta.nx = nx;

  if (is4D) {
    m_meta.start_offset = gostart;
    m_meta.end_offset = goend;
    int no = (goend - gostart) / ostep + 1;
    if (no < 0 || no > kMaxSizeOneDimemsion) {
      std::ostringstream oss;
      oss << "Error when scan this file. The size of offset is error. no = "
          << no << ", max size = " << kMaxSizeOneDimemsion
          << ". If no < 0, maybe ostep is error. goend = " << goend
          << ", gostart = " << gostart << ", ostep = " << ostep;
      throw std::runtime_error(oss.str());
    }
    m_meta.no = no;
  }

  // if line or xline is not continouse, we record their idx for fast indexing
  size_t rcount = 0; // We only read 10 lines? OPTIMIZE: How many lines?
  for (auto& linfo : m_iinfos) {
    g_check_signals_callback();
    if (rcount > 10) {
      break;
    }
    if (!is4D) {
      if (!(linfo.count == kInvalid)) {
        continue;
      }

      // we set count = kInvalid for not continous line
      linfo.idx.resize(m_meta.nx, kInvalid);
      for (size_t xt = linfo.itstart; xt < linfo.itend + 1; xt++) {
        linfo.idx[itx2ix(xt)] = xt;
      }
      linfo.count = 0;
      rcount++;
    } else {
      for (auto& xinfo : linfo.xinfos) {
        if (!(xinfo.count == kInvalid)) {
          continue;
        }

        // we set count = -1 for not continous xline
        xinfo.idx.resize(m_meta.no, kInvalid);
        for (size_t ot = xinfo.itstart; ot < xinfo.itend + 1; ot++) {
          xinfo.idx[itx2io(ot)] = ot;
        }
        xinfo.count = 0;
        rcount++;
      }
    }
  }

  { /********** calculate interval **********/

    // cal x, y interval
    // find a line contained more than 5 traces
    size_t s = m_meta.ni / 8;
    while (m_iinfos[s].count < 6 && s < m_meta.ni) {
      s++;
    }

    istep = std::abs(istep);
    xstep = std::abs(xstep);
    // Y_interval = ((x2-x1)^2 + (y2-y1)^2)^0.5 / abs(xline2 - xline1) * xstep
    size_t te = m_iinfos[s].itend;
    size_t ts = m_iinfos[s].itstart;
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

void SegyRW::read4d(float *dst, size_t is, size_t ie, size_t xs, size_t xe,
                    size_t os, size_t oe, size_t ts, size_t te) {
  if (!isScan || m_ndim != 4) {
    throw std::runtime_error("Not scan or is not a 4d pretrack SEG-Y");
  }

  size_t nx = xe - xs;
  size_t no = oe - os;
  size_t nt = te - ts;
  uint64_t sizeOT = no * nt;
  uint64_t sizeXOT = nx * sizeOT;

  int step = cal_progress_steps(ie - is, show_progress_, 50, 50);
  for (size_t ii = is; ii < ie; ii++) {
    g_check_signals_callback();
    if (step && (ii - is) % step == 0) {
      g_progress_callback(ii - is, ie - is);
    }

    LineInfo &linfo = m_iinfos[ii];
    float *dstiline = dst + (ii - is) * sizeXOT;
    _read4d_xo(dstiline, linfo, xs, xe, os, oe, ts, te);
  }
  if (step) {
    g_progress_callback(-1, ie - is);
  }
}

void SegyRW::read3d(float *dst, size_t is, size_t ie, size_t xs, size_t xe,
                    size_t ts, size_t te) {
  if (!isScan || m_ndim != 3) {
    throw std::runtime_error("Not scan or is not a 3d SEG-Y");
  }

  size_t nx = xe - xs;
  size_t nt = te - ts;
  uint64_t sizeXT = nx * nt;

  int step = cal_progress_steps(ie - is, show_progress_, 50, 50);
  for (size_t ii = is; ii < ie; ii++) {
    g_check_signals_callback();
    if (step && (ii - is) % step == 0) {
      g_progress_callback(ii - is, ie - is);
    }

    LineInfo &linfo = m_iinfos[ii];
    float *dstiline = dst + (ii - is) * sizeXT;
    _read_inner(dstiline, linfo, xs, xe, ts, te);
  }
  if (step) {
    g_progress_callback(-1, ie - is);
  }
}

// This function is for accelerated reading a time slice of a 3D SEG-Y.
// We set setpi and setpx to read a part of the data,
// and remove lots of if conditions.
void SegyRW::read_tslice(float *dst, size_t t, size_t stepi, size_t stepx) {
  if (!isScan || m_ndim != 3) {
    throw std::runtime_error("Not scan or is not a 3d SEG-Y");
  }
  if (t > m_meta.nt) {
    throw std::runtime_error("t > nt");
  }
  uint64_t sizeXT = (m_meta.nx + stepx - 1) / stepx;

  for (size_t ii = 0; ii < m_meta.ni; ii += stepi) {
    g_check_signals_callback();

    LineInfo &linfo = m_iinfos[ii];
    float *dstl = dst + ii / stepi * sizeXT;
    // no need to check NoOverlap
    if (linfo.count == 0 && linfo.idx.size() == 0) {
      std::fill(dstl, dstl + sizeXT, m_meta.fillNoValue);
      continue;
    }

    // line is not continuous
    else if (linfo.count == 0 && linfo.idx.size() > 0) {

      for (size_t ix = 0; ix < m_meta.nx; ix += stepx) {
        size_t it = linfo.idx[ix];
        if (it == kInvalid) {
          *dstl = m_meta.fillNoValue;
        } else {
          *dstl = m_readfuncone(trDataStart(it, t));
        }
        dstl += 1;
      }

    }

    // not continuous, but we don't set idx
    else if (linfo.count == kInvalid) {
      std::fill(dstl, dstl + sizeXT, m_meta.fillNoValue);
      // find start & end
      size_t its, ite;
      find_nearest_idx(linfo, 0, m_meta.nx, its, ite);
      for (size_t i = its; i < ite; i++) {
        size_t ix = linfo.isline ? itx2ix(i) : itx2io(i);
        if (ix % stepx == 0) {
          *(dstl + ix / stepx) = m_readfuncone(trDataStart(i, t));
        }
      }
      if (ite > its) {
        dstl += sizeXT;
      }
    }

    // normal line
    else if (linfo.count == m_meta.nx) {

      for (size_t ix = 0; ix < m_meta.nx; ix += stepx) {
        size_t it = linfo.itstart + ix;
        *dstl = m_readfuncone(trDataStart(it, t));
        dstl += 1;
      }
    }

    // continuous line
    else {
      float *odst = dstl;
      size_t left = (xl2ix(linfo.lstart) + stepx - 1) / stepx;
      if (left > 0) {
        std::fill(dstl, dstl + left, m_meta.fillNoValue);
        dstl += left;
      }
      size_t right = xl2ix(linfo.lend) + 1;
      for (size_t ix = left * stepx; ix < right; ix += stepx) {
        *dstl = m_readfuncone(trDataStart(linfo.itstart + ix, t));
        dstl += 1;
      }
      size_t fill = sizeXT - (dstl - odst);
      if (fill > 0) {
        std::fill(dstl, dstl + fill, m_meta.fillNoValue);
      }
    }
  }
}

void SegyRW::read(float *dst) {
  if (m_ndim == 2) {
    collect(dst, (size_t)0, m_meta.ntrace, 0, m_meta.nt);
  } else if (m_ndim == 3) {
    read3d(dst, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.nt);
  } else if (m_ndim == 4) {
    read4d(dst, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.no, 0, m_meta.nt);
  }
}

void SegyRW::write(const float *data) {
  check_write(m_w);
  if (m_ndim == 2) {
    write_traces(data, (size_t)0, m_meta.ntrace, 0, m_meta.nt);
  } else if (m_ndim == 3) {
    write3d(data, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.nt);
  } else if (m_ndim == 4) {
    write4d(data, 0, m_meta.ni, 0, m_meta.nx, 0, m_meta.no, 0, m_meta.nt);
  }
}

void SegyRW::write3d(const float *data, size_t is, size_t ie, size_t xs,
                     size_t xe, size_t ts, size_t te) {
  check_write(m_w);
  if (!isScan || m_ndim != 3) {
    throw std::runtime_error("Not scan or is not a 3d SEG-Y");
  }

  size_t nx = xe - xs;
  size_t nt = te - ts;
  uint64_t sizeXT = nx * nt;

  int step = cal_progress_steps(ie - is, show_progress_, 50, 50);
  for (size_t ii = is; ii < ie; ii++) {
    g_check_signals_callback();
    if (step && (ii - is) % step == 0) {
      g_progress_callback(ii - is, ie - is);
    }

    LineInfo &linfo = m_iinfos[ii];
    const float *srcline = data + (ii - is) * sizeXT;
    _write_inner(srcline, linfo, xs, xe, ts, te);
  }
  if (step) {
    g_progress_callback(-1, ie - is);
  }
}

void SegyRW::write4d(const float *data, size_t is, size_t ie, size_t xs,
                     size_t xe, size_t os, size_t oe, size_t ts, size_t te) {
  check_write(m_w);
  if (!isScan || m_ndim != 4) {
    throw std::runtime_error("Not scan or is not a 4d pretrack SEG-Y");
  }

  size_t nx = xe - xs;
  size_t no = oe - os;
  size_t nt = te - ts;
  uint64_t sizeOT = no * nt;
  uint64_t sizeXOT = nx * sizeOT;

  int step = cal_progress_steps(ie - is, show_progress_, 50, 50);
  for (size_t ii = is; ii < ie; ii++) {
    g_check_signals_callback();
    if (step && (ii - is) % step == 0) {
      g_progress_callback(ii - is, ie - is);
    }

    LineInfo &linfo = m_iinfos[ii];
    const float *srcline = data + (ii - is) * sizeXOT;
    _write4d_xo(srcline, linfo, xs, xe, os, oe, ts, te);
  }
  if (step) {
    g_progress_callback(-1, ie - is);
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
    collect(dst, (size_t)0, m_meta.ntrace, 0, m_meta.nt);
  }

  rw_mmap.unmap();
}

/*   priviate function for reading */
void SegyRW::_read_inner(float *dst, LineInfo &linfo, size_t ks, size_t ke,
                         size_t ts, size_t te) {
  size_t nt = te - ts;
  size_t nk = ke - ks;
  uint64_t sizeKT = nt * nk;

  /* no data to copy, just fill */
  if ((linfo.count == 0 && linfo.idx.size() == 0) || NoOverlap(linfo, ks, ke)) {
    std::fill(dst, dst + sizeKT, m_meta.fillNoValue);
    return;
  }

  /* line is not continuous */
  else if (linfo.count == 0 && linfo.idx.size() > 0) {
    for (size_t ik = ks; ik < ke; ik++) {
      size_t it = linfo.idx[ik];
      if (it == kInvalid) {
        std::fill(dst, dst + nt, m_meta.fillNoValue);
      } else {
        m_readfunc(dst, trDataStart(it, ts), nt);
      }
      dst += nt;
    }
    return;
  }

  /* line is not continuous, but we don't set idx */
  else if (linfo.count == kInvalid) {
    std::fill(dst, dst + sizeKT, m_meta.fillNoValue);
    // find start & end
    size_t its, ite;
    find_nearest_idx(linfo, ks, ke, its, ite);
    if (linfo.isline) {
      for (size_t i = its; i < ite; i++) {
        m_readfunc(dst + (itx2ix(i) - ks) * nt, trDataStart(i, ts), nt);
      }
    } else {
      for (size_t i = its; i < ite; i++) {
        m_readfunc(dst + (itx2io(i) - ks) * nt, trDataStart(i, ts), nt);
      }
    }
  }

  // continuous line
  else {
    // float *odst = dst;

    std::array<size_t, 4> idx;
    find_idx(idx, linfo, ks, ke);

    if (idx[0] > 0) {
      std::fill(dst, dst + idx[0] * nt, m_meta.fillNoValue);
      dst += idx[0] * nt;
    }
    collect(dst, idx[2], idx[3], ts, te);
    dst += (idx[3] - idx[2]) * nt;

    if (idx[1] > 0) {
      std::fill(dst, dst + idx[1] * nt, m_meta.fillNoValue);
      dst += idx[1] * nt;
    }
    // assert((dst - odst) == sizeKT);
  }
}

void SegyRW::_read4d_xo(float *dst, LineInfo &linfo, size_t xs, size_t xe,
                        size_t os, size_t oe, size_t ts, size_t te) {
  size_t nt = te - ts;
  size_t no = oe - os;
  size_t nx = xe - xs;
  uint64_t sizeOT = nt * no;
  uint64_t sizeXOT = nx * sizeOT;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, xs, xe)) {
    std::fill(dst, dst + sizeXOT, m_meta.fillNoValue);
    return;
  }

  // float *odst = dst;

  size_t start = xl2ix(linfo.lstart);
  size_t end = xl2ix(linfo.lend) + 1;

  // fill the left
  if (xs < start) {
    std::fill(dst, dst + (start - xs) * sizeOT, m_meta.fillNoValue);
    dst += (start - xs) * sizeOT;
    xs = start;
  }
  size_t xend = xe, tofill = 0;
  if (xe > end) {
    xend = end;
    tofill = xe - end;
  }

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = xs; ix < xend; ix++) {
    _read_inner(dst, linfo.xinfos[ix - start], os, oe, ts, te);
    dst += sizeOT;
  }

  // fill the right
  if (tofill) {
    std::fill(dst, dst + tofill * sizeOT, m_meta.fillNoValue);
    dst += tofill * sizeOT;
  }
}

void SegyRW::_write_inner(const float *src, LineInfo &linfo, size_t ks,
                          size_t ke, size_t ts, size_t te) {
  size_t nt = te - ts;

  /* no data to copy, just fill */
  if ((linfo.count == 0 && linfo.idx.size() == 0) || NoOverlap(linfo, ks, ke)) {
    return;
  }

  /* line is not continuous */
  else if (linfo.count == 0 && linfo.idx.size() > 0) {
    for (size_t ik = ks; ik < ke; ik++) {
      size_t it = linfo.idx[ik];
      if (it != kInvalid) {
        m_wfunc(twDataStart(it, ts), src + (ik - ks) * nt, nt);
      }
    }
    return;
  }

  /* line is not continuous, not set idx */
  else if (linfo.count == kInvalid) {
    // find start & end
    size_t its, ite;
    find_nearest_idx(linfo, ks, ke, its, ite);
    if (linfo.isline) {
      for (size_t it = its; it < ite; it++) {
        m_wfunc(twDataStart(it, ts), src + (itx2ix(it) - ks) * nt, nt);
      }
    } else {
      for (size_t it = its; it < ite; it++) {
        m_wfunc(twDataStart(it, ts), src + (itx2io(it) - ks) * nt, nt);
      }
    }
  }

  // continuous
  else {
    std::array<size_t, 4> idx;
    find_idx(idx, linfo, ks, ke);

    write_traces(src + idx[0] * nt, idx[2], idx[3], ts, te);
  }
}

void SegyRW::_write4d_xo(const float *src, LineInfo &linfo, size_t xs,
                         size_t xe, size_t os, size_t oe, size_t ts,
                         size_t te) {
  size_t nt = te - ts;
  size_t no = oe - os;
  uint64_t sizeOT = nt * no;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, xs, xe)) {
    return;
  }

  size_t start = xl2ix(linfo.lstart);
  size_t end = xl2ix(linfo.lend) + 1;
  size_t skip = 0;

  if (xs < start) {
    skip = start - xs;
    xs = start;
  }
  if (xe > end) {
    xe = end;
  }

  const float *srcf = src + skip * sizeOT;

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = xs; ix < xe; ix++) {
    _write_inner(srcf + (ix - xs) * sizeOT, linfo.xinfos[ix - start], os, oe,
                 ts, te);
  }
}

void SegyRW::find_nearest_idx(LineInfo &linfo, size_t xs, size_t xe,
                                     size_t &its, size_t &ite) {
  size_t start = linfo.isline ? xl2ix(linfo.lstart) : of2io(linfo.lstart);
  size_t end = linfo.isline ? xl2ix(linfo.lend) + 1 : of2io(linfo.lend) + 1;
  if (xs < start) {
    xs = start;
  }
  if (xe > end) {
    xe = end;
  }

  auto iindex = [&](auto it) {
    if (linfo.isline) {
      return itx2ix(it);
    } else {
      return itx2io(it);
    }
  };

  // find start idx
  its = linfo.itstart + (xs - start);
  its = its > linfo.itend ? linfo.itend : its;


  if (iindex(its) > xs) {
    while (iindex(its - 1) >= xs && (its - 1) >= linfo.itstart) {
      its--;
    } // its == linfo.itstart?
  } else if (its + 1 <= linfo.itend && iindex(its + 1) <= xs) {
    while (iindex(its + 1) <= xs && (its + 1) <= linfo.itend) {
      its++;
    }
  }

  // find end idx
  ite = its + (xe - xs);
  ite = ite > linfo.itend ? linfo.itend : ite;

  if (iindex(ite) > xe) {
    while (iindex(ite - 1) >= xe && (ite - 1) >= linfo.itstart) {
      ite--;
    }
  } else if (ite + 1 <= linfo.itend && iindex(ite + 1) <= xe) {
    while (iindex(ite + 1) <= xe && (ite + 1) <= linfo.itend) {
      ite++;
    }
  }
  if (ite == linfo.itend && iindex(ite) < xe) {
    ite++;
  }
}


uint64_t SegyRW::_copy_inner(char *dst, const float *src, LineInfo &linfo,
                             size_t ks, size_t ke, size_t ts, size_t te,
                             bool fromsrc) {
  char *odst = dst;
  size_t nt = te - ts;
  bool tchanged = nt == m_meta.nt ? false : true;

  /* no data to copy, just fill */
  if ((linfo.count == 0 && linfo.idx.size() == 0) || NoOverlap(linfo, ks, ke)) {
    return 0;
  }

  /* line is not continuous */
  else if (linfo.count == 0 && linfo.idx.size() > 0) {
    for (size_t ik = ks; ik < ke; ik++) {
      size_t it = linfo.idx[ik];
      if (it != kInvalid) {
        // copy trace header
        memcpy(dst, trheader(it), kTraceHeaderSize);
        if (tchanged) {
          if (ts > 0) { // : ts need *1000 ?
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

  else if (linfo.count == kInvalid) {
    // find start & end
    size_t its, ite;
    find_nearest_idx(linfo, ks, ke, its, ite);
    for (size_t it = its; it < ite; it++) {
      // copy trace header
      memcpy(dst, trheader(it), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) { // : ts need *1000 ?
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
        if (linfo.isline) {
          m_wfunc(dst, src + (itx2ix(it) - ks) * nt, nt);
        } else {
          m_wfunc(dst, src + (itx2io(it) - ks) * nt, nt);
        }
      }
      dst += nt * m_meta.esize;
    }
    return dst - odst;
  }

  else {
    std::array<size_t, 4> idx;
    find_idx(idx, linfo, ks, ke);
    const float *srcf = src + idx[0] * nt;

    for (size_t tx = idx[2]; tx < idx[3]; tx++) {
      // copy trace header
      memcpy(dst, trheader(tx), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) {
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
                            size_t xs, size_t xe, size_t os, size_t oe,
                            size_t ts, size_t te, bool fromsrc) {
  size_t nt = te - ts;
  size_t no = oe - os;
  // size_t nx = xe - xs;
  uint64_t sizeOT = nt * no;
  // uint64_t sizeXOT = nx * sizeOT;

  /* no data to copy, just fill */
  if (linfo.count == 0 || NoOverlap(linfo, xs, xe)) {
    return 0;
  }

  uint64_t jump = 0;

  size_t start = xl2ix(linfo.lstart);
  size_t end = xl2ix(linfo.lend) + 1;
  size_t skip = 0;

  if (xs < start) {
    skip = start - xs;
    xs = start;
  }
  if (xe > end) {
    xe = end;
  }

  char *odst = dst;
  const float *srcf = nullptr;
  if (!fromsrc) {
    srcf = src + skip * sizeOT;
  }

  // read, when is 4D, xlines are always continuous as we filled
  for (size_t ix = xs; ix < xe; ix++) {
    jump = _copy_inner(dst, srcf + (ix - xs) * sizeOT,
                       linfo.xinfos[ix - start], os, oe, ts, te, fromsrc);
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
  if (m_w) {
    m_meta.ntrace = (m_sink.size() - kTraceHeaderStart) / m_meta.tracesize;
  } else {
    m_meta.ntrace = (m_src.size() - kTraceHeaderStart) / m_meta.tracesize;
  }

  m_meta.trace_sorting_code = bkeyi2(kBTraceSortingCodeField);
  setRWFunc(m_meta.dformat);
}

void SegyRW::_create_from_segy(const std::string &outname, const float *src,
                               const std::vector<size_t> &ranges, bool is2d,
                               const std::string &textual, bool fromsrc,
                               uint64_t check_size) {
  size_t rsize = ranges.size();
  if (rsize < 4) {
    throw std::runtime_error(
        "Error ranges size, ranges.size() must >= 4, but got " +
        std::to_string(rsize));
  }
  size_t ts = ranges[rsize - 2], te = ranges[rsize - 1];
  if (ts > te || te > m_meta.nt) {
    throw std::runtime_error("ts and te index out of range");
  }

  size_t tstart = 0, tend = 0, is = 0, ie = 0, xs = 0, xe = 0, os = 0, oe = 0;
  uint64_t maxsize = 0, tracesize = 0;
  if (is2d || m_ndim == 2) {
    if (ranges.size() != 4) {
      throw std::runtime_error("Error ranges. you set (is2d || ndim==2), but "
                               "ranges.size() != 4, ranges.size() = " +
                               std::to_string(ranges.size()));
    }
    tstart = ranges[0];
    tend = ranges[1];
    if (check_size > 0 &&
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
      if (check_size > 0 && check_size != (uint64_t)(ie - is) * (xe - xs) *
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
      if (check_size > 0 && check_size != (uint64_t)(ie - is) * (xe - xs) *
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
  size_t nt = te - ts;

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
    int step = cal_progress_steps(tend - tstart, show_progress_, 10000, 50);
    for (size_t it = tstart; it < tend; it++) {
      g_check_signals_callback();
      if (step && (it - tstart) % step == 0) {
        g_progress_callback(it - tstart, tend - tstart);
      }
      memcpy(outptr, trheader(it), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) {
          segy::set_keyi2(outptr, kTStartTimeField, ts);
        }
        segy::set_keyi2(outptr, kTSampleCountField, nt);
      }
      outptr += kTraceHeaderSize;
      memcpy(outptr, trDataStart(it, ts), nt * m_meta.esize);
      outptr += nt * m_meta.esize;
    }
    if (step) {
      g_progress_callback(-1, tend - tstart);
    }
  } else {
    uint64_t jump = 0;
    int step = cal_progress_steps(ie - is, show_progress_, 10000, 50);
    for (size_t ii = is; ii < ie; ii++) {
      g_check_signals_callback();
      if (step && (ii - is) % step == 0) {
        g_progress_callback(ii - is, ie - is);
      }
      LineInfo &linfo = m_iinfos[ii];
      if (linfo.count == 0 && linfo.idx.size() == 0) {
        continue;
      }
      if (m_ndim == 3) {
        jump = _copy_inner(outptr, src, linfo, xs, xe, ts, te, fromsrc);
        outptr += jump;
        if (!fromsrc) {
          src += (xe - xs) * nt;
        }
      } else {
        jump = _copy4d_xo(outptr, src, linfo, xs, xe, os, oe, ts, te, fromsrc);
        outptr += jump;
        if (!fromsrc) {
          src += (xe - xs) * (oe - os) * nt;
        }
      }
    }
    if (step) {
      g_progress_callback(-1, ie - is);
    }
  }

  int64_t to_remove = maxsize - (outptr - rw_mmap.data());
  rw_mmap.unmap();
  if (to_remove > 0) {
    truncate_file(outname, to_remove);
  }
}

void SegyRW::cut(const std::string &segy_name,
                 const std::vector<size_t> &ranges, bool is2d,
                 const std::string &textual) {
  _create_from_segy(segy_name, nullptr, ranges, is2d, textual, true);
}

void SegyRW::create_by_sharing_header(const std::string &segy_name,
                                      const float *src,
                                      const std::vector<size_t> &shape,
                                      const std::vector<size_t> &start,
                                      bool is2d, const std::string &textual) {
  size_t nd = is2d ? 2 : m_ndim;

  if (shape.size() != start.size() || shape.size() != nd) {
    throw std::runtime_error("shape's size must be equal to the start's size, "
                             "and shape's size must be equal to ndim");
  }
  std::vector<size_t> ranges(nd * 2);
  for (size_t i = 0; i < nd; ++i) {
    ranges[i * 2] = start[i];
    ranges[i * 2 + 1] = start[i] + shape[i];
  }
  _create_from_segy(segy_name, src, ranges, is2d, textual, false);
}

void SegyRW::create_by_sharing_header(const std::string &segy_name,
                                      const std::string &src_name,
                                      const std::vector<size_t> &shape,
                                      const std::vector<size_t> &start,
                                      bool is2d, const std::string &textual) {
  std::error_code error;
  mio::mmap_source float_file;
  float_file.map(src_name, error);
  if (error) {
    throw std::runtime_error("Cannot mmap binary file in 'r' mode: " +
                             src_name);
  }

  size_t nd = is2d ? 2 : m_ndim;
  if (shape.size() != start.size() || shape.size() != nd) {
    throw std::runtime_error("shape's size must be equal to the start's size, "
                             "and shape's size must be equal to ndim");
  }
  std::vector<size_t> ranges(nd * 2);
  for (size_t i = 0; i < nd; ++i) {
    ranges[i * 2] = start[i];
    ranges[i * 2 + 1] = start[i] + shape[i];
  }

  const float *src = reinterpret_cast<const float *>(float_file.data());
  _create_from_segy(segy_name, src, ranges, is2d, textual, false,
                    float_file.size());
}

inline static void modify_keys(char *dst, const int32_t *keys, size_t keysize) {
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
                 const int32_t *keys, const std::vector<size_t> &shape,
                 const std::string &textual, const uchar *bheader,
                 const uchar *theader, size_t keysize) {
  size_t ndim = shape.size();
  size_t ni, nx, no;
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
  size_t esize = 0;

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

  int step = cal_progress_steps(ntrace, true, 10000, 50);
  for (size_t i = 0; i < ntrace; i++) {
    g_check_signals_callback();
    if (step && i % step == 0) {
      g_progress_callback(i, ntrace);
    }
    memcpy(dst, theader, kTraceHeaderSize);
    modify_keys(dst, keys + i * (uint64_t)keysize, keysize);
    dst += kTraceHeaderSize;
    wfunc(dst, src + i * nt, nt);
    dst += nt * esize;
  }
  if (step) {
    g_progress_callback(-1, ntrace);
  }

  assert(dst - rw_mmap.data() == needsize);
  rw_mmap.unmap();
}

} // namespace segy
