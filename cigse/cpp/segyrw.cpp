/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

#include "segyrw.h"
#include <iostream>

namespace segy {

// modify header keys function, for cut, create_by_sharing_header
inline static void set_bkeyi2(char *bheader, int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(bheader + loc - 1) = swap_endian<int16_t>(val);
}
inline static void set_keyi2(char *theader, int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(theader + loc - 1) = swap_endian<int16_t>(val);
}

void SegyRW::set_segy_type(int ndim) {
  if (ndim < 2 || ndim > 4) {
    std::runtime_error("Error SEG-Y type, 2 for 2D, 3 for poststack, 4 for "
                       "prestack. But now is " +
                       std::to_string(ndim));
  }
  m_ndim = ndim;
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

    if (iline(it) == iiline || iline(it - 1) != iiline) {
      std::cout << iline(it) << " " << iiline << " " << iline(it - 1)
                << std::endl;
      throw std::runtime_error("error1");
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
      // assert (m_iinfos[ii].xend - xs) % xstep == 0
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

        if (xline(xt) == xxline || xline(xt - 1) != xxline) {
          std::cout << xline(xt) << " " << xxline << " " << xline(xt - 1)
                    << std::endl;
          throw std::runtime_error("error2");
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
            std::cout << xline(xt) << " " << xxline << " " << xstep
                      << std::endl;
            std::cout << "skipx: " << skipx << " xstep: " << xstep << std::endl;
            throw std::runtime_error(
                "Cannot analysis this xline, maybe xstep is wrong");
          }
          skipx /= xstep;
        }
      }
    }

    // missing line
    if (iline(it) != (iiline + istep)) {
      skipi = xline(it) - (iiline + istep);
      if (skipi % istep != 0) {
        std::cout << iline(it) << " " << iiline << " " << istep << std::endl;
        std::cout << "skipi: " << skipi << " istep: " << istep << std::endl;
        throw std::runtime_error("Cannot analysis line, maybe istep is wrong");
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
      if (linfo.count == (linfo.lend - linfo.lstart) / m_keys.xstep + 1) {
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
        if (xinfo.count == (xinfo.lend - xinfo.lstart) / m_keys.ostep + 1) {
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

std::vector<int> SegyRW::shape() const {
  if (m_ndim == 2) {
    return {(int)m_meta.ntrace, m_meta.nt};
  } else if (m_ndim == 3) {
    return {m_meta.ni, m_meta.nx, m_meta.nt};
  } else {
    return {m_meta.ni, m_meta.nx, m_meta.no, m_meta.nt};
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

void SegyRW::tofile(const std::string &binary_out_name, bool is2d) {
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
void SegyRW::_read_inner(float *dst, LineInfo &linfo, int ks, int ke, int ts,
                         int te) {
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
        m_readfunc(dst, trDataStart(it) + ts, nt);
      }
      dst += nt;
    }
    return;
  }

  /* continuous line TODO: remove if confition?*/
  else if (linfo.count > 0 && isCtnL(linfo)) {
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
    }

    // assert
    return;
  } else {
    throw std::runtime_error("error count");
  }
}

void SegyRW::_read4d_xo(float *dst, LineInfo &linfo, int xs, int xe, int os,
                        int oe, int ts, int te) {
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
  // assert
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
            segy::set_keyi2(dst, kTStartTimeField, ts);
          }
          segy::set_keyi2(dst, kTSampleCountField, nt);
        }
        dst += kTraceHeaderSize;

        // copy trace
        if (fromsrc) {
          // copy from src, i.e., cut
          memcpy(dst, trDataStart(it) + ts, nt * m_meta.esize);
        } else {
          // copy from data, i.e., create_by_sharing_header
          m_wfunc(dst, src + (ik - ks) * nt, nt);
        }
        dst += nt * m_meta.esize;
      }
    }
    return dst - odst;
  }

  /* continuous line TODO: remove if confition?*/
  else if (linfo.count > 0 && isCtnL(linfo)) {
    std::array<int32_t, 4> idx;
    find_idx(idx, linfo, ks, ke);
    const float *srcf = src + idx[0] * nt;

    for (size_t tx = idx[2]; tx < idx[3]; tx++) {
      // copy trace header
      memcpy(dst, trheader(tx), kTraceHeaderSize);
      if (tchanged) {
        if (ts > 0) { // TODO: ts need *1000 ?
          segy::set_keyi2(dst, kTStartTimeField, ts);
        }
        segy::set_keyi2(dst, kTSampleCountField, nt);
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

bool SegyRW::isPreStack() {
  int o0 = offset(0);
  int n = m_meta.ntrace < 500 ? m_meta.ntrace : 500;
  for (size_t i = 1; i < n; i++) {
    if (o0 != offset(i)) {
      return true;
    }
  }
  return false;
}

bool SegyRW::isCtnL(LineInfo &linfo) {
  if (linfo.isline) {
    return linfo.count == ((linfo.lend - linfo.lstart) / m_keys.xstep + 1);
  } else {
    return linfo.count == ((linfo.lend - linfo.lstart) / m_keys.ostep + 1);
  }
}

uint64_t SegyRW::_need_size_b(int ndim) {
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

uint64_t SegyRW::_need_size_s(uint64_t ntrace) {
  return kTraceHeaderStart + m_meta.tracesize * ntrace;
}

bool SegyRW::NoOverlap(LineInfo &linfo, int s, int e) {
  if (linfo.isline) {
    return s > xl2ix(linfo.lend) || e <= xl2ix(linfo.lstart);
  } else {
    return s > of2io(linfo.lend) || e <= of2io(linfo.lstart);
  }
}

// Use only if xs/xe and xinfo have overlapped parts
void SegyRW::find_idx(std::array<int32_t, 4> &idx, LineInfo &linfo, int xs,
                      int xe) {
  memset(&idx, 0, 4 * 4);
  int start, end;
  if (linfo.isline) {
    start = xl2ix(linfo.lstart);
    end = xl2ix(linfo.lend) + 1;
  } else {
    start = of2io(linfo.lstart);
    end = of2io(linfo.lend) + 1;
  }
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

void SegyRW::scanBinaryHeader() {
  int dformat = bkeyi2(kBSampleFormatField);
  auto it = kElementSize.find(dformat);
  if (it != kElementSize.end()) {
    m_meta.esize = it->second;
  } else {
    throw std::runtime_error("Unknown data format");
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
                               const std::string &textual, bool fromsrc) {

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
    tracesize = kTraceHeaderSize + (te - ts) * m_meta.esize;
    maxsize = kTraceHeaderStart + tracesize * (tend - tstart);

  } else {
    is = ranges[0];
    ie = ranges[1];
    xs = ranges[2];
    xe = ranges[3];

    if (m_ndim == 3) {
      if (ranges.size() != 8) {
        throw std::runtime_error("error ranges");
      }
      ts = ranges[4];
      te = ranges[5];
      tracesize = kTraceHeaderSize + (te - ts) * m_meta.esize;
      maxsize = kTraceHeaderStart + tracesize * (ie - is) * (xe - xs);
    } else {
      if (ranges.size() != 8) {
        throw std::runtime_error("error ranges");
      }
      os = ranges[4];
      oe = ranges[5];
      ts = ranges[6];
      te = ranges[7];
      tracesize = kTraceHeaderSize + (te - ts) * m_meta.esize;
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
                                      const std::vector<int> &ranges, bool is2d,
                                      const std::string &textual) {
  _create_from_segy(segy_name, src, ranges, is2d, textual, false);
}

void SegyRW::create_by_sharing_header(const std::string &segy_name,
                                      const std::string &src_name,
                                      const std::vector<int> &shape,
                                      const std::vector<int> &ranges, bool is2d,
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
