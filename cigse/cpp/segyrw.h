/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

/*
NOTE:
Be careful when using this class, we don't check the boundary of the input
parameters. Because we check it in wraping file or in python code.
So, if you want to use this class in other place, you should check the boundary.
*/

#ifndef CIG_SEGY_RW_H
#define CIG_SEGY_RW_H

#include "mio.hpp"
#include "segybase.hpp"
#include "utils.hpp"
#include <array>
#include <stdexcept>
#include <vector>

#define MMAX(a, b) ((a) > (b) ? (a) : (b))
#define MMIN(a, b) ((a) < (b) ? (a) : (b))

namespace segy {

struct LineInfo {

  /*
  (count == 0 && idx.size() == 0) ==> this line is missing
  (count == 0 && idx.size() > 0) ==> this line is not continuous and we fill idx
  for fast reading (count == Invalid) ==> this line is continuous, but we don't
  set idx (count > 0) ==> this line is continuous and we set count
  */

  bool isline;
  size_t line;    // iline/xline number for this line
  size_t itstart; // for this line, start trace index
  size_t itend;   // for this line, end trace index
  size_t count;   // the number of traces for this line
  size_t lstart;  // the start xline/offset
  size_t lend;    // the end xline/offset

  // if not continuous, record all idx
  std::vector<size_t> idx;

  std::vector<LineInfo> xinfos; // XLineInfos for Line

  LineInfo(bool t)
      : isline(t), itstart(kInvalid), itend(kInvalid), count(kInvalid),
        lstart(kInvalid), lend(kInvalid) {}

  void set_xinfos(size_t size) { xinfos.resize(size, LineInfo(false)); }
};

class SegyRW : public SegyBase {
public:
  SegyRW(const std::string &segyname, bool write = false) : SegyBase() {
    this->m_w = write;
    std::error_code error;

    if (write) {
      this->m_sink.map(segyname, error);
      if (error) {
        throw std::runtime_error("Cannot mmap the file as RW mode: " +
                                 segyname + ", message: " + error.message());
      }
      m_data_ptr = this->m_sink.data();
    } else {
      this->m_src.map(segyname, error);
      if (error) {
        throw std::runtime_error("Cannot mmap the file as R mode: " + segyname +
                                 ", message: " + error.message());
      }
      m_data_ptr = this->m_src.data();
    }

    scanBinaryHeader();
  }

  void set_segy_type(size_t ndim);
  void scan();
  std::vector<size_t> shape() const;
  inline size_t ndim() const { return m_ndim; }
  inline bool is_scanned() const { return isScan; }

  // R mode
  void read4d(float *dst, size_t ib, size_t ie, size_t xb, size_t xe, size_t ob,
              size_t oe, size_t tb, size_t te);
  void read3d(float *dst, size_t ib, size_t ie, size_t xb, size_t xe, size_t tb,
              size_t te);
  void read(float *dst);
  void read_tslice(float *dst, size_t it, size_t stepi = 1, size_t stepx = 1);
  void tofile(const std::string &binary_out_name, bool is2d = false);
  void cut(const std::string &outname, const std::vector<size_t> &ranges,
           bool is2d = false, const std::string &textual = "");
  void create_by_sharing_header(const std::string &segy_name, const float *src,
                                const std::vector<size_t> &shape,
                                const std::vector<size_t> &start,
                                bool is2d = false,
                                const std::string &textual = "");
  void create_by_sharing_header(const std::string &segy_name,
                                const std::string &src_name,
                                const std::vector<size_t> &shape,
                                const std::vector<size_t> &start,
                                bool is2d = false,
                                const std::string &textual = "");

  // write mode
  void write(const float *data);
  void write3d(const float *data, size_t ib, size_t ie, size_t xb, size_t xe,
               size_t tb, size_t te);
  void write4d(const float *data, size_t ib, size_t ie, size_t xb, size_t xe,
               size_t ob, size_t oe, size_t tb, size_t te);

  // Geometry view?

protected:
  std::vector<LineInfo> m_iinfos;

private:
  size_t m_ndim = 2;
  bool isScan = false;

  bool isPreStack();
  inline size_t xl2ix(size_t xl) {
    return (xl - m_meta.start_xline) / m_keys.xstep;
  }
  inline size_t of2io(size_t of) {
    return (of - m_meta.start_offset) / m_keys.ostep;
  }
  inline size_t itx2ix(size_t it) { return xl2ix(xline(it)); }
  inline size_t itx2io(size_t it) { return of2io(offset(it)); }
  bool NoOverlap(LineInfo &linfo, size_t s, size_t e);
  void find_idx(std::array<size_t, 4> &idx, LineInfo &linfo, size_t xs,
                size_t xe);
  void find_nearest_idx(LineInfo &linfo, size_t xs, size_t xe, size_t &its,
                        size_t &ite);

  void _read_inner(float *dst, LineInfo &linfo, size_t ks, size_t ke, size_t ts,
                   size_t te);
  void _read4d_xo(float *dst, LineInfo &linfo, size_t xs, size_t xe, size_t os,
                  size_t oe, size_t ts, size_t te);
  uint64_t _copy_inner(char *dst, const float *src, LineInfo &linfo, size_t ks,
                       size_t ke, size_t ts, size_t te, bool fromsrc);
  uint64_t _copy4d_xo(char *dst, const float *src, LineInfo &linfo, size_t xs,
                      size_t xe, size_t os, size_t oe, size_t ts, size_t te,
                      bool fromsrc);
  void _create_from_segy(const std::string &outname, const float *src,
                         const std::vector<size_t> &ranges, bool is2d,
                         const std::string &textual, bool fromsrc,
                         uint64_t check_size = 0);

  void _write_inner(const float *src, LineInfo &linfo, size_t ks, size_t ke,
                    size_t ts, size_t te);
  void _write4d_xo(const float *src, LineInfo &linfo, size_t xs, size_t xe,
                   size_t os, size_t oe, size_t ts, size_t te);

  uint64_t _need_size_b(size_t ndim = kInvalid);
  void scanBinaryHeader();
};

inline void SegyRW::set_segy_type(size_t ndim) {
  if (ndim < 2 || ndim > 4) {
    std::runtime_error("Error SEG-Y type, 2 for 2D, 3 for poststack, 4 for "
                       "prestack. But now is " +
                       std::to_string(ndim));
  }
  m_ndim = ndim;
}

inline std::vector<size_t> SegyRW::shape() const {
  if (m_ndim == 2) {
    return {static_cast<size_t>(m_meta.ntrace), m_meta.nt};
  } else if (m_ndim == 3) {
    return {m_meta.ni, m_meta.nx, m_meta.nt};
  } else {
    return {m_meta.ni, m_meta.nx, m_meta.no, m_meta.nt};
  }
}

// This function is too simple, so we recommend to
// set dim by `set_segy_type` function
inline bool SegyRW::isPreStack() {
  int o0 = offset(0);
  size_t n = m_meta.ntrace < 500 ? m_meta.ntrace : 500;
  for (size_t i = 1; i < n; i++) {
    if (o0 != offset(i)) {
      return true;
    }
  }
  return false;
}

inline uint64_t SegyRW::_need_size_b(size_t ndim) {
  uint64_t need_size = 0;
  if (ndim > 4) {
    ndim = m_ndim;
  }

  // priviate function, so we don't check the ndim
  if (ndim == 2) {
    need_size = m_meta.ntrace * m_meta.nt;
  } else if (ndim > 2) {
    need_size = static_cast<uint64_t>(m_meta.ni) * m_meta.nx * m_meta.nt;
  }

  if (ndim == 4) {
    need_size = need_size * m_meta.no;
  }

  return need_size * sizeof(float);
}

inline bool SegyRW::NoOverlap(LineInfo &linfo, size_t s, size_t e) {
  if (linfo.isline) {
    return s > xl2ix(linfo.lend) || e <= xl2ix(linfo.lstart);
  } else {
    return s > of2io(linfo.lend) || e <= of2io(linfo.lstart);
  }
}

// Use only if xs/xe and xinfo have overlapped parts
inline void SegyRW::find_idx(std::array<size_t, 4> &idx, LineInfo &linfo,
                             size_t xs, size_t xe) {
  memset(&idx, 0, 4 * 4);
  size_t start = linfo.isline ? xl2ix(linfo.lstart) : of2io(linfo.lstart);
  size_t end = linfo.isline ? xl2ix(linfo.lend) + 1 : of2io(linfo.lend) + 1;
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



void create_segy(const std::string &segyname, const float *src,
                 const int32_t *keys, const std::vector<size_t> &shape,
                 const std::string &textual, const uchar *bheader,
                 const uchar *theader, size_t keysize);

} // namespace segy
#endif