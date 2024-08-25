/*********************************************************************
** Copyright (c) 2023 Jintao Li.
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
#include "sutils.h"

#include <fstream>
#include <stdexcept>
#include <vector>

#define MMAX(a, b) ((a) > (b) ? (a) : (b))
#define MMIN(a, b) ((a) < (b) ? (a) : (b))

namespace segy {

struct LineInfo {

  bool isline;
  int line;    // iline/xline number for this line
  int itstart; // for this line, start trace index
  int itend;   // for this line, end trace index
  int count;   // the number of traces for this line
  int lstart;  // the start xline/offset
  int lend;    // the end xline/offset

  // if not continuous, record all idx
  std::vector<int> idx;

  std::vector<LineInfo> xinfos; // XLineInfos for Line

  LineInfo(bool t)
      : isline(t), itstart(-1), itend(-1), count(-1), lstart(-1), lend(-1) {}

  void set_xinfos(size_t size) { xinfos.resize(size, LineInfo(false)); }
};

class SegyRW : public SegyBase {
public:
  explicit SegyRW(const std::string &segyname);
  ~SegyRW() override;

  void set_segy_type(int ndim);
  void scan();
  std::vector<int> shape();

  // R mode
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

  // W mode
  void set_bkeyi2(int loc, int16_t val);
  void set_bkeyi4(int loc, int32_t val);
  // void set_bkeyi8(int loc, int64_t val);

  void set_keyi2(int n, int loc, int16_t val);
  void set_keyi4(int n, int loc, int32_t val);
  // void set_keyi8(int n, int loc, int64_t val);

  inline void set_iline(int n, int32_t val) { set_keyi4(n, m_keys.iline, val); }
  inline void set_xline(int n, int32_t val) { set_keyi4(n, m_keys.xline, val); }
  inline void set_offset(int n, int32_t val) {
    set_keyi4(n, m_keys.offset, val);
  }
  inline void set_coordx(int n, int32_t val) { set_keyi4(n, m_keys.xloc, val); }
  inline void set_coordy(int n, int32_t val) { set_keyi4(n, m_keys.yloc, val); }

  void set_trace_keys(const int *dst, const std::vector<int> &keys,
                      const std::vector<int> &length, int beg, int end);

  void write_itrace(const float *data, int n);
  void write_traces(const float *data, int beg, int end, int tbeg, int tend);
  void write_traces(const float *data, const int32_t *index, int n, int tbeg,
                    int tend);

  void write(const float *data);
  void write3d(const float *data, int is, int ie, int xs, int xe, int ts,
               int te);
  void write4d(const float *data, int is, int ie, int xs, int xe, int os,
               int oe, int ts, int te);

protected:
  mio::mmap_source m_src;
  mio::mmap_sink m_sink;
  std::vector<LineInfo> m_iinfos;
  int m_ndim = 2;
  void scanBinaryHeader();

  inline char *bwheader() { return m_sink.data() + kTextualHeaderSize; }
  inline char *twheader(int n) {
    return m_sink.data() + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline char *twDataStart(int n) { return twheader(n) + kTraceHeaderSize; }

private:
  bool isScan = false;
  uint64_t m_sizeL = 1;

  bool isPreStack();
  // inline int il2ii(int il) { return (il - m_meta.start_iline) / m_keys.istep;
  // } inline int ii2il(int ii) { return ii * m_keys.istep + m_meta.start_iline;
  // }
  inline int xl2ix(int xl) { return (xl - m_meta.start_xline) / m_keys.xstep; }
  // inline int ix2xl(int ix) { return ix * m_keys.xstep + m_meta.start_xline; }
  inline int of2io(int of) { return (of - m_meta.start_offset) / m_keys.ostep; }
  // inline int io2of(int io) { return io * m_keys.ostep + m_meta.start_offset;
  // }
  bool isCtnL(LineInfo &linfo);
  bool NoOverlap(LineInfo &linfo, int s, int e);
  void find_idx(std::array<int32_t, 4> &idx, LineInfo &linfo, int xs, int xe);

  void _read_inner(float *dst, LineInfo &linfo, int ks, int ke, int ts, int te);
  void _read4d_xo(float *dst, LineInfo &linfo, int xs, int xe, int os, int oe,
                  int ts, int te);
  uint64_t _copy_inner(char *dst, const float *src, LineInfo &linfo, int ks,
                       int ke, int ts, int te, bool fromsrc);
  uint64_t _copy4d_xo(char *dst, const float *src, LineInfo &linfo, int xs,
                      int xe, int os, int oe, int ts, int te, bool fromsrc);
  void _create_from_segy(const std::string &outname, const float *src,
                         const std::vector<int> &ranges, bool is2d,
                         const std::string &textual, bool fromsrc);

  uint64_t _need_size_b(int ndim = -1);
  uint64_t _need_size_s(uint64_t ntrace);
};

} // namespace segy
#endif