/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

#ifndef CIG_SEGY_BASE_H
#define CIG_SEGY_BASE_H

#include "mio.hpp"
#include "sutils.hpp"

#include <_types/_uint64_t.h>
#include <vector>

namespace segy {

struct KeyLocs {
  int iline = 189; // inline number location
  int xline = 193; // crossline number location
  int offset = 37; // offset number location
  int xloc = 181;  // X location
  int yloc = 185;  // Y location
  int istep = 1;   // inline step
  int xstep = 1;   // crossline step
  int ostep = 1;   // offset step
};

struct MetaInfo {
  // count information
  int32_t nt = 0;         // number of samples in time dimension
  int32_t no = 1;         // number of samples in offset dimension
  int32_t nx = 0;         // number of samples in crosline dimension
  int32_t ni = 0;         // number of samples in inline dimension
  int64_t ntrace = 0;     // number of traces
  uint64_t tracesize = 0; // trace size in char

  int16_t dt = 0; // interval of time/depth
  float dx = 0;   // interval of crossline
  float di = 0;   // interval of inline

  int16_t scalar = 0;
  int16_t dformat = 0;

  int16_t start_time = 0;
  int32_t start_iline = 0;
  int32_t end_iline = 0;
  int32_t start_xline = 0;
  int32_t end_xline = 0;
  int32_t start_offset = 0;
  int32_t end_offset = 0;

  int32_t trace_sorting_code = 0;
  int32_t esize = 4;
  float fillNoValue = 0;
};

class SegyBase {
public:
  using ReadFunc = std::function<void(float *, const char *, int)>;
  using WriteFunc = std::function<void(char *, const float *, int)>;

  SegyBase(const std::string &segyname)
      : segyname(segyname), m_data_ptr(nullptr) {}

  ~SegyBase() { this->close_file(); }

  void close_file() {
    if (m_sink.is_mapped()) {
      m_sink.unmap();
    }
  }

  inline int64_t ntrace() const { return m_meta.ntrace; }
  inline int nt() const { return m_meta.nt; }

  void setLocations(int iline, int xline, int offset = 37);
  void setInlineLocation(int loc);
  void setCrosslineLocation(int loc);
  void setOffsetLocation(int loc);

  void setXLocation(int loc);
  void setYLocation(int loc);
  void setXYLocations(int xloc, int yloc);

  void setInlineStep(int step);
  void setCrosslineStep(int step);
  void setOffsetStep(int step);
  void setSteps(int istep, int xstep, int ostep = 1);
  void setFill(float fill) { m_meta.fillNoValue = fill; }

  std::string textual_header(char coding = 'u');

  // binary header
  int16_t bkeyi2(int loc);
  int32_t bkeyi4(int loc);

  // trace header
  int16_t keyi2(int n, int loc);
  int32_t keyi4(int n, int loc);
  inline int32_t iline(int n) { return keyi4(n, m_keys.iline); }
  inline int32_t xline(int n) { return keyi4(n, m_keys.xline); }
  inline int32_t offset(int n) { return keyi4(n, m_keys.offset); }
  inline int32_t coordx(int n) { return keyi4(n, m_keys.xloc); }
  inline int32_t coordy(int n) { return keyi4(n, m_keys.yloc); }

  void get_binary_header(uchar *binheader);
  void get_trace_header(uchar *traceheader, int n);

  void get_trace_keys(int *dst, const std::vector<int> &keys,
                      const std::vector<int> &length, int beg, int end);

  void itrace(float *data, int n);
  void collect(float *data, int beg, int end, int tbeg, int tend);
  void collect(float *data, const int32_t *index, int n, int tbeg, int tend);

  // W mode
  void set_bkeyi2(int loc, int16_t val);
  void set_bkeyi4(int loc, int32_t val);

  void set_keyi2(int n, int loc, int16_t val);
  void set_keyi4(int n, int loc, int32_t val);

  inline void set_iline(int n, int32_t val) { set_keyi4(n, m_keys.iline, val); }
  inline void set_xline(int n, int32_t val) { set_keyi4(n, m_keys.xline, val); }
  inline void set_offset(int n, int32_t val) {
    set_keyi4(n, m_keys.offset, val);
  }
  inline void set_coordx(int n, int32_t val) { set_keyi4(n, m_keys.xloc, val); }
  inline void set_coordy(int n, int32_t val) { set_keyi4(n, m_keys.yloc, val); }

  // void set_trace_keys(const int *dst, const std::vector<int> &keys,
  //                     const std::vector<int> &length, int beg, int end);

  void write_itrace(const float *data, int n);
  void write_traces(const float *data, int beg, int end, int tbeg, int tend);
  void write_traces(const float *data, const int32_t *index, int n, int tbeg,
                    int tend);

protected:
  std::string segyname;
  const char *m_data_ptr;
  mio::mmap_sink m_sink;
  KeyLocs m_keys;
  MetaInfo m_meta;
  ReadFunc m_readfunc;
  WriteFunc m_wfunc;

  inline const char *brheader() { return m_data_ptr + kTextualHeaderSize; }
  inline const char *trheader(int n) {
    return m_data_ptr + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline const char *trDataStart(int n) {
    return trheader(n) + kTraceHeaderSize;
  }

  inline char *bwheader() { return m_sink.data() + kTextualHeaderSize; }
  inline char *twheader(int n) {
    return m_sink.data() + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline char *twDataStart(int n) { return twheader(n) + kTraceHeaderSize; }

  void setRWFunc(int dformat);
};

} // namespace segy
#endif