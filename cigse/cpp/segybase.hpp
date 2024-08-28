/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

#ifndef CIG_SEGY_BASE_H
#define CIG_SEGY_BASE_H

#include "mio.hpp"
#include "utils.hpp"

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
  SegyBase() : m_data_ptr(nullptr) {}

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
  const char *m_data_ptr;
  mio::mmap_sink m_sink;
  KeyLocs m_keys;
  MetaInfo m_meta;
  ReadFunc m_readfunc;
  WriteFunc m_wfunc;

  inline const char *brheader() const {
    return m_data_ptr + kTextualHeaderSize;
  }
  inline const char *trheader(int n) const {
    return m_data_ptr + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline const char *trDataStart(int n, int tbeg = 0) const {
    return trheader(n) + kTraceHeaderSize + tbeg * m_meta.esize;
  }

  inline char *bwheader() { return m_sink.data() + kTextualHeaderSize; }
  inline char *twheader(int n) {
    return m_sink.data() + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline char *twDataStart(int n, int tbeg = 0) {
    return twheader(n) + kTraceHeaderSize + tbeg * m_meta.esize;
  }

  void setRWFunc(int dformat);
};

/******************************************************/
/**************** For Public  *************************/
/******************************************************/

/**************** set locations  *************************/

inline void SegyBase::setLocations(int iloc, int xloc, int oloc) {
  m_keys.iline = iloc;
  m_keys.xline = xloc;
  m_keys.offset = oloc;
}
inline void SegyBase::setInlineLocation(int loc) { m_keys.iline = loc; }
inline void SegyBase::setCrosslineLocation(int loc) { m_keys.xline = loc; }
inline void SegyBase::setOffsetLocation(int loc) { m_keys.offset = loc; }

inline void SegyBase::setXLocation(int loc) { m_keys.xloc = loc; }
inline void SegyBase::setYLocation(int loc) { m_keys.yloc = loc; }
inline void SegyBase::setXYLocations(int xloc, int yloc) {
  m_keys.xloc = xloc;
  m_keys.yloc = yloc;
}

inline void SegyBase::setInlineStep(int step) { m_keys.istep = step; }
inline void SegyBase::setCrosslineStep(int step) { m_keys.xstep = step; }
inline void SegyBase::setOffsetStep(int step) { m_keys.ostep = step; }
inline void SegyBase::setSteps(int istep, int xstep, int ostep) {
  m_keys.istep = istep;
  m_keys.xstep = xstep;
  m_keys.ostep = ostep;
}

/******************************************************/
/**************** For read mode  **********************/
/******************************************************/

inline std::string SegyBase::textual_header(char coding) {
  const char *src = m_data_ptr;
  char out[kTextualHeaderSize + kTextualRows];
  bool isEBCDIC = true;
  if (coding == 'a') {
    isEBCDIC = false;
  } else if (coding == 'u') {
    isEBCDIC = isTextInEBCDICFormat(src, kTextualHeaderSize);
  } else {
    throw std::invalid_argument("Only support 'a' and 'u'");
  }
  for (int iRow = 0; iRow < kTextualRows; iRow++) {
    int offset = iRow * kTextualColumns;
    for (int iCol = 0; iCol < kTextualColumns; iCol++) {
      if (isEBCDIC) {
        out[iCol + offset + iRow] = getASCIIfromEBCDIC(src[iCol + offset]);
      } else {
        out[iCol + offset + iRow] = src[iCol + offset];
      }
    }
    if (iRow < kTextualRows - 1) {
      out[(iRow + 1) * (kTextualColumns + 1) - 1] = '\n';
    }
  }

  return std::string(out);
}

inline int16_t SegyBase::bkeyi2(int loc) {
  return swap_endian<int16_t>(brheader() + loc - 1);
}
inline int32_t SegyBase::bkeyi4(int loc) {
  return swap_endian<int32_t>(brheader() + loc - 1);
}

inline int16_t SegyBase::keyi2(int n, int loc) {
  return swap_endian<int16_t>(trheader(n) + loc - 1);
}
inline int32_t SegyBase::keyi4(int n, int loc) {
  return swap_endian<int32_t>(trheader(n) + loc - 1);
}

inline void SegyBase::get_binary_header(uchar *binheader) {
  memcpy(binheader, brheader(), kBinaryHeaderSize);
}

inline void SegyBase::get_trace_header(uchar *traceheader, int n) {
  memcpy(traceheader, trheader(n), kTraceHeaderSize);
}

inline void SegyBase::get_trace_keys(int *dst, const std::vector<int> &keys,
                                     const std::vector<int> &length, int beg,
                                     int end) {

  for (int i = beg; i < end; i++) {
    CHECK_SIGNALS();
    for (int j = 0; j < keys.size(); j++) {
      if (length[j] == 4) {
        *dst = keyi4(i, keys[j]);
      } else if (length[j] == 2) {
        *dst = (int)keyi2(i, keys[j]);
      } else {
        throw std::runtime_error(
            "only spport int32 and int16 type, but got length: " +
            std::to_string(length[j]));
      }
      dst++;
    }
  }
}

inline void SegyBase::itrace(float *data, int n) {
  m_readfunc(data, trDataStart(n), m_meta.nt);
}

inline void SegyBase::collect(float *data, int beg, int end, int tbeg,
                              int tend) {
  int nt = tend - tbeg;
  for (size_t i = beg; i < end; i++) {
    CHECK_SIGNALS();
    m_readfunc(data, trDataStart(i, tbeg), nt);
    data += nt;
  }
}

inline void SegyBase::collect(float *data, const int32_t *index, int n,
                              int tbeg, int tend) {
  int nt = tend - tbeg;
  for (size_t i = 0; i < n; i++) {
    CHECK_SIGNALS();
    // TODO: remove this?
    if (index[i] >= m_meta.ntrace) {
      throw std::runtime_error("Index out of bound." +
                               std::to_string(index[i]));
    }
    if (index[i] < 0) {
      std::fill(data, data + nt, 0);
    } else {
      m_readfunc(data, trDataStart(index[i], tbeg), nt);
    }
    data += nt;
  }
}

/******************************************************/
/**************** For write mode  *********************/
/******************************************************/

inline void SegyBase::set_bkeyi2(int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(bwheader() + loc - 1) =
      swap_endian<int16_t>(val);
}
inline void SegyBase::set_bkeyi4(int loc, int32_t val) {
  *reinterpret_cast<int32_t *>(bwheader() + loc - 1) =
      swap_endian<int32_t>(val);
}

inline void SegyBase::set_keyi2(int n, int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(twheader(n) + loc - 1) =
      swap_endian<int16_t>(val);
}
inline void SegyBase::set_keyi4(int n, int loc, int32_t val) {
  *reinterpret_cast<int32_t *>(twheader(n) + loc - 1) =
      swap_endian<int32_t>(val);
}

inline void SegyBase::write_itrace(const float *data, int n) {
  m_wfunc(twDataStart(n), data, m_meta.nt);
}

inline void SegyBase::write_traces(const float *data, int beg, int end,
                                   int tbeg, int tend) {
  int n = end - beg;
  for (size_t i = beg; i < end; i++) {
    CHECK_SIGNALS();
    m_wfunc(twDataStart(i, tbeg), data + (uint64_t)(i - beg) * n, n);
  }
}

inline void SegyBase::write_traces(const float *data, const int32_t *index,
                                   int n, int tbeg, int tend) {
  int len = tend - tbeg;
  for (size_t i = 0; i < n; i++) {
    CHECK_SIGNALS();
    m_wfunc(twDataStart(index[i], tbeg), data + i * (uint64_t)len, len);
  }
}

/******************************************************/
/**************** For Priviate  ***********************/
/******************************************************/

inline void SegyBase::setRWFunc(int dformat) {
  setRFunc(m_readfunc, dformat);
  setWFunc(m_wfunc, dformat);
}

} // namespace segy
#endif