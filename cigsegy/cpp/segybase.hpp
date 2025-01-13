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
  size_t iline = 189; // inline number location
  size_t xline = 193; // crossline number location
  size_t offset = 37; // offset number location
  size_t xloc = 181;  // X location
  size_t yloc = 185;  // Y location
  int istep = 1;      // inline step
  int xstep = 1;      // crossline step
  int ostep = 1;      // offset step
};

struct MetaInfo {
  // count information
  size_t nt = 0;          // number of samples in time dimension
  size_t no = 1;          // number of samples in offset dimension
  size_t nx = 0;          // number of samples in crosline dimension
  size_t ni = 0;          // number of samples in inline dimension
  uint64_t ntrace = 0;    // number of traces
  uint64_t tracesize = 0; // trace size in char

  size_t dt = 0; // interval of time/depth
  float dx = 0;  // interval of crossline
  float di = 0;  // interval of inline

  int scalar = 0;
  int dformat = 0;

  size_t start_time = 0;
  size_t start_iline = 0;
  size_t end_iline = 0;
  size_t start_xline = 0;
  size_t end_xline = 0;
  size_t start_offset = 0;
  size_t end_offset = 0;

  int trace_sorting_code = 0;
  size_t esize = 4;
  float fillNoValue = 0;
};

class SegyBase {
public:
  SegyBase() : m_data_ptr(nullptr) {}

  ~SegyBase() { this->close_file(); }

  KeyLocs m_keys;
  MetaInfo m_meta;

  void close_file() {
    if (m_sink.is_mapped()) {
      m_sink.unmap();
    }
    if (m_src.is_mapped()) {
      m_src.unmap();
    }
  }

  inline uint64_t ntrace() const { return m_meta.ntrace; }
  inline size_t nt() const { return m_meta.nt; }

  void setLocations(size_t iline, size_t xline, size_t offset = 37);
  void setInlineLocation(size_t loc);
  void setCrosslineLocation(size_t loc);
  void setOffsetLocation(size_t loc);

  void setXLocation(size_t loc);
  void setYLocation(size_t loc);
  void setXYLocations(size_t xloc, size_t yloc);

  void setInlineStep(int step);
  void setCrosslineStep(int step);
  void setOffsetStep(int step);
  void setSteps(int istep, int xstep, int ostep = 1);
  void setFill(float fill) { m_meta.fillNoValue = fill; }

  std::string textual_header(char coding = 'u');

  void show_progress(bool show) { show_progress_ = show; }

  // binary header
  int16_t bkeyi2(size_t loc);
  int32_t bkeyi4(size_t loc);

  // trace header
  int16_t keyi2(size_t n, size_t loc);
  int32_t keyi4(size_t n, size_t loc);
  inline int32_t iline(size_t n) { return keyi4(n, m_keys.iline); }
  inline int32_t xline(size_t n) { return keyi4(n, m_keys.xline); }
  inline int32_t offset(size_t n) { return keyi4(n, m_keys.offset); }
  inline int32_t coordx(size_t n) { return keyi4(n, m_keys.xloc); }
  inline int32_t coordy(size_t n) { return keyi4(n, m_keys.yloc); }

  void get_binary_header(uchar *binheader);
  void get_trace_header(uchar *traceheader, size_t n);

  void get_trace_keys(int32_t *dst, const std::vector<size_t> &keys,
                      const std::vector<size_t> &length, size_t beg,
                      size_t end);

  void itrace(float *data, size_t n);
  void collect(float *data, size_t beg, size_t end, size_t tbeg, size_t tend);
  void collect(float *data, const int32_t *index, size_t n, size_t tbeg,
               size_t tend);

  // W mode
  void set_bkeyi2(size_t loc, int16_t val);
  void set_bkeyi4(size_t loc, int32_t val);

  void set_keyi2(size_t n, size_t loc, int16_t val);
  void set_keyi4(size_t n, size_t loc, int32_t val);

  inline void set_iline(size_t n, int32_t val) {
    set_keyi4(n, m_keys.iline, val);
  }
  inline void set_xline(size_t n, int32_t val) {
    set_keyi4(n, m_keys.xline, val);
  }
  inline void set_offset(size_t n, int32_t val) {
    set_keyi4(n, m_keys.offset, val);
  }
  inline void set_coordx(size_t n, int32_t val) {
    set_keyi4(n, m_keys.xloc, val);
  }
  inline void set_coordy(size_t n, int32_t val) {
    set_keyi4(n, m_keys.yloc, val);
  }

  // void set_trace_keys(const int *dst, const std::vector<int> &keys,
  //                     const std::vector<int> &length, int beg, int end);

  void write_itrace(const float *data, size_t n);
  void write_traces(const float *data, size_t beg, size_t end, size_t tbeg,
                    size_t tend);
  void write_traces(const float *data, const int32_t *index, size_t n,
                    size_t tbeg, size_t tend);

protected:
  const char *m_data_ptr;
  mio::mmap_sink m_sink;
  mio::mmap_source m_src;
  ReadFunc m_readfunc;
  ReadFuncOne m_readfuncone;
  WriteFunc m_wfunc;
  bool m_w = true;
  bool show_progress_ = true;

  inline const char *brheader() const {
    return m_data_ptr + kTextualHeaderSize;
  }
  inline const char *trheader(size_t n) const {
    return m_data_ptr + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline const char *trDataStart(size_t n, size_t tbeg = 0) const {
    return trheader(n) + kTraceHeaderSize + tbeg * m_meta.esize;
  }

  inline char *bwheader() { return m_sink.data() + kTextualHeaderSize; }
  inline char *twheader(size_t n) {
    return m_sink.data() + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline char *twDataStart(size_t n, size_t tbeg = 0) {
    return twheader(n) + kTraceHeaderSize + tbeg * m_meta.esize;
  }

  void setRWFunc(int dformat);
};

/******************************************************/
/**************** For Public  *************************/
/******************************************************/

/**************** set locations  *************************/

inline void SegyBase::setLocations(size_t iloc, size_t xloc, size_t oloc) {
  m_keys.iline = iloc;
  m_keys.xline = xloc;
  m_keys.offset = oloc;
}
inline void SegyBase::setInlineLocation(size_t loc) { m_keys.iline = loc; }
inline void SegyBase::setCrosslineLocation(size_t loc) { m_keys.xline = loc; }
inline void SegyBase::setOffsetLocation(size_t loc) { m_keys.offset = loc; }

inline void SegyBase::setXLocation(size_t loc) { m_keys.xloc = loc; }
inline void SegyBase::setYLocation(size_t loc) { m_keys.yloc = loc; }
inline void SegyBase::setXYLocations(size_t xloc, size_t yloc) {
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
  for (size_t iRow = 0; iRow < kTextualRows; iRow++) {
    size_t offset = iRow * kTextualColumns;
    for (size_t iCol = 0; iCol < kTextualColumns; iCol++) {
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

inline int16_t SegyBase::bkeyi2(size_t loc) {
  return swap_endian<int16_t>(brheader() + loc - 1);
}
inline int32_t SegyBase::bkeyi4(size_t loc) {
  return swap_endian<int32_t>(brheader() + loc - 1);
}

inline int16_t SegyBase::keyi2(size_t n, size_t loc) {
  return swap_endian<int16_t>(trheader(n) + loc - 1);
}
inline int32_t SegyBase::keyi4(size_t n, size_t loc) {
  return swap_endian<int32_t>(trheader(n) + loc - 1);
}

inline void SegyBase::get_binary_header(uchar *binheader) {
  memcpy(binheader, brheader(), kBinaryHeaderSize);
}

inline void SegyBase::get_trace_header(uchar *traceheader, size_t n) {
  memcpy(traceheader, trheader(n), kTraceHeaderSize);
}

inline void SegyBase::get_trace_keys(int32_t *dst,
                                     const std::vector<size_t> &keys,
                                     const std::vector<size_t> &length,
                                     size_t beg, size_t end) {

  int step = cal_progress_steps(end - beg, show_progress_, 10000, 50);
  for (size_t i = beg; i < end; i++) {
    g_check_signals_callback();
    if (step && (i-beg) % step == 0) {
      g_progress_callback(i-beg, end-beg);
    }
    for (size_t j = 0; j < keys.size(); j++) {
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
  if (step) {
    g_progress_callback(-1, end-beg);
  }
}

inline void SegyBase::itrace(float *data, size_t n) {
  m_readfunc(data, trDataStart(n), m_meta.nt);
}

inline void SegyBase::collect(float *data, size_t beg, size_t end, size_t tbeg,
                              size_t tend) {
  size_t nt = tend - tbeg;
  int step = cal_progress_steps(end - beg, show_progress_, 10000, 50);
  for (size_t i = beg; i < end; i++) {
    g_check_signals_callback();
    if (step && (i-beg) % step == 0) {
      g_progress_callback(i-beg, end-beg);
    }
    m_readfunc(data, trDataStart(i, tbeg), nt);
    data += nt;
  }
  if (step) {
    g_progress_callback(-1, end-beg);
  }
}

inline void SegyBase::collect(float *data, const int32_t *index, size_t n,
                              size_t tbeg, size_t tend) {
  size_t nt = tend - tbeg;
  int step = cal_progress_steps(n, show_progress_, 10000, 50);
  for (size_t i = 0; i < n; i++) {
    g_check_signals_callback();
    if (step && i % step == 0) {
      g_progress_callback(i, n);
    }
    // TODO: remove this?
    if (index[i] >= static_cast<int32_t>(m_meta.ntrace)) {
      throw std::runtime_error("Index out of bound. Index: " +
                               std::to_string(index[i]));
    }
    if (index[i] < 0) {
      std::fill(data, data + nt, m_meta.fillNoValue);
    } else {
      m_readfunc(data, trDataStart(index[i], tbeg), nt);
    }
    data += nt;
  }
  if (step) {
    g_progress_callback(-1, n);
  }
}

/******************************************************/
/**************** For write mode  *********************/
/******************************************************/

inline void check_write(bool w){
  if (!w) {
    throw std::runtime_error("You set write=false, so you can't access write functions.");
  }
}

inline void SegyBase::set_bkeyi2(size_t loc, int16_t val) {
  check_write(m_w);
  *reinterpret_cast<int16_t *>(bwheader() + loc - 1) =
      swap_endian<int16_t>(val);
}
inline void SegyBase::set_bkeyi4(size_t loc, int32_t val) {
  check_write(m_w);
  *reinterpret_cast<int32_t *>(bwheader() + loc - 1) =
      swap_endian<int32_t>(val);
}

inline void SegyBase::set_keyi2(size_t n, size_t loc, int16_t val) {
  check_write(m_w);
  *reinterpret_cast<int16_t *>(twheader(n) + loc - 1) =
      swap_endian<int16_t>(val);
}
inline void SegyBase::set_keyi4(size_t n, size_t loc, int32_t val) {
  check_write(m_w);
  *reinterpret_cast<int32_t *>(twheader(n) + loc - 1) =
      swap_endian<int32_t>(val);
}

inline void SegyBase::write_itrace(const float *data, size_t n) {
  check_write(m_w);
  m_wfunc(twDataStart(n), data, m_meta.nt);
}

inline void SegyBase::write_traces(const float *data, size_t beg, size_t end,
                                   size_t tbeg, size_t tend) {
  check_write(m_w);
  size_t n = end - beg;
  int step = cal_progress_steps(n, show_progress_, 10000, 50);
  for (size_t i = beg; i < end; i++) {
    g_check_signals_callback();
    if (step && (i-beg) % step == 0) {
      g_progress_callback(i-beg, n);
    }
    m_wfunc(twDataStart(i, tbeg), data + (uint64_t)(i - beg) * n, n);
  }
  if (step) {
    g_progress_callback(-1, n);
  }
}

inline void SegyBase::write_traces(const float *data, const int32_t *index,
                                   size_t n, size_t tbeg, size_t tend) {
  check_write(m_w);
  size_t len = tend - tbeg;
  int step = cal_progress_steps(n, show_progress_, 10000, 50);
  for (size_t i = 0; i < n; i++) {
    g_check_signals_callback();
    if (step && i % step == 0) {
      g_progress_callback(i, n);
    }
    if (index[i] >= static_cast<int32_t>(m_meta.ntrace)) {
      throw std::runtime_error("Index out of bound." +
                               std::to_string(index[i]));
    }
    if (index[i] < 0) {
      continue;
    }
    m_wfunc(twDataStart(index[i], tbeg), data + i * (uint64_t)len, len);
  }
  if (step) {
    g_progress_callback(-1, n);
  }
}

/******************************************************/
/**************** For Priviate  ***********************/
/******************************************************/

inline void SegyBase::setRWFunc(int dformat) {
  setRFunc(m_readfunc, dformat);
  setWFunc(m_wfunc, dformat);
  setRFuncOne(m_readfuncone, dformat);
}

} // namespace segy
#endif