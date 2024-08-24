#ifndef CIG_SEGY_BASE_H
#define CIG_SEGY_BASE_H
#define FMT_HEADER_ONLY

#include "mio.hpp"
#include "progressbar.hpp"
#include "sutils.h"

#include <fmt/format.h>
#include <stdexcept>
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
  int32_t nt = 0;         // same as time
  int32_t nx = 0;         // same as crossline
  int32_t ni = 0;         // same as inline
  int32_t no = 1;         // same as offset
  int64_t ntrace = 0;     // trace_count
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

  bool isNormalSegy = true;

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

  virtual ~SegyBase() {}

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
  int64_t bkeyi8(int loc);

  // trace header
  int16_t keyi2(int n, int loc);
  int32_t keyi4(int n, int loc);
  int64_t keyi8(int n, int loc);
  inline int32_t iline(int n) { return keyi4(n, m_keys.iline); }
  inline int32_t xline(int n) { return keyi4(n, m_keys.xline); }
  inline int32_t offset(int n) { return keyi4(n, m_keys.offset); }
  inline int32_t coordx(int n) { return keyi4(n, m_keys.xloc); }
  inline int32_t coordy(int n) { return keyi4(n, m_keys.yloc); }

  void get_binary_header(uchar *binheader, bool raw = false);
  void get_trace_header(uchar *traceheader, int n, bool raw = false);
  void get_trace_full(uchar *trace, int n, bool raw = false);
  void get_trace_keys(int *dst, const std::vector<int> &keys,
                      const std::vector<int> &length, int beg, int end);

  void itrace(float *data, int n);
  void collect(float *data, int beg, int end, int tbeg, int tend);
  void collect(float *data, const int32_t *index, int n, int tbeg, int tend);

protected:
  std::string segyname;
  const char *m_data_ptr;
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

  void setRWFunc(int dformat);
};

/** For public **/

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

inline std::string SegyBase::textual_header(char coding) {
  const char *src = m_data_ptr;
  char out[kTextualHeaderSize + kTextualRows];
  bool isEBCDIC = true;
  if (coding == 'a') {
    isEBCDIC = false;
  } else if (coding == 'u') {
    bool isEBCDIC = isTextInEBCDICFormat(src, kTextualHeaderSize);
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
inline int64_t SegyBase::bkeyi8(int loc) {
  return swap_endian<int64_t>(brheader() + loc - 1);
}

inline int16_t SegyBase::keyi2(int n, int loc) {
  return swap_endian<int16_t>(trheader(n) + loc - 1);
}
inline int32_t SegyBase::keyi4(int n, int loc) {
  return swap_endian<int32_t>(trheader(n) + loc - 1);
}
inline int64_t SegyBase::keyi8(int n, int loc) {
  return swap_endian<int64_t>(trheader(n) + loc - 1);
}

inline void SegyBase::get_binary_header(uchar *binheader, bool raw) {
  const char *src = m_data_ptr + kTextualHeaderSize;
  if (raw) {
    memcpy(binheader, src, kBinaryHeaderSize);
  } else {
    read_binary_header(binheader, src);
  }
}

inline void SegyBase::get_trace_header(uchar *traceheader, int n, bool raw) {
  if (raw) {
    memcpy(traceheader, trheader(n), kTraceHeaderSize);
  } else {
    read_one_trace_header(traceheader, trheader(n));
  }
}

inline void SegyBase::get_trace_full(uchar *trace, int n, bool raw) {
  const char *src = trheader(n);
  if (raw) {
    memcpy(trace, src, m_meta.tracesize);
  } else {
    read_one_trace_header(trace, src);
    float *data = reinterpret_cast<float *>(trace + kTraceHeaderSize);
    m_readfunc(data, src + kTraceHeaderSize, m_meta.nt);
  }
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
        throw std::runtime_error("only spport int32 and int16 type.");
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
  int n = end - beg;
  for (size_t i = beg; i < end; i++) {
    CHECK_SIGNALS();
    m_readfunc(data, trDataStart(i) + tbeg, n);
    data += n;
  }
}

inline void SegyBase::collect(float *data, const int32_t *index, int n,
                              int tbeg, int tend) {
  int len = tend - tbeg;
  for (size_t i = 0; i < n; i++) {
    CHECK_SIGNALS();
    m_readfunc(data + i * (uint64_t)len, trDataStart(index[i]) + tbeg, len);
  }
}

inline void SegyBase::setRWFunc(int dformat) {
  switch (dformat) {
  case 1:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npibm(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyibm(dst, src, size);
    };
    break;
  case 2:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<int32_t>(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<int32_t>(dst, src, size);
    };
    break;
  case 3:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<int16_t>(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<int16_t>(dst, src, size);
    };
    break;
  case 5:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<float>(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<float>(dst, src, size);
    };
    break;
  case 8:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<int8_t>(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<int8_t>(dst, src, size);
    };
    break;
  case 10:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<uint32_t>(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<uint32_t>(dst, src, size);
    };
    break;
  case 11:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<uint16_t>(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<uint16_t>(dst, src, size);
    };
    break;
  case 16:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<uint8_t>(dst, src, size);
    };
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<uint8_t>(dst, src, size);
    };
    break;
  default:
    throw std::invalid_argument("Unsupported dformat value");
  }
}

} // namespace segy
#endif