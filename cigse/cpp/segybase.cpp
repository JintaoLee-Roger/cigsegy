/*********************************************************************
** Copyright (c) 2024 Jintao Li.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
** All rights reserved.
*********************************************************************/

#include "segybase.h"

namespace segy {

/******************************************************/
/**************** For Public  *************************/
/******************************************************/

/**************** set locations  *************************/

void SegyBase::setLocations(int iloc, int xloc, int oloc) {
  m_keys.iline = iloc;
  m_keys.xline = xloc;
  m_keys.offset = oloc;
}
void SegyBase::setInlineLocation(int loc) { m_keys.iline = loc; }
void SegyBase::setCrosslineLocation(int loc) { m_keys.xline = loc; }
void SegyBase::setOffsetLocation(int loc) { m_keys.offset = loc; }

void SegyBase::setXLocation(int loc) { m_keys.xloc = loc; }
void SegyBase::setYLocation(int loc) { m_keys.yloc = loc; }
void SegyBase::setXYLocations(int xloc, int yloc) {
  m_keys.xloc = xloc;
  m_keys.yloc = yloc;
}

void SegyBase::setInlineStep(int step) { m_keys.istep = step; }
void SegyBase::setCrosslineStep(int step) { m_keys.xstep = step; }
void SegyBase::setOffsetStep(int step) { m_keys.ostep = step; }
void SegyBase::setSteps(int istep, int xstep, int ostep) {
  m_keys.istep = istep;
  m_keys.xstep = xstep;
  m_keys.ostep = ostep;
}

/******************************************************/
/**************** For read mode  **********************/
/******************************************************/

std::string SegyBase::textual_header(char coding) {
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

int16_t SegyBase::bkeyi2(int loc) {
  return swap_endian<int16_t>(brheader() + loc - 1);
}
int32_t SegyBase::bkeyi4(int loc) {
  return swap_endian<int32_t>(brheader() + loc - 1);
}

int16_t SegyBase::keyi2(int n, int loc) {
  return swap_endian<int16_t>(trheader(n) + loc - 1);
}
int32_t SegyBase::keyi4(int n, int loc) {
  return swap_endian<int32_t>(trheader(n) + loc - 1);
}

void SegyBase::get_binary_header(uchar *binheader) {
  memcpy(binheader, brheader(), kBinaryHeaderSize);
}

void SegyBase::get_trace_header(uchar *traceheader, int n) {
  memcpy(traceheader, trheader(n), kTraceHeaderSize);
}

void SegyBase::get_trace_keys(int *dst, const std::vector<int> &keys,
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

void SegyBase::itrace(float *data, int n) {
  m_readfunc(data, trDataStart(n), m_meta.nt);
}

void SegyBase::collect(float *data, int beg, int end, int tbeg, int tend) {
  int nt = tend - tbeg;
  for (size_t i = beg; i < end; i++) {
    CHECK_SIGNALS();
    m_readfunc(data, trDataStart(i) + tbeg, nt);
    data += nt;
  }
}

void SegyBase::collect(float *data, const int32_t *index, int n, int tbeg,
                       int tend) {
  int nt = tend - tbeg;
  for (size_t i = 0; i < n; i++) {
    CHECK_SIGNALS();
    m_readfunc(data + i * (uint64_t)nt, trDataStart(index[i]) + tbeg, nt);
  }
}

/******************************************************/
/**************** For write mode  *********************/
/******************************************************/

void SegyBase::set_bkeyi2(int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(bwheader() + loc - 1) =
      swap_endian<int16_t>(val);
}
void SegyBase::set_bkeyi4(int loc, int32_t val) {
  *reinterpret_cast<int32_t *>(bwheader() + loc - 1) =
      swap_endian<int32_t>(val);
}

void SegyBase::set_keyi2(int n, int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(twheader(n) + loc - 1) =
      swap_endian<int16_t>(val);
}
void SegyBase::set_keyi4(int n, int loc, int32_t val) {
  *reinterpret_cast<int32_t *>(twheader(n) + loc - 1) =
      swap_endian<int32_t>(val);
}

void SegyBase::write_itrace(const float *data, int n) {
  m_wfunc(twDataStart(n), data, m_meta.nt);
}

void SegyBase::write_traces(const float *data, int beg, int end, int tbeg,
                            int tend) {
  int n = end - beg;
  for (size_t i = beg; i < end; i++) {
    CHECK_SIGNALS();
    m_wfunc(twDataStart(i) + tbeg, data + (uint64_t)(i - beg) * n, n);
  }
}

void SegyBase::write_traces(const float *data, const int32_t *index, int n,
                            int tbeg, int tend) {
  int len = tend - tbeg;
  for (size_t i = 0; i < n; i++) {
    CHECK_SIGNALS();
    m_wfunc(twDataStart(index[i]) + tbeg, data + i * (uint64_t)len, len);
  }
}

/******************************************************/
/**************** For Priviate  ***********************/
/******************************************************/

void SegyBase::setRWFunc(int dformat) {
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
