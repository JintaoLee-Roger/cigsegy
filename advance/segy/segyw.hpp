#ifndef CIG_SEGY_WRITE_H
#define CIG_SEGY_WRITE_H
#include "segybase.hpp"
#include "segynd.hpp"
#include "sutils.h"
#include <vector>

namespace segy {

class SegyRW : public SegyND {
public:
  explicit SegyRW(const std::string &segyname);
  ~SegyRW() override;

  inline char *bwheader() { return m_sink.data() + kTextualHeaderSize; }
  inline char *twheader(int n) {
    return m_sink.data() + kTraceHeaderStart + n * m_meta.tracesize;
  }
  inline char *twDataStart(int n) { return twheader(n) + kTraceHeaderSize; }

  void set_bkeyi2(int loc, int16_t val);
  void set_bkeyi4(int loc, int32_t val);
  void set_bkeyi8(int loc, int64_t val);

  void set_keyi2(int n, int loc, int16_t val);
  void set_keyi4(int n, int loc, int32_t val);
  void set_keyi8(int n, int loc, int64_t val);

  inline void set_iline(int n, int32_t val) { set_keyi4(n, m_keys.iline, val); }
  inline void set_xline(int n, int32_t val) { set_keyi4(n, m_keys.xline, val); }
  inline void set_offset(int n, int32_t val) {
    set_keyi4(n, m_keys.offset, val);
  }
  inline void set_coordx(int n, int32_t val) { set_keyi4(n, m_keys.xloc, val); }
  inline void set_coordy(int n, int32_t val) { set_keyi4(n, m_keys.yloc, val); }

  // void set_binary_header(const char *binheader, bool raw = false);
  // void set_trace_header(const char *traceheader, int n, bool raw = false);
  // void set_trace_full(const char *trace, int n, bool raw = false);
  void set_trace_keys(const int *dst, const std::vector<int> &keys,
                      const std::vector<int> &length, int beg, int end);

  void write_itrace(const float *data, int n);
  void write_traces(const float *data, int beg, int end, int tbeg, int tend);
  void write_traces(const float *data, const int32_t *index, int n, int tbeg,
                    int tend);

  void write(const float *data);
  void write2d(const float *data, int is, int ie, int ts, int te);
  void write3d(const float *data, int is, int ie, int xs, int xe, int ts,
               int te);
  void write4d(const float *data, int is, int ie, int xs, int xe, int os,
               int oe, int ts, int te);

protected:
  mio::mmap_sink m_sink;
};

inline SegyRW::SegyRW(const std::string &segyname) : SegyND(segyname) {
  std::error_code error;
  this->m_sink.map(segyname, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
}

inline SegyRW::~SegyRW() {
  if (m_src.is_mapped()) {
    m_src.unmap();
  }
  if (m_sink.is_mapped()) {
    m_sink.unmap();
  }
}

inline void SegyRW::set_bkeyi2(int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(bwheader() + loc - 1) =
      swap_endian<int16_t>(val);
}
inline void SegyRW::set_bkeyi4(int loc, int32_t val) {
  *reinterpret_cast<int32_t *>(bwheader() + loc - 1) =
      swap_endian<int32_t>(val);
}
inline void SegyRW::set_bkeyi8(int loc, int64_t val) {
  *reinterpret_cast<int64_t *>(bwheader() + loc - 1) =
      swap_endian<int64_t>(val);
}

inline void SegyRW::set_keyi2(int n, int loc, int16_t val) {
  *reinterpret_cast<int16_t *>(twheader(n) + loc - 1) =
      swap_endian<int16_t>(val);
}
inline void SegyRW::set_keyi4(int n, int loc, int32_t val) {
  *reinterpret_cast<int32_t *>(twheader(n) + loc - 1) =
      swap_endian<int32_t>(val);
}
inline void SegyRW::set_keyi8(int n, int loc, int64_t val) {
  *reinterpret_cast<int64_t *>(twheader(n) + loc - 1) =
      swap_endian<int64_t>(val);
}

} // namespace segy

#endif