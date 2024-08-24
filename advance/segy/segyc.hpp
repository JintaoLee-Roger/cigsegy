#ifndef CIG_SEGY_CREATE_H
#define CIG_SEGY_CREATE_H
#include "segybase.hpp"
#include "sutils.h"
#include <vector>

namespace segy {

class SegyC : public SegyBase {
public:
  explicit SegyC(const std::string &segyname, int ntrace, int nt);
  explicit SegyC(const std::string &segyname, int ni, int nx, int nt);
  explicit SegyC(const std::string &segyname, int ni, int nx, int no, int nt);
  ~SegyC() override;

  void setSampleInterval(int interval) { m_meta.dt = interval; }
  void setDataFormatCode(int dformat) {
    m_meta.dformat = dformat;
    setRWFunc(dformat);
  }
  void setStartTime(int start_time) { m_meta.start_time = start_time; };
  void setInlineInterval(float di) { m_meta.di = di; }
  void setCrosslineInterval(float dx) { m_meta.dx = dx; }
  void setMinInline(int il) { m_meta.start_iline = il; }
  void setMinCrossline(int xl) { m_meta.start_offset = xl; }
  void setMinOffset(int of) { m_meta.start_offset = of; }

  void copy_textual_from(const std::string &segyname);
  void copy_bheader_from(const std::string &segyname);
  void copy_theader_from(const std::string &segyname);

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

  void add_trace(const float *src, int idx, int iline, int xline,
                 int offset = 0);
  void create(const float *src, const int32_t *xyico);

protected:
  mio::mmap_sink m_sink;

private:
  int m_ndim;
  std::string textual;
  std::vector<char> th = std::vector<char>(kTraceHeaderSize, 0);
  uint64_t _need_size(int ni, int nx, int no, int nt);
  void init_info(int ni, int nx, int no, int nt);
  void update_bheader();
  void update_theader();
};

inline SegyC::SegyC(const std::string &segyname, int ntrace, int nt)
    : SegyBase(segyname) {
  uint64_t needsize = _need_size(ntrace, 1, 1, nt);
  create_file(segyname, needsize);
  std::error_code error;
  this->m_sink.map(segyname, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
  m_data_ptr = this->m_sink.data();
  m_ndim = 2;
  init_info(ntrace, 1, 1, nt);
}

inline SegyC::SegyC(const std::string &segyname, int ni, int nx, int nt)
    : SegyBase(segyname) {
  uint64_t needsize = _need_size(ni, nx, 1, nt);
  create_file(segyname, needsize);
  std::error_code error;
  this->m_sink.map(segyname, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
  m_data_ptr = this->m_sink.data();
  m_ndim = 3;
  init_info(ni, nx, 1, nt);
}

inline SegyC::SegyC(const std::string &segyname, int ni, int nx, int no, int nt)
    : SegyBase(segyname) {
  uint64_t needsize = _need_size(ni, nx, no, nt);
  create_file(segyname, needsize);
  std::error_code error;
  this->m_sink.map(segyname, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
  m_data_ptr = this->m_sink.data();
  m_ndim = 4;
  init_info(ni, nx, no, nt);
}

inline SegyC::~SegyC() {
  if (m_sink.is_mapped()) {
    m_sink.unmap();
  }
}

inline uint64_t SegyC::_need_size(int ni, int nx, int no, int nt) {
  return kTraceHeaderStart +
         static_cast<uint64_t>(ni) * nx * no * (nt * 4 + kTraceHeaderSize);
}

inline void SegyC::init_info(int ni, int nx, int no, int nt) {
  m_meta.ni = ni;
  m_meta.nx = nx;
  m_meta.no = no;
  m_meta.nt = nt;
  m_meta.ntrace = static_cast<int64_t>(ni) * nx * no;
  m_meta.dt = 2000;
  m_meta.dx = 25;
  m_meta.di = 25;
  m_meta.scalar = 10;
  m_meta.dformat = 5;
  setRWFunc(5);
  m_meta.start_iline = 1000;
  m_meta.start_xline = 2000;
  m_meta.start_offset = 100;

  if (m_ndim == 2) {
    m_meta.trace_sorting_code = 0;
  } else if (m_ndim == 3) {
    m_meta.trace_sorting_code = 4;
  } else {
    m_meta.trace_sorting_code = 2;
  }
}

inline void SegyC::update_bheader() {
  set_bkeyi2(kBSampleCountField, (int16_t)m_meta.nt);
  set_bkeyi2(kBSampleIntervalField, (int16_t)m_meta.dt);
}

} // namespace segy

#endif