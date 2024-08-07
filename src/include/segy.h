/*********************************************************************
** Copyright (c) 2023 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
**
** @File: segy.h
** @Description :
*********************************************************************/

#ifndef CIG_SEGY_H
#define CIG_SEGY_H
#define FMT_HEADER_ONLY

#include <stdexcept>
#include <vector>
// #include <omp.h>
#include "mio.hpp"
#include "utils.h"
#include <fmt/format.h>

#ifdef USE_PYBIND11
#include <pybind11/pybind11.h>

inline void checkSignals() {
    if (PyErr_CheckSignals() != 0) {
        throw pybind11::error_already_set();
    }
}
#define CHECK_SIGNALS() checkSignals()

#else
#define CHECK_SIGNALS() do {} while (0)
#endif

namespace segy {

extern bool showpbar;

struct MetaInfo {
  // count information
  int32_t sizeX; // same as time
  int32_t sizeY; // same as crossline
  int32_t sizeZ; // same as inline
  int64_t trace_count;
  int16_t sample_interval; // dt
  int16_t data_format;     // 1 or 5
  float Y_interval;        // interval of crossline
  float Z_interval;        // interval of inline
  int16_t start_time;
  int16_t scalar;

  int min_inline;
  int max_inline;
  int min_crossline;
  int max_crossline;

  bool isNormalSegy;

  float fillNoValue;

  // field information
  int inline_field = kDefaultInlineField;
  int crossline_field = kDefaultCrosslineField;
  int X_field = kDefaultXField;
  int Y_field = kDefaultYField;

  int inline_step = 1;
  int crossline_step = 1;

  int trace_sorting_code;
  int esize=4;
};

struct LineInfo {
  int line_num;
  int crossline_start;
  int crossline_end;
  int trace_start;
  int trace_end;
  int count;
};

struct TraceInfo {
  int inline_num;
  int crossline_num;
  int X;
  int Y;
};

class SegyIO {
public:
  // read segy mode
  explicit SegyIO(const std::string &segyname);
  // create segy from memory
  SegyIO(int sizeX, int sizeY, int sizeZ);
  // create segy file from binary file
  SegyIO(const std::string &binaryname, int sizeX, int sizeY, int sizeZ);

  ~SegyIO();

  inline int shape(int dimension) {
    if (dimension == 0) {
      return m_metaInfo.sizeX;
    } else if (dimension == 1) {
      return m_metaInfo.sizeY;
    } else if (dimension == 2) {
      return m_metaInfo.sizeZ;
    } else {
      throw std::runtime_error("shape(dim), dim can be only {0, 1, 2}");
    }
  }

  int64_t trace_count();
  int nt();
  void set_size(int x, int y, int z);
  MetaInfo get_metaInfo();

  inline std::vector<LineInfo> line_info() { return m_lineInfo; }

  // read segy
  void setInlineLocation(int loc);
  void setCrosslineLocation(int loc);
  void setXLocation(int loc);
  void setYLocation(int loc);
  void setInlineStep(int step);
  void setCrosslineStep(int step);
  void setSteps(int istep, int xstep);
  void setFillNoValue(float noValue);

  // create segy
  void setSampleInterval(int interval);
  void setDataFormatCode(int fdormat);
  void setStartTime(int start_time);
  void setInlineInterval(float dz);
  void setCrosslineInterval(float dy);
  void setMinInline(int in);
  void setMinCrossline(int cross);

  // 'u' is unknown (default), 'e' is EBCDIC, 'a' is ASCII
  std::string textual_header(char coding = 'u');
  std::string metaInfo();

  void get_binary_header_full(uchar *binheader, bool raw = false);
  void get_trace_header_full(int n, uchar *traceheader, bool raw = false);
  void get_trace_full(int n, uchar *trace, bool raw = false);
  void get_trace_keys_c(int *dst, const std::vector<int>& keys, const std::vector<int>& length, int beg, int end);

  void collect(float *data, int beg = -1, int end = 0);

  // read segy
  void scan();
  void read(float *dst, int startX, int endX, int startY, int endY, int startZ,
            int endZ);
  void read(float *dst);
  void read(float *dst, int sizeY, int sizeZ, int minY, int minZ); // for read unsorted
  void read_inline_slice(float *dst, int iZ);
  void read_cross_slice(float *dst, int iY);
  void read_time_slice(float *dst, int iX);
  void read_trace(float *dst, int iY, int iZ);
  void tofile(const std::string &binary_out_name);
  void cut(const std::string &outname, int startX, int endX, int startY, int endY,
      int startZ, int endZ,
      const std::vector<std::string> &custom_info = std::vector<std::string>());
  void cut(const std::string &outname, int startY, int endY, int startZ, int endZ,
      const std::vector<std::string> &custom_info = std::vector<std::string>());
  void cut(const std::string &outname, int startX, int endX,
      const std::vector<std::string> &custom_info = std::vector<std::string>());

  // create
  void create(
      const std::string &segy_out_name, const float *src,
      const std::vector<std::string> &custom_info = std::vector<std::string>());
  void create(
      const std::string &segy_out_name,
      const std::vector<std::string> &custom_info = std::vector<std::string>());

  void close_file();

private:
  bool isReadSegy;
  bool isScan = false;
  mio::mmap_source m_source;
  mio::mmap_sink m_sink;
  std::vector<LineInfo> m_lineInfo;
  MetaInfo m_metaInfo;

  void scanBinaryHeader();
  void initMetaInfo();
  void initTraceHeader(TraceHeader *trace_header);
  void write_textual_header(
      char *dst, const std::string &segy_out_name,
      const std::vector<std::string> &custom_info = std::vector<std::string>());
  void write_binary_header(char *dst);
  void write_trace_header(char *dst, TraceHeader *trace_header, int32_t iY,
                          int32_t iZ, int32_t x, int32_t y);
  void read_all_fast(float *dst);

  inline void _get_TraceInfo(uint64_t n, TraceInfo &tmetaInfo) {
    const char *field =
        m_source.data() + kTextualHeaderSize + kBinaryHeaderSize +
        n * (kTraceHeaderSize + m_metaInfo.sizeX * m_metaInfo.esize);
    tmetaInfo.inline_num =
        swap_endian<int32_t>(field + m_metaInfo.inline_field - 1);
    tmetaInfo.crossline_num =
        swap_endian<int32_t>(field + m_metaInfo.crossline_field - 1);
    tmetaInfo.X = swap_endian<int32_t>(field + m_metaInfo.X_field - 1);
    tmetaInfo.Y = swap_endian<int32_t>(field + m_metaInfo.Y_field - 1);
  }
};

inline void disable_progressbar() { showpbar = false; }

void read(const std::string &segy_name, float *dst,
          int iline = kDefaultInlineField, int xline = kDefaultCrosslineField,
          int istep = 1, int xstep = 1);

void tofile(const std::string &segy_name, const std::string &out_name,
            int iline = kDefaultInlineField, int xline = kDefaultCrosslineField,
            int istep = 1, int xstep = 1);

void create_by_sharing_header(
    const std::string &segy_name, const std::string &header_segy,
    const float *src, int sizeX, int sizeY, int sizeZ, int iline = 189,
    int xline = 193, int istep = 1, int xstep = 1, int offsetX = -1,
    int offsetY = -1, int offsetZ = -1,
    const std::vector<std::string> &custom_info = std::vector<std::string>());

void create_by_sharing_header(
    const std::string &segy_name, const std::string &header_segy,
    const std::string &src_file, int sizeX, int sizeY, int sizeZ,
    int iline = 189, int xline = 193, int istep = 1, int xstep = 1,
    int offsetX = -1, int offsetY = -1, int offsetZ = -1,
    const std::vector<std::string> &custom_info = std::vector<std::string>());

void load_prestack3D(float *dst, const std::string &segy_name, int sizeX,
                     int min_inline, int max_inline, int min_crossline,
                     int max_crossline, int min_offset, int max_offset,
                     int inline_step, int crossline_step, int offset_step,
                     int inline_field, int crossline_field, int offset_field,
                     float fill);

template <typename T>
void modify_bin_key(const std::string &segy_name, int loc, T value) {
  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(segy_name, error);
  if (error) {
    throw std::runtime_error("mmap fail when write data");
  }

  T *data =
      reinterpret_cast<T *>(rw_mmap.data() + kTextualHeaderSize + loc - 1);
  *data = swap_endian(value);
  rw_mmap.unmap();
}

template <typename T>
void modify_trace_key(const std::string &segy_name, int loc, T value, int idx) {
  SegyIO segy_data(segy_name);
  int sizeX = segy_data.nt();
  int esize = segy_data.get_metaInfo().esize;
  int tracecount = segy_data.trace_count();
  int tracesize = sizeX * esize + kTraceHeaderSize;
  segy_data.close_file();
  if (idx >= tracecount) {
    throw std::runtime_error(
        fmt::format("trace index out of range. (index={}, trace_count={})", idx,
                    tracecount));
  }

  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(segy_name, error);
  if (error) {
    throw std::runtime_error("mmap fail when write data");
  }

  if (idx >= 0) {
    T *data =
        reinterpret_cast<T *>(rw_mmap.data() + kTextualHeaderSize +
                              kBinaryHeaderSize + idx * tracesize + loc - 1);
    *data = swap_endian(value);
  } else {
    for (int i = 0; i < tracecount; i++) {
      T *data =
          reinterpret_cast<T *>(rw_mmap.data() + kTextualHeaderSize +
                                kBinaryHeaderSize + i * tracesize + loc - 1);
      *data = swap_endian(value);
    }
  }
  rw_mmap.unmap();
}

} // namespace segy

#endif