/*********************************************************************
** Copyright (c) 2023 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
**
** @File: segy.cpp
** @Description :
*********************************************************************/

#include "segy.h"
#include <cassert>
#include <chrono>
#include <cstring>
#include <numeric>
#include <vector>
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fstream>
#include <stdexcept>

// #ifdef WIN32
// #include <fcntl.h>
// #include <io.h>
// #define open _open
// #define lseek _lseek
// #define close _close
// #define write _write
// #define O_RDWR _O_RDWR
// #define O_CREAT _O_CREAT
// #define O_TRUNC _O_TRUNC
// #endif

#include "mio.hpp"
#include "progressbar.hpp"
#include "utils.h"

namespace segy {

/******************** extern ******************/
bool showpbar = true;

/********************* static functions ***************************/

static inline void updatebar(progressbar &bar, int n) {
  for (int i = 0; i < n; i++) {
    bar.update();
  }
}

static inline std::string field_str(int field, int len = 2) {
  return fmt::format("{}-{}", field, field + len - 1);
}

static inline int32_t getCrossline(const char *source, int field) {
  return swap_endian<int32_t>(source + field - 1);
}

static inline std::string align_80(const std::string &s) {
  return fmt::format("{:<76}", s);
}

static inline std::string catstr(const std::string &s1, const std::string &s2) {
  return s1 + s2;
}

static inline std::string create_textual_header(
    const MetaInfo &meta_info, const std::string &segyname = "",
    const std::vector<std::string> &custom_info = std::vector<std::string>()) {
  std::vector<std::string> header;

  auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::string time_string(25, ' ');
  std::strftime(&time_string[0], time_string.size(), "%Y-%m-%d %H:%M:%S",
                std::localtime(&now));

  float Y_interval, Z_interval;
  if (meta_info.scalar != 0) {
    Y_interval = meta_info.scalar > 0
                     ? meta_info.Y_interval * meta_info.scalar
                     : meta_info.Y_interval / -meta_info.scalar;
    Z_interval = meta_info.scalar > 0
                     ? meta_info.Z_interval * meta_info.scalar
                     : meta_info.Z_interval / -meta_info.scalar;
  }

  std::string dformat = meta_info.data_format == 1
                            ? "4-bytes IBM floating-point"
                            : "4-bytes IEEE floating-point";

  if (custom_info.size() > 12) {
    fmt::print("The max rows of custom_info header are 12 rows, but the vector "
               "contian {} strings. We will drop the last strings.",
               custom_info.size());
  }

  if (custom_info.size() > 0) {
    for (int i = 0; i < 12; i++) {
      if (i >= custom_info.size()) {
        header.push_back(align_80(""));
      } else {
        header.push_back(align_80(custom_info[i]));
      }
    }
  } else {
    header.push_back(
        align_80("Create By CIGSEGY software (CIG, USTC, 2023), see:"));
    header.push_back(
        align_80("github: https://github.com/JintaoLee-Roger/cigsegy"));
    header.push_back(align_80(""));
    header.push_back(align_80(fmt::format("Name: {}", segyname)));
    header.push_back(align_80(
        fmt::format("Type: 3D seismic  Created Time: {}", time_string)));
    header.push_back(align_80(""));
    header.push_back(align_80(""));
    header.push_back(align_80(""));
    header.push_back(align_80(""));
    header.push_back(align_80(""));
    header.push_back(align_80(""));
    header.push_back(align_80(""));
  }

  header.push_back(align_80(fmt::format(
      "Inline range:    {} - {}", meta_info.min_inline, meta_info.max_inline)));
  header.push_back(
      align_80(fmt::format("Crossline range: {} - {}", meta_info.min_crossline,
                           meta_info.max_crossline)));
  header.push_back(
      align_80(fmt::format("Inline step: {}, Crossline step: {}",
                           meta_info.inline_step, meta_info.crossline_step)));
  header.push_back(align_80(
      fmt::format("shape: (n-time, n-crossline, n-time) = ({}, {}, {})",
                  meta_info.sizeX, meta_info.sizeY, meta_info.sizeZ)));
  header.push_back(align_80(
      fmt::format("Number of samples per data trace: {}", meta_info.sizeX)));
  header.push_back(
      align_80(fmt::format("Sample nterval: {} ", meta_info.sample_interval)));
  header.push_back(align_80(fmt::format("Data sample format: {}", dformat)));
  header.push_back(align_80(fmt::format(
      "interval of inline: {:.6g} meters, interval of crossline: {:.6g} meters",
      Z_interval, Y_interval)));
  header.push_back(
      align_80(fmt::format("Time start: {}", meta_info.start_time)));

  header.push_back(align_80(""));
  header.push_back(align_80(""));
  header.push_back(align_80(""));
  header.push_back(align_80("Binary header locations:"));
  header.push_back(
      align_80(fmt::format("Sample interval             : bytes {}",
                           field_str(kBSampleIntervalField))));
  header.push_back(
      align_80(fmt::format("Number of samples per trace : bytes {}",
                           field_str(kBSampleCountField))));
  header.push_back(
      align_80(fmt::format("Data sample format code     : bytes {}",
                           field_str(kBSampleFormatField))));
  header.push_back(align_80(""));
  header.push_back(align_80(""));
  header.push_back(align_80("Trace header locations:"));
  header.push_back(
      align_80(fmt::format("Inline number               : bytes {}",
                           field_str(meta_info.inline_field, 4))));
  header.push_back(
      align_80(fmt::format("Crossline number            : bytes {}",
                           field_str(meta_info.crossline_field, 4))));
  header.push_back(
      align_80(fmt::format("X coordinate                : bytes {}",
                           field_str(meta_info.X_field, 4))));
  header.push_back(
      align_80(fmt::format("Y coordinate                : bytes {}",
                           field_str(meta_info.Y_field, 4))));
  header.push_back(align_80(fmt::format(
      "Trace start time/depth      : bytes {}", field_str(kTStartTimeField))));
  header.push_back(
      align_80(fmt::format("Number of samples per trace : bytes {}",
                           field_str(kTSampleCountField))));
  header.push_back(
      align_80(fmt::format("Sample interval             : bytes {}",
                           field_str(kTSampleIntervalField))));
  header.push_back(align_80(""));
  header.push_back(align_80(""));

  assert(header.size() == kTextualRows);

  for (int i = 0; i < header.size(); i++) {
    header[i] = fmt::format("C{:02} ", i + 1) + header[i];
  }

  std::string out =
      std::accumulate(header.begin(), header.end(), std::string{}, catstr);

  assert(out.length() == kTextualHeaderSize);

  for (auto &s : out) {
    s = getEBCIDfromASCII(s);
  }

  return out;
}

static inline void modify_traceheader(char *trace_header, int sizeX,
                                      int start_time = -1) {
  int16_t *int_header = reinterpret_cast<int16_t *>(trace_header);
  int_header[(kTSampleCountField - 1) / 2] =
      swap_endian(static_cast<int16_t>(sizeX));
  if (start_time > 0) {
    int_header[(kTStartTimeField - 1) / 2] =
        swap_endian(static_cast<int16_t>(start_time));
    int_header[(kTDelayTimeField - 1) / 2] =
        swap_endian(static_cast<int16_t>(start_time));
  }
}

static int64_t copy_traces(const mio::mmap_source &header,
                           const std::vector<LineInfo> &line_info,
                           const MetaInfo &meta_info, const float *src,
                           char *outptr, int sizeX, int sizeY, int sizeZ,
                           int offsetX, int offsetY, int offsetZ,
                           int start_time) {
  // trace header and data
  char *start_outptr = outptr;

  int trace_count = 0;

  int step = sizeZ >= 100 ? 1 : 100 / sizeZ + 1;
  int nbar = sizeZ >= 100 ? sizeZ : step * sizeZ;
  progressbar bar(nbar);

  uint64_t trace_size = kTraceHeaderSize + sizeX * sizeof(float);
  uint64_t trace_size_ori = kTraceHeaderSize + meta_info.sizeX * sizeof(float);

  // iz: numpy sub volume
  // izo: numpy full volume
  for (int iz = 0; iz < sizeZ; iz++) {
    if (showpbar) {
      updatebar(bar, step);
    }

    const float *srcopy = src + iz * sizeY * sizeX;

    int izo = iz + offsetZ;
    const char *header_ptr = header.data() + kTextualHeaderSize +
                             kBinaryHeaderSize +
                             trace_size_ori * line_info[izo].trace_start;

    int count_in_line = 0;
    int xyzspace_start = 0;
    int xyzspace_end = 0;
    // iyo: segy index
    // iy: numpy sub volume
    // srct: numpy full volume
    for (int iyo = 0; iyo < line_info[izo].count; iyo++) {
      int srct = iyo;
      if (line_info[izo].count != meta_info.sizeY) {
        const int32_t *tmp = reinterpret_cast<const int32_t *>(header_ptr);
        srct = swap_endian(
                   tmp[(iyo * trace_size_ori + meta_info.crossline_field - 1) /
                       4]) -
               meta_info.min_crossline;

        srct /= meta_info.crossline_step;
      }
      if (iyo == 0) {
        xyzspace_start = srct;
      }
      if (iyo == line_info[izo].count) {
        xyzspace_end = srct + 1;
      }

      if (srct >= offsetY && srct < (offsetY + sizeY)) {
        count_in_line++;
        trace_count++;
        int iy = srct - offsetY;

        // copy header & modify sample count and start time
        memcpy(outptr, header_ptr + iyo * trace_size_ori, kTraceHeaderSize);
        if (sizeX != meta_info.sizeX) {
          modify_traceheader(outptr, sizeX, start_time);
        }
        outptr += kTraceHeaderSize;

        // copy data
        memcpy(outptr, srcopy + iy * sizeX, sizeX * sizeof(float));
        float *floatptr = reinterpret_cast<float *>(outptr);
        for (int ix = 0; ix < sizeX; ix++) {
          if (meta_info.data_format == 1) {
            floatptr[ix] = ieee_to_ibm(floatptr[ix], true);
          }
          floatptr[ix] = swap_endian(floatptr[ix]);
        }
        floatptr = nullptr;
        outptr += (sizeX * sizeof(float));
      }
    }

    if (count_in_line == 0) {
      throw std::runtime_error(fmt::format(
          "Can not copy headers through the index, becuase the line {} "
          "contains 0 trace if use the offsetY and sizeY."
          "In line {}, the header file's crossline range: "
          "[{}, {}], in the 3D volume space, y axis is range:"
          "[{}, {}), while your binary is range: [offsetY, offsetY+sizeY] "
          "= [{}, {}).",
          izo, izo, line_info[izo].crossline_start,
          line_info[izo].crossline_end, xyzspace_start, xyzspace_end, offsetY,
          offsetY + sizeY));
    }
  }
  assert((outptr - start_outptr) == (trace_count * trace_size));
  if (showpbar) {
    fmt::print("\n");
  }

  return trace_count * trace_size;
}

/********************* SegyIO public ***************************/

SegyIO::SegyIO(const std::string &segyname) {
  this->isReadSegy = true;
  memset(&this->m_metaInfo, 0, sizeof(MetaInfo));
  std::error_code error;
  this->m_source.map(segyname, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
  scanBinaryHeader();
}

SegyIO::SegyIO(int sizeX, int sizeY, int sizeZ) {
  this->isReadSegy = false;
  this->m_metaInfo.sizeX = sizeX;
  this->m_metaInfo.sizeY = sizeY;
  this->m_metaInfo.sizeZ = sizeZ;
  this->m_metaInfo.trace_count = sizeY * sizeZ;
  this->initMetaInfo();
}

SegyIO::SegyIO(const std::string &binaryname, int sizeX, int sizeY, int sizeZ) {
  this->isReadSegy = false;
  std::error_code error;
  this->m_source.map(binaryname, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
  this->m_metaInfo.sizeX = sizeX;
  this->m_metaInfo.sizeY = sizeY;
  this->m_metaInfo.sizeZ = sizeZ;
  this->m_metaInfo.trace_count = sizeY * sizeZ;
  this->initMetaInfo();
}

SegyIO::~SegyIO() {
  if (m_source.is_mapped()) {
    m_source.unmap();
  }
}

int64_t SegyIO::trace_count() { return m_metaInfo.trace_count; }

void SegyIO::set_size(int x, int y, int z) {
  m_metaInfo.sizeX = x;
  m_metaInfo.sizeY = y;
  m_metaInfo.sizeZ = z;
  if (isReadSegy) {
    m_metaInfo.isNormalSegy = true;
    isScan = true;
    int64_t trace_count =
        (m_source.size() - kTextualHeaderSize - kBinaryHeaderSize) /
        (kTraceHeaderSize + x * sizeof(float));
    if ((int64_t)y * z != (trace_count)) {
      throw std::runtime_error("invalid shape. inline * crossline != "
                               "total_trace_count");
    }
  }
}

MetaInfo SegyIO::get_metaInfo() { return m_metaInfo; }

void SegyIO::setInlineLocation(int loc) {
  if (loc <= 0) {
    throw std::runtime_error("Invalid location (must > 0)");
  }
  m_metaInfo.inline_field = loc;
  isScan = false;
}

void SegyIO::setCrosslineLocation(int loc) {
  if (loc <= 0) {
    throw std::runtime_error("Invalid location (must > 0)");
  }
  m_metaInfo.crossline_field = loc;
  isScan = false;
}

void SegyIO::setXLocation(int loc) {
  if (loc != 73 && loc != 181) {
    fmt::print("[Warning]: You set a unusual X field: {}, the best choice "
               "is set it as 73 or 181.\n",
               loc);
  }
  if (loc <= 0) {
    throw std::runtime_error("Invalid location (must > 0)");
  }
  m_metaInfo.X_field = loc;
  isScan = false;
}

void SegyIO::setYLocation(int loc) {
  if (loc != 77 && loc != 185) {
    fmt::print("[Warning]: You set a unusual Y field: {}, the best choice "
               "is set it as 77 or 185.\n",
               loc);
  }
  if (loc <= 0) {
    throw std::runtime_error("Invalid location (must > 0)");
  }
  m_metaInfo.Y_field = loc;
  isScan = false;
}

void SegyIO::setInlineStep(int step) {
  if (step <= 0) {
    throw std::runtime_error("Invalid inline step (must > 0)");
  }
  m_metaInfo.inline_step = step;
}

void SegyIO::setCrosslineStep(int step) {
  if (step <= 0) {
    throw std::runtime_error("Invalid crossline step (must > 0)");
  }
  m_metaInfo.crossline_step = step;
}

void SegyIO::setSteps(int istep, int xstep) {
  setInlineStep(istep);
  setCrosslineStep(xstep);
}

void SegyIO::setFillNoValue(float noValue) {
  m_metaInfo.fillNoValue = noValue;
  isScan = false;
}

void SegyIO::setSampleInterval(int dt) {
  if (dt <= 0) {
    throw std::runtime_error("Invalid Interval (must > 0)");
  }
  m_metaInfo.sample_interval = dt;
  isScan = false;
}

void SegyIO::setDataFormatCode(int dformat) {
  if (dformat != 1 && dformat != 5) {
    throw std::runtime_error("Don't support this data format now.");
  }
  m_metaInfo.data_format = dformat;
  isScan = false;
}

void SegyIO::setStartTime(int start_time) {
  if (start_time < 0) {
    throw std::runtime_error("Invalid start time (must >= 0)");
  }
  m_metaInfo.start_time = start_time;
  isScan = false;
}

void SegyIO::setInlineInterval(float dz) {
  if (dz <= 0) {
    throw std::runtime_error("Invalid interval (must > 0)");
  }
  m_metaInfo.Z_interval = dz * static_cast<float>(-m_metaInfo.scalar);
  isScan = false;
}

void SegyIO::setCrosslineInterval(float dy) {
  if (dy <= 0) {
    throw std::runtime_error("Invalid interval (must > 0)");
  }
  m_metaInfo.Y_interval = dy * static_cast<float>(-m_metaInfo.scalar);
  isScan = false;
}

void SegyIO::setMinInline(int in) {
  if (in <= 0) {
    throw std::runtime_error("Invalid line number (must > 0)");
  }
  m_metaInfo.min_inline = in;
  m_metaInfo.max_inline = in + m_metaInfo.sizeZ - 1;
  isScan = false;
}

void SegyIO::setMinCrossline(int cross) {
  if (cross <= 0) {
    throw std::runtime_error("Invalid crossline number (must > 0)");
  }
  m_metaInfo.min_crossline = cross;
  m_metaInfo.max_crossline = cross + m_metaInfo.sizeY - 1;
  isScan = false;
}

std::string SegyIO::textual_header(char coding) {
  if (!isReadSegy && m_sink.size() < kTextualHeaderSize) {
    throw std::runtime_error(
        "No textual header, because this is not a segy "
        "file (read mode) or you don't create textual header (create mode)");
  }
  const char *textual_header = nullptr;
  if (isReadSegy) {
    textual_header = m_source.data();
  } else {
    textual_header = m_sink.data();
  }
  char out[kTextualHeaderSize + kTextualRows];
  bool isEBCDIC = true;
  if (coding == 'a') {
    isEBCDIC = false;
  } else if (coding == 'u') {
    bool isEBCDIC = isTextInEBCDICFormat(textual_header, kTextualHeaderSize);
  }
  for (int iRow = 0; iRow < kTextualRows; iRow++) {
    int offset = iRow * kTextualColumns;
    for (int iCol = 0; iCol < kTextualColumns; iCol++) {
      if (isEBCDIC) {
        out[iCol + offset + iRow] =
            getASCIIfromEBCDIC(textual_header[iCol + offset]);
      } else {
        out[iCol + offset + iRow] = textual_header[iCol + offset];
      }
    }
    if (iRow < kTextualRows - 1) {
      out[(iRow + 1) * (kTextualColumns + 1) - 1] = '\n';
    }
  }

  return std::string(out);
}

std::string SegyIO::metaInfo() {
  if (!isScan && isReadSegy) {
    scan();
  }

  float Y_interval = 0;
  float Z_interval = 0;

  if (m_metaInfo.scalar != 0) {

    Y_interval = m_metaInfo.scalar > 0
                     ? m_metaInfo.Y_interval * m_metaInfo.scalar
                     : m_metaInfo.Y_interval / -m_metaInfo.scalar;
    Z_interval = m_metaInfo.scalar > 0
                     ? m_metaInfo.Z_interval * m_metaInfo.scalar
                     : m_metaInfo.Z_interval / -m_metaInfo.scalar;
  }

  if (!is_crossline_fast) {
    float t = Y_interval;
    Y_interval = Z_interval;
    Z_interval = t;
  }

  int sizeY, sizeZ, ifield, xfield, istart, iend, xstart, xend, istep, xstep;
  sizeY = is_crossline_fast ? m_metaInfo.sizeY : m_metaInfo.sizeZ;
  sizeZ = is_crossline_fast ? m_metaInfo.sizeZ : m_metaInfo.sizeY;
  ifield =
      is_crossline_fast ? m_metaInfo.inline_field : m_metaInfo.crossline_field;
  xfield =
      is_crossline_fast ? m_metaInfo.crossline_field : m_metaInfo.inline_field;
  istart = is_crossline_fast ? m_metaInfo.min_inline : m_metaInfo.min_crossline;
  iend = is_crossline_fast ? m_metaInfo.max_inline : m_metaInfo.max_crossline;
  xstart = is_crossline_fast ? m_metaInfo.min_crossline : m_metaInfo.min_inline;
  xend = is_crossline_fast ? m_metaInfo.max_crossline : m_metaInfo.max_inline;
  istep =
      is_crossline_fast ? m_metaInfo.inline_step : m_metaInfo.crossline_step;
  xstep =
      is_crossline_fast ? m_metaInfo.crossline_step : m_metaInfo.inline_step;

  std::string shapeinfo =
      is_crossline_fast
          ? fmt::format(
                "In python, the shape is (n-inline, n-crossline, n-time) "
                "= ({}, {}, {}).\n\n",
                sizeZ, sizeY, m_metaInfo.sizeX)
          : fmt::format(
                "In python, the shape is (n-crossline, n-inline, n-time) "
                "= ({}, {}, {}), as the fast order is inline\n"
                "You need transpose it manully, "
                "such as, in numpy, `d = d.transpose(1, 0, 2)`\n\n",
                sizeY, sizeZ, m_metaInfo.sizeX);

  std::string trace_sorting_string =
      fmt::format("{} - {}", m_metaInfo.trace_sorting_code,
                  kTraceSortingHelp.at(m_metaInfo.trace_sorting_code));

  std::string dformat = m_metaInfo.data_format == 1
                            ? "4-bytes IBM floating-point"
                            : "4-bytes IEEE floating-point";
  return fmt::format("{}shape: (n-time, n-crossline, n-inline) = ({}, {}, {})\n"
                     "trace sorting code: {}\n"
                     "sample interval: {}, data format code: {}\n"
                     "inline range: {} - {}, crossline range: {} - {}\n"
                     "interval of inline: {:.1f}, interval of crossline: "
                     "{:.1f}, time start: {}\n"
                     "inline field: {}, crossline field: {}\n"
                     "inline step: {}, crossline step: {}\n"
                     "Is regular file (no missing traces): {}",
                     shapeinfo, m_metaInfo.sizeX, sizeY, sizeZ,
                     trace_sorting_string, m_metaInfo.sample_interval, dformat,
                     istart, iend, xstart, xend, Z_interval, Y_interval,
                     m_metaInfo.start_time, ifield, xfield, istep, xstep,
                     m_metaInfo.isNormalSegy);
}

void SegyIO::get_binary_header_full(uchar *binheader, bool raw) {
  if (!isReadSegy) {
    throw std::runtime_error("get_binary_full func is used when read mode");
  }
  const char *src = m_source.data() + kTextualHeaderSize;
  if (raw) {
    memcpy(binheader, src, kBinaryHeaderSize);
  } else {
    read_binary_header(binheader, src);
  }
}

void SegyIO::get_trace_header_full(int n, uchar *traceheader, bool raw) {
  if (!isReadSegy) {
    throw std::runtime_error("get_binary_full func is used when read mode");
  }
  const char *src = m_source.data() + kTextualHeaderSize + kBinaryHeaderSize +
                    (m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize) *
                        static_cast<uint64_t>(n);
  if (raw) {
    memcpy(traceheader, src, kTraceHeaderSize);
  } else {
    read_one_trace_header(traceheader, src);
  }
}

void SegyIO::get_trace_full(int n, uchar *trace, bool raw) {
  if (!isReadSegy) {
    throw std::runtime_error("get_binary_full func is used when read mode");
  }
  const char *src = m_source.data() + kTextualHeaderSize + kBinaryHeaderSize +
                    (m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize) *
                        static_cast<uint64_t>(n);
  int tracesize = m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize;

  if (raw) {
    memcpy(trace, src, tracesize);
  } else {
    read_one_trace_header(trace, src);
    float *data = reinterpret_cast<float *>(trace + kTraceHeaderSize);
    memcpy(data, src + kTraceHeaderSize, m_metaInfo.sizeX * sizeof(float));
    for (int i = 0; i < m_metaInfo.sizeX; i++) {
      if (m_metaInfo.data_format == 1) {
        data[i] = ibm_to_ieee(data[i], true);
      } else {
        data[i] = swap_endian(data[i]);
      }
    }
  }
}

void SegyIO::collect(float *data, int beg, int end) {
  if (beg < 0) { // collect all traces
    beg = 0;
    end = m_metaInfo.trace_count;
  }
  if (end == 0) { // collect beg-th trace
    end = beg + 1;
  }
  if (end < 0) { // collect beg to the last traces
    end = m_metaInfo.trace_count;
  }
  if (beg >= end) {
    throw std::runtime_error("invalid range: beg >= end");
  }
  if (end > m_metaInfo.trace_count) {
    throw std::runtime_error("invalid range: end > trace_count");
  }

  int total = end - beg;

  uint64_t trace_size = m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize;
  const char *source = m_source.data() + kTextualHeaderSize + kBinaryHeaderSize;
  progressbar bar(100);
  for (int i = beg; i < end; i++) {
    if (i % (total / 100) == 0 && showpbar) {
      bar.update();
    }
    memcpy(data, source + i * trace_size + kTraceHeaderSize,
           m_metaInfo.sizeX * sizeof(float));
    for (int j = 0; j < m_metaInfo.sizeX; j++) {
      if (m_metaInfo.data_format == 1) {
        data[j] = ibm_to_ieee(data[j], true);
      } else {
        data[j] = swap_endian(data[j]);
      }
    }
    data += m_metaInfo.sizeX;
  }
}

void SegyIO::scan() {
  if (!isReadSegy) {
    throw std::runtime_error(
        "'scan()' function only used in reading segy mode.");
  }

  isScan = true;
  if (m_metaInfo.inline_field == 0) {
    m_metaInfo.inline_field = kDefaultInlineField;
  }

  if (m_metaInfo.crossline_field == 0) {
    m_metaInfo.crossline_field = kDefaultCrosslineField;
  }

  if (m_metaInfo.X_field == 0) {
    m_metaInfo.X_field = kDefaultXField;
  }
  if (m_metaInfo.Y_field == 0) {
    m_metaInfo.Y_field = kDefaultYField;
  }

  if (m_metaInfo.inline_step <= 0) {
    m_metaInfo.inline_step = 1;
  }
  if (m_metaInfo.crossline_step <= 0) {
    m_metaInfo.crossline_step = 1;
  }

  // get sizeZ, i.e. line_count
  uint64_t trace_size = m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize;
  const char *start = m_source.data() + kTextualHeaderSize + kBinaryHeaderSize;
  m_metaInfo.start_time = swap_endian<int16_t>(start + kTStartTimeField - 1);
  m_metaInfo.scalar = swap_endian<int16_t>(start + kTScalarField - 1);

  // check order, is the fast order inline or crossline?
  // the default is crossline
  // if the fast order is inline, then exchange the location
  // of inline and crossline
  check_order();

  // line x: ... trace1
  // line x+1: trace2 ...
  TraceInfo trace1{}, trace2{};
  _get_TraceInfo(0, trace1);
  _get_TraceInfo(m_metaInfo.trace_count - 1, trace2);

  m_metaInfo.sizeZ =
      (trace2.inline_num - trace1.inline_num) / m_metaInfo.inline_step + 1;
  m_metaInfo.min_inline = trace1.inline_num;
  m_metaInfo.max_inline = trace2.inline_num;

  // init this two field
  m_metaInfo.min_crossline = trace1.crossline_num;
  m_metaInfo.max_crossline = trace2.crossline_num;

  if (m_metaInfo.sizeZ > kMaxSizeOneDimemsion ||
      m_metaInfo.trace_count / m_metaInfo.sizeZ == 0 || m_metaInfo.sizeZ < 2) {
    throw std::runtime_error(
        "Size Z (inline number) is invalid, don't support. Maybe the "
        "inline location is wrong, use 'setInlineLocation(loc)' to set.");
  }

  m_metaInfo.isNormalSegy =
      m_metaInfo.trace_count % m_metaInfo.sizeZ == 0 ? true : false;
  m_lineInfo.resize(m_metaInfo.sizeZ);
  m_lineInfo[m_metaInfo.sizeZ - 1].crossline_end = trace2.crossline_num;

  // fill m_lineInfo
  int jump = m_metaInfo.trace_count / m_metaInfo.sizeZ;
  m_metaInfo.sizeY = jump;
  int itrace = 0;
  _get_TraceInfo(0, trace2);
  for (int i = 0; i < m_metaInfo.sizeZ - 1; i++) {
    m_lineInfo[i].trace_start = itrace;
    m_lineInfo[i].line_num = trace2.inline_num;
    m_lineInfo[i].count = 1;
    m_lineInfo[i].crossline_start = trace2.crossline_num;

    if (trace2.crossline_num < m_metaInfo.min_crossline) {
      m_metaInfo.min_crossline = trace2.crossline_num;
    }

    itrace += jump;
    if (itrace >= m_metaInfo.trace_count) {
      jump -= (itrace - m_metaInfo.trace_count + 1);
      itrace = m_metaInfo.trace_count - 1;
    }

    _get_TraceInfo(itrace, trace2);
    _get_TraceInfo(itrace - 1, trace1);

    if (trace2.inline_num == m_lineInfo[i].line_num) {
      m_metaInfo.isNormalSegy = false;
      while (trace2.inline_num !=
                 m_lineInfo[i].line_num + m_metaInfo.inline_step &&
             itrace < m_metaInfo.trace_count &&
             m_metaInfo.sizeY <= kMaxSizeOneDimemsion) {
        itrace++;
        jump++;
        if (jump > kMaxSizeOneDimemsion || itrace >= m_metaInfo.trace_count) {
          throw std::runtime_error(
              "inline/crossline location is wrong, use "
              "'setInlineLocation(loc)'/'setCrosslineLocation(loc)' to set");
        }
        _get_TraceInfo(itrace, trace2);
      }
      _get_TraceInfo(itrace - 1, trace1);
      // if (jump > m_metaInfo.sizeY) {
      //   m_metaInfo.sizeY = jump;
      // }
    } else if (trace1.inline_num > m_lineInfo[i].line_num) {
      m_metaInfo.isNormalSegy = false;
      while (trace1.inline_num != m_lineInfo[i].line_num && itrace > 0 &&
             jump > 0) {
        itrace--;
        jump--;
        if (jump <= 0 || itrace <= 0) {
          throw std::runtime_error(
              "inline/crossline location is wrong, use "
              "'setInlineLocation(loc)'/'setCrosslineLocation(loc)' to set");
        }
        trace2 = trace1;
        _get_TraceInfo(itrace - 1, trace1);
      }
      _get_TraceInfo(itrace, trace2);
    }

    m_metaInfo.sizeY = (m_metaInfo.max_crossline - m_metaInfo.min_crossline) /
                           m_metaInfo.crossline_step +
                       1;

    if (trace2.inline_num == m_lineInfo[i].line_num + m_metaInfo.inline_step &&
        trace1.inline_num == m_lineInfo[i].line_num) {
      if (trace1.crossline_num > m_metaInfo.max_crossline) {
        m_metaInfo.max_crossline = trace1.crossline_num;
      }
      m_lineInfo[i].trace_end = itrace - 1;
      m_lineInfo[i].count = itrace - m_lineInfo[i].trace_start;
      m_lineInfo[i].crossline_end = trace1.crossline_num;
    } else {
      throw std::runtime_error("Cannot analysis this segy file, "
                               "may inline step != 1");
    }
  }

  if (m_metaInfo.isNormalSegy &&
      m_metaInfo.trace_count != (m_metaInfo.sizeY * m_metaInfo.sizeZ)) {
    throw std::runtime_error(
        "isNormal is true while sizeY * sizeZ != trace_count, "
        "maybe inline/crossline location is wrong");
  }

  // the last line
  m_lineInfo[m_metaInfo.sizeZ - 1].trace_start = itrace;
  m_lineInfo[m_metaInfo.sizeZ - 1].line_num = trace2.inline_num;
  m_lineInfo[m_metaInfo.sizeZ - 1].trace_end = m_metaInfo.trace_count - 1;
  m_lineInfo[m_metaInfo.sizeZ - 1].count = m_metaInfo.trace_count - itrace;
  m_lineInfo[m_metaInfo.sizeZ - 1].crossline_start = trace2.crossline_num;

  // cal x, y interval
  _get_TraceInfo(0, trace1);
  _get_TraceInfo(m_lineInfo[0].trace_end, trace2);
  // float scale =
  //     m_metaInfo.scalar < 0 ? -1.0 / m_metaInfo.scalar : m_metaInfo.scalar;
  m_metaInfo.Y_interval = std::sqrt(std::pow(trace2.X - trace1.X, 2) +
                                    std::pow(trace2.Y - trace1.Y, 2)) /
                          (m_lineInfo[0].count - 1);
  int num = m_metaInfo.trace_count > 10 ? 10 : m_metaInfo.trace_count - 1;
  _get_TraceInfo(m_lineInfo[num].trace_start, trace2);
  m_metaInfo.Z_interval =
      std::sqrt(std::pow(trace2.X - trace1.X, 2) +
                std::pow(trace2.Y - trace1.Y, 2) -
                std::pow(float(trace2.crossline_num - trace1.crossline_num) /
                             m_metaInfo.crossline_step * m_metaInfo.Y_interval,
                         2)) /
      ((trace2.inline_num - trace1.inline_num) / float(m_metaInfo.inline_step));
}

void SegyIO::read(float *dst, int startX, int endX, int startY, int endY,
                  int startZ, int endZ) {
  if (!isReadSegy) {
    throw std::runtime_error(
        "'read()' function used only in reading segy mode");
  }
  if (startX >= endX || startY >= endY || startZ >= endZ) {
    throw std::runtime_error("Index 'end' must large than 'start'");
  }
  if (startX < 0 || endX > m_metaInfo.sizeX || startY < 0 ||
      endY > m_metaInfo.sizeY || startZ < 0 || endZ > m_metaInfo.sizeZ) {
    throw std::runtime_error("Index out of range");
  }

  const char *source = m_source.data() + kTextualHeaderSize + kBinaryHeaderSize;
  uint64_t trace_size = m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize;

  int sizeX = endX - startX;
  int sizeY = endY - startY;
  int sizeZ = endZ - startZ;
  int offset = startX * sizeof(float) + kTraceHeaderSize;

  int step = sizeZ >= 100 ? 1 : 100 / sizeZ + 1;
  int nbar = sizeZ >= 100 ? sizeZ : step * sizeZ;
  progressbar bar(nbar);

  auto time_start = std::chrono::high_resolution_clock::now();

  // #pragma omp parallel for
  for (int iZ = startZ; iZ < endZ; iZ++) {
    // init start trace in iZ line
    int istart = startY;
    float *dstline = dst + static_cast<uint64_t>(iZ - startZ) * sizeX * sizeY;
    uint64_t trace_start = m_metaInfo.isNormalSegy ? iZ * m_metaInfo.sizeY
                                                   : m_lineInfo[iZ].trace_start;
    const char *sourceline = source + trace_start * trace_size;

    bool normal = true;
    // find istart
    if (!m_metaInfo.isNormalSegy) {
      normal = m_lineInfo[iZ].count == m_metaInfo.sizeY ? true : false;
      if (!normal) {
        if (istart > m_lineInfo[iZ].count) {
          // fill all zeros in this inline
          istart = m_lineInfo[iZ].count + 100;
        } else {
          int dst_crossline =
              m_metaInfo.min_crossline + startY * m_metaInfo.crossline_step;
          while (getCrossline(sourceline + istart * trace_size,
                              m_metaInfo.crossline_field) > dst_crossline &&
                 istart > 0) {
            istart--;
          }
        }
      }
    }

    for (int iY = startY; iY < endY; iY++) {
      float *dsttrace = dstline + (iY - startY) * sizeX;
      if (normal ||
          (istart <= m_lineInfo[iZ].count &&
           (getCrossline(sourceline + istart * trace_size,
                         m_metaInfo.crossline_field) ==
            (m_metaInfo.min_crossline + iY * m_metaInfo.crossline_step)))) {
        memcpy(dsttrace, sourceline + istart * trace_size + offset,
               sizeX * sizeof(float));
        for (int iX = 0; iX < sizeX; iX++) {
          if (m_metaInfo.data_format == 1) {
            dsttrace[iX] = ibm_to_ieee(dsttrace[iX], true);
          } else if (m_metaInfo.data_format == 5) {
            dsttrace[iX] = swap_endian(dsttrace[iX]);
          } else {
            throw std::runtime_error("Unsuport sample format");
          }
        }
        istart++;
      } else {
        std::fill(dsttrace, dsttrace + sizeX, m_metaInfo.fillNoValue);
      }
    }
    // #pragma omp critical
    if (showpbar) {
      updatebar(bar, step);
    }
  }
  if (showpbar) {
    fmt::print("\n");
  }

  auto time_end = std::chrono::high_resolution_clock::now();

  fmt::print("need time: {}s\n",
             std::chrono::duration_cast<std::chrono::nanoseconds>(time_end -
                                                                  time_start)
                     .count() *
                 1e-9);
}

void SegyIO::read(float *dst) {
  if (!isScan) {
    scan();
  }
  read(dst, 0, m_metaInfo.sizeX, 0, m_metaInfo.sizeY, 0, m_metaInfo.sizeZ);
}

void SegyIO::read_inline_slice(float *dst, int iZ) {
  if (!isScan) {
    scan();
  }
  read(dst, 0, m_metaInfo.sizeX, 0, m_metaInfo.sizeY, iZ, iZ + 1);
}

void SegyIO::read_cross_slice(float *dst, int iY) {
  if (!isScan) {
    scan();
  }
  read(dst, 0, m_metaInfo.sizeX, iY, iY + 1, 0, m_metaInfo.sizeZ);
}

void SegyIO::read_time_slice(float *dst, int iX) {
  if (!isScan) {
    scan();
  }
  read(dst, iX, iX + 1, 0, m_metaInfo.sizeY, 0, m_metaInfo.sizeZ);
}

void SegyIO::read_trace(float *dst, int iY, int iZ) {
  if (!isScan) {
    scan();
  }
  read(dst, 0, m_metaInfo.sizeX, iY, iY + 1, iZ, iZ + 1);
}

void SegyIO::tofile(const std::string &binary_out_name) {
  if (!isScan) {
    scan();
  }
  uint64_t need_size = static_cast<uint64_t>(m_metaInfo.sizeX) *
                       m_metaInfo.sizeY * m_metaInfo.sizeZ * sizeof(float);
  std::ofstream ffst(binary_out_name, std::ios::binary);
  if (!ffst) {
    throw std::runtime_error("create file failed");
  }
  for (int i = 0; i < int(need_size / kMaxLSeekSize) + 1; i++) {
    uint64_t move_point = need_size > kMaxLSeekSize ? kMaxLSeekSize : need_size;
    ffst.seekp(move_point - 1, std::ios_base::cur);
    ffst.put(0);
    need_size -= move_point;
  }
  if (need_size != 0) {
    throw std::runtime_error("create file failed");
  }
  ffst.close();

  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(binary_out_name, error);
  if (error) {
    throw std::runtime_error("mmap fail when write data");
  }
  // or need split into serveral chunks?
  read(reinterpret_cast<float *>(rw_mmap.data()));
  rw_mmap.unmap();
}

void SegyIO::cut(const std::string &outname, int startX, int endX, int startY,
                 int endY, int startZ, int endZ,
                 const std::vector<std::string> &custom_info) {
  if (!isScan) {
    scan();
  }

  if (startX >= endX || startY >= endY || startZ >= endZ) {
    throw std::runtime_error("Index 'end' must large than 'start'");
  }
  if (startX < 0 || endX > m_metaInfo.sizeX || startY < 0 ||
      endY > m_metaInfo.sizeY || startZ < 0 || endZ > m_metaInfo.sizeZ) {
    throw std::runtime_error("Index out of range");
  }

  int sizeX, sizeY, sizeZ;
  sizeX = endX - startX;
  sizeY = endY - startY;
  sizeZ = endZ - startZ;
  int start_time =
      m_metaInfo.start_time + startX * m_metaInfo.sample_interval / 1000;

  std::fstream out_(outname, std::ios::binary | std::ios::out);
  if (!out_.is_open()) {
    throw std::runtime_error("create file failed");
  }

  // textual header
  MetaInfo subinfo = m_metaInfo;
  subinfo.sizeX = sizeX;
  subinfo.sizeY = sizeY;
  subinfo.sizeZ = sizeZ;
  subinfo.start_time = start_time;
  subinfo.min_inline = m_metaInfo.min_inline + startZ * m_metaInfo.inline_step;
  subinfo.min_crossline =
      m_metaInfo.min_crossline + startY * m_metaInfo.crossline_step;
  subinfo.max_inline =
      subinfo.min_inline + (sizeZ - 1) * m_metaInfo.inline_step;
  subinfo.max_crossline =
      subinfo.min_crossline + (sizeY - 1) * m_metaInfo.crossline_step;

  std::string textual_header =
      create_textual_header(subinfo, outname, custom_info);
  out_.write(textual_header.data(), kTextualHeaderSize);

  // binary header
  std::vector<char> tmp(kBinaryHeaderSize, 0);
  memcpy(tmp.data(), m_source.data() + kTextualHeaderSize, kBinaryHeaderSize);
  if (sizeX != m_metaInfo.sizeX) {
    int16_t *int16ptr = reinterpret_cast<int16_t *>(tmp.data());
    int16ptr[(kBSampleCountField - 1) / 2] =
        swap_endian(static_cast<int16_t>(sizeX));
    int16ptr = nullptr;
  }
  out_.write(tmp.data(), kBinaryHeaderSize);

  int step = sizeZ >= 100 ? 1 : 100 / sizeZ + 1;
  int nbar = sizeZ >= 100 ? sizeZ : step * sizeZ;
  progressbar bar(nbar);

  uint64_t trace_size = kTraceHeaderSize + sizeX * sizeof(float);
  uint64_t trace_size_ori = kTraceHeaderSize + m_metaInfo.sizeX * sizeof(float);
  tmp.resize(trace_size);
  tmp.assign(trace_size, 0);
  int tracecount = 0;

  // iy: sub-voume dat space
  // iyo: full segy trace space
  // srct: full-volume dat space
  for (int iz = 0; iz < sizeZ; iz++) {
    if (showpbar) {
      updatebar(bar, step);
    }

    int izo = iz + startZ;
    const char *srcptr = m_source.data() + kTextualHeaderSize +
                         kBinaryHeaderSize +
                         trace_size_ori * m_lineInfo[izo].trace_start;
    int count_in_line = 0;
    int xyzspace_start = 0;
    int xyzspace_end = 0;
    for (int iyo = 0; iyo < m_lineInfo[izo].count; iyo++) {
      const char *traceptr = srcptr + iyo * trace_size_ori;
      int srct = iyo;
      if (m_lineInfo[izo].count != m_metaInfo.sizeY) {
        const int32_t *trace_tmp = reinterpret_cast<const int32_t *>(traceptr);
        srct = swap_endian(trace_tmp[(m_metaInfo.crossline_field - 1) / 4]) -
               m_metaInfo.min_crossline;
        srct /= m_metaInfo.crossline_step; // in dat space
      }
      if (iyo == 0) {
        xyzspace_start = srct;
      }
      if (iyo == m_lineInfo[izo].count - 1) {
        xyzspace_end = srct + 1;
      }

      if (srct >= startY && srct < endY) {
        tracecount++;
        count_in_line++;
        if (sizeX == m_metaInfo.sizeX) {
          memcpy(tmp.data(), traceptr, trace_size);
          out_.write(tmp.data(), trace_size);
        } else {
          memcpy(tmp.data(), traceptr, kTraceHeaderSize);
          memcpy(tmp.data() + kTraceHeaderSize,
                 traceptr + kTraceHeaderSize + startX * sizeof(float),
                 sizeX * sizeof(float));
          modify_traceheader(tmp.data(), sizeX, start_time);
          out_.write(tmp.data(), trace_size);
        }
      }
    } // end iyo

    if (count_in_line == 0) {
      throw std::runtime_error(
          fmt::format("Can not cut through the index, becuase the line {} "
                      "contains 0 trace if use the startY and endY. "
                      "In line {}, the original crossline range: "
                      "[{}, {}], in the 3D volume space, y axis is range:"
                      "[{}, {}), while your cut is range of: [{}, {}).",
                      izo, izo, m_lineInfo[izo].crossline_start,
                      m_lineInfo[izo].crossline_end, xyzspace_start,
                      xyzspace_end, startY, endY));
    }
  } // end iz
  auto pos = out_.tellp();
  assert(pos ==
         (kTextualHeaderSize + kBinaryHeaderSize + tracecount * trace_size));
  out_.close();

  if (showpbar) {
    fmt::print("\n");
  }
}

void SegyIO::cut(const std::string &outname, int startY, int endY, int startZ,
                 int endZ, const std::vector<std::string> &custom_info) {
  if (!isScan) {
    scan();
  }

  cut(outname, 0, m_metaInfo.sizeX, startY, endY, startZ, endZ, custom_info);
}

void SegyIO::cut(const std::string &outname, int startX, int endX,
                 const std::vector<std::string> &custom_info) {
  if (!isScan) {
    scan();
  }

  cut(outname, startX, endX, 0, m_metaInfo.sizeY, 0, m_metaInfo.sizeZ,
      custom_info);
}

void SegyIO::create(const std::string &segy_out_name, const float *src,
                    const std::vector<std::string> &custom_info) {
  if (isReadSegy) {
    throw std::runtime_error(
        "'create() function only can be used for creating segy file.'");
  }
  uint64_t need_size =
      kTextualHeaderSize + kBinaryHeaderSize +
      static_cast<uint64_t>(m_metaInfo.sizeY) * m_metaInfo.sizeZ *
          (m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize);
  std::ofstream ffst(segy_out_name, std::ios::binary);
  if (!ffst) {
    throw std::runtime_error("create file failed");
  }
  // int fd = open(binary_out_name.c_str(), O_RDWR | O_CREAT | O_TRUNC, 00644);
  for (int i = 0; i < int(need_size / kMaxLSeekSize) + 1; i++) {
    uint64_t move_point = need_size > kMaxLSeekSize ? kMaxLSeekSize : need_size;
    // if (lseek(fd, move_point - 1, SEEK_END) < 0) {
    //   throw std::runtime_error("create file failed");
    // }
    // if (write(fd, "", 1) < 0) {
    //   throw std::runtime_error("create file failed");
    // }
    ffst.seekp(move_point - 1, std::ios_base::cur);
    ffst.put(0);
    need_size -= move_point;
  }
  if (need_size != 0) {
    throw std::runtime_error("create file failed");
  }
  ffst.close();

  std::error_code error;
  mio::mmap_sink rw_mmap = mio::make_mmap_sink(segy_out_name, error);
  if (error) {
    throw std::runtime_error("mmap fail when write data");
  }

  write_textual_header(rw_mmap.data(), segy_out_name, custom_info);
  write_binary_header(rw_mmap.data() + kTextualHeaderSize);
  TraceHeader trace_header{};
  initTraceHeader(&trace_header);
  char *dst = rw_mmap.data() + kTextualHeaderSize + kBinaryHeaderSize;
  char *dstline = dst;

  int step = m_metaInfo.sizeZ >= 100 ? 1 : 100 / m_metaInfo.sizeZ + 1;
  int nbar =
      m_metaInfo.sizeZ >= 100 ? m_metaInfo.sizeZ : step * m_metaInfo.sizeZ;
  progressbar bar(nbar);

  uint64_t trace_size = m_metaInfo.sizeX * sizeof(float) + kTraceHeaderSize;
  // #pragma omp parallel for
  for (int iZ = 0; iZ < m_metaInfo.sizeZ; iZ++) {
    for (int iY = 0; iY < m_metaInfo.sizeY; iY++) {
      // write header
      int64_t x = iY * m_metaInfo.Y_interval + 5200;
      int64_t y = iZ * m_metaInfo.Z_interval + 5200;
      write_trace_header(dstline, &trace_header, iY + m_metaInfo.min_crossline,
                         iZ + m_metaInfo.min_inline, x, y);

      // copy data
      float *dstdata = reinterpret_cast<float *>(dstline + kTraceHeaderSize);
      const float *srcline =
          src + static_cast<uint64_t>(iY) * m_metaInfo.sizeX +
          static_cast<uint64_t>(iZ) * m_metaInfo.sizeX * m_metaInfo.sizeY;
      memcpy(dstdata, srcline, m_metaInfo.sizeX * sizeof(float));
      for (int iX = 0; iX < m_metaInfo.sizeX; iX++) {
        if (m_metaInfo.data_format == 1) {
          dstdata[iX] = ieee_to_ibm(dstdata[iX], true);
        }
        dstdata[iX] = swap_endian(dstdata[iX]);
      }

      dstline += trace_size;
    }
    // #pragma omp critical
    if (showpbar) {
      updatebar(bar, step);
    }
  }
  if (showpbar) {
    fmt::print("\n");
  }
  rw_mmap.unmap();
}

void SegyIO::create(const std::string &segy_out_name,
                    const std::vector<std::string> &custom_info) {
  if (isReadSegy) {
    throw std::runtime_error("Now is read segy mode, cannot create a segy");
  }
  if (!m_source.is_mapped()) {
    throw std::runtime_error("You need to read a binary file before create, or "
                             "you can create from memory");
  }
  create(segy_out_name, (float *)m_source.data(), custom_info);
}

/********************** SegyIO private *******************************/

void SegyIO::initMetaInfo() {
  m_metaInfo.isNormalSegy = true;
  m_metaInfo.crossline_field = kDefaultCrosslineField;
  m_metaInfo.inline_field = kDefaultInlineField;
  m_metaInfo.min_inline = 1;
  m_metaInfo.max_inline = m_metaInfo.min_inline + m_metaInfo.sizeZ - 1;
  m_metaInfo.min_crossline = 1;
  m_metaInfo.max_crossline = m_metaInfo.min_crossline + m_metaInfo.sizeY - 1;
  m_metaInfo.data_format = 5;
  m_metaInfo.sample_interval = 2000;
  m_metaInfo.Y_interval = 25 * 100;
  m_metaInfo.Z_interval = 25 * 100;
  m_metaInfo.scalar = -100;
  m_metaInfo.start_time = 0;
  m_metaInfo.X_field = kDefaultXField;
  m_metaInfo.Y_field = kDefaultYField;
}

void SegyIO::scanBinaryHeader() {
  const auto *binary_header = reinterpret_cast<const BinaryHeader *>(
      m_source.data() + kTextualHeaderSize);
  m_metaInfo.data_format = swap_endian(binary_header->data_format);
  m_metaInfo.sizeX = swap_endian(binary_header->trace_length);
  m_metaInfo.sample_interval = swap_endian(binary_header->sample_interval);
  m_metaInfo.trace_count =
      (m_source.size() - kTextualHeaderSize - kBinaryHeaderSize) /
      (kTraceHeaderSize + m_metaInfo.sizeX * sizeof(float));

  m_metaInfo.trace_sorting_code =
      swap_endian(binary_header->trace_sorting_code);
}

void SegyIO::check_order() {
  int i1, x1, i2, x2, i3, x3;
  TraceInfo trace{};
  _get_TraceInfo(0, trace);
  i1 = trace.inline_num;
  x1 = trace.crossline_num;
  _get_TraceInfo(1, trace);
  i2 = trace.inline_num;
  x2 = trace.crossline_num;
  _get_TraceInfo(2, trace);
  i3 = trace.inline_num;
  x3 = trace.crossline_num;

  if (x1 == x2 && x1 == x3) {
    if (i1 != i2 && i1 != i3 && i2 != i3) {
      is_crossline_fast = false;
      int t = m_metaInfo.inline_field;
      int t2 = m_metaInfo.inline_step;
      m_metaInfo.inline_field = m_metaInfo.crossline_field;
      m_metaInfo.crossline_field = t;
      m_metaInfo.inline_step = m_metaInfo.crossline_step;
      m_metaInfo.crossline_step = t2;

      fmt::print(
          "[Warining] The fast order of you segy file "
          "is inline order (default is crossline order). This means "
          "the file you obtain (numpy array or a binary file in disk) "
          "is shape as (n-time, n-inline, n-crossline) (in python numpy array"
          ", the shape is (n-crossline, n-inline, n-time)). "
          "You need to transpose it manully, such as "
          "`d = d.transpose(1, 0, 2)`\n\n");
    } else {
      throw std::runtime_error(fmt::format(
          "Cannot check the fast order (default is crossline). "
          "We check the first three traces, now the crossline numbers "
          "are constant, the inline numbers must be different if "
          "the file is valid (fast order is inline). "
          "But the inline number is {}, {}, {}. May be the locations "
          "is wrong.\n\n",
          i1, i2, i3));
    }
  } else if (i1 != i2 || i1 != i3 || i2 != i3) {
    throw std::runtime_error(fmt::format(
        "Cannot check the fast order (default is crossline). "
        "We check the first three traces, and can not evaluate "
        "the fast order, becuase the inline and crossline numbers are "
        "both different. The inline numbers: {}, {}, {}"
        "The crossline numbers: {}, {}, {}. Maybe the locations is "
        "wrong, or the file is small\n\n",
        i1, i2, i3, x1, x2, x3));
  }
}

void SegyIO::write_textual_header(char *dst, const std::string &segy_out_name,
                                  const std::vector<std::string> &custom_info) {
  std::string textual_header =
      create_textual_header(m_metaInfo, segy_out_name, custom_info);

  memcpy(dst, textual_header.c_str(), kTextualHeaderSize);
}

void SegyIO::write_binary_header(char *dst) {
  memset(dst, 0, kBinaryHeaderSize);
  BinaryHeader binary_header{};
  // auto &binary_header = (BinaryHeader &)dst;
  binary_header.jobID = swap_endian(int32_t(1));
  binary_header.line_number = swap_endian(m_metaInfo.min_inline);
  binary_header.num_traces_per_ensemble =
      swap_endian(int16_t(m_metaInfo.min_crossline));
  binary_header.sample_interval = swap_endian(m_metaInfo.sample_interval);
  binary_header.trace_length = swap_endian(int16_t(m_metaInfo.sizeX));
  binary_header.data_format = swap_endian(m_metaInfo.data_format);
  binary_header.ensemble_fold = swap_endian(int16_t(1));
  binary_header.trace_sorting_code = swap_endian(int16_t(4));
  binary_header.measurement_system = swap_endian(int16_t(1));
  binary_header.fixed_length_trace = swap_endian(int16_t(1));
  memcpy(dst, &binary_header, kBinaryHeaderSize);
}

void SegyIO::initTraceHeader(TraceHeader *trace_header) {
  memset(trace_header, 0, kTraceHeaderSize);
  trace_header->trace_sequence_number_in_line = swap_endian(int32_t(1));
  trace_header->trace_num_in_orig = swap_endian(int32_t(1));
  trace_header->trace_num_in_ensemble = swap_endian(int32_t(1));
  trace_header->trace_ID_code = swap_endian(int16_t(1));
  trace_header->data_used_for = swap_endian(int16_t(1));
  trace_header->scalar_for_elev_and_depth = swap_endian(int16_t(1));
  trace_header->scalar_for_coord = swap_endian(m_metaInfo.scalar);
  trace_header->coord_units = swap_endian(int16_t(1));
  trace_header->lag_time_A = swap_endian(m_metaInfo.start_time);
  trace_header->delay_record_time = swap_endian(m_metaInfo.start_time);
  trace_header->num_sample = swap_endian(int16_t(m_metaInfo.sizeX));
  trace_header->sample_interval = swap_endian(m_metaInfo.sample_interval);
  trace_header->correlated = swap_endian(int16_t(1));
  trace_header->sweep_type_code = swap_endian(int16_t(1));
  trace_header->taper_type = swap_endian(int16_t(1));
}

void SegyIO::write_trace_header(char *dst, TraceHeader *trace_header,
                                int32_t iY, int32_t iZ, int32_t x, int32_t y) {
  trace_header->trace_sequence_number_in_file = swap_endian(iZ);
  trace_header->orig_field_num = swap_endian(iZ);
  trace_header->source_point_num = swap_endian(iY);
  trace_header->ensemble_num = swap_endian(iY);
  trace_header->source_coord_X = swap_endian(x);
  trace_header->source_coord_Y = swap_endian(y);
  trace_header->X = swap_endian(x);
  trace_header->Y = swap_endian(y);
  trace_header->inline_num = swap_endian(iZ);
  trace_header->crossline_num = swap_endian(iY);
  memcpy(dst, trace_header, kTraceHeaderSize);
}

void SegyIO::close_file() {
  if (m_source.is_mapped()) {
    m_source.unmap();
  }
}

/************** Usefual function, need to be binded ******************/

void read(const std::string &segy_name, float *dst, int iline, int xline,
          int istep, int xstep) {
  SegyIO segy_data(segy_name);
  segy_data.setInlineLocation(iline);
  segy_data.setCrosslineLocation(xline);
  segy_data.setSteps(istep, xstep);
  segy_data.scan();
  segy_data.read(dst);
  segy_data.close_file();
}

void tofile(const std::string &segy_name, const std::string &out_name,
            int iline, int xline, int istep, int xstep) {
  SegyIO segy_data(segy_name);
  segy_data.setInlineLocation(iline);
  segy_data.setCrosslineLocation(xline);
  segy_data.setSteps(istep, xstep);
  segy_data.scan();
  segy_data.tofile(out_name);
  segy_data.close_file();
}

void read_ignore_header(const std::string &segy_name, float *dst, int sizeX,
                        int sizeY, int sizeZ, int format) {
  SegyIO segy_data(segy_name);
  segy_data.setDataFormatCode(format);
  segy_data.set_size(sizeX, sizeY, sizeZ);
  segy_data.read(dst, 0, sizeX, 0, sizeY, 0, sizeZ);
  segy_data.close_file();
}

void tofile_ignore_header(const std::string &segy_name,
                          const std::string &out_name, int sizeX, int sizeY,
                          int sizeZ, int format) {
  SegyIO segy_data(segy_name);
  segy_data.setDataFormatCode(format);
  segy_data.set_size(sizeX, sizeY, sizeZ);
  segy_data.tofile(out_name);
  segy_data.close_file();
}

void create_by_sharing_header(const std::string &segy_name,
                              const std::string &header_segy, const float *src,
                              int sizeX, int sizeY, int sizeZ, int iline,
                              int xline, int istep, int xstep, int offsetX,
                              int offsetY, int offsetZ,
                              const std::vector<std::string> &custom_info) {
  SegyIO header(header_segy);
  header.setInlineLocation(iline);
  header.setCrosslineLocation(xline);
  header.setSteps(istep, xstep);
  header.scan();
  auto line_info = header.line_info();
  auto meta_info = header.get_metaInfo();
  auto trace_count = header.trace_count();
  header.close_file();

  std::error_code error;
  mio::mmap_source m_source;
  m_source.map(header_segy, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }

  if (offsetX < 0 || offsetY < 0 || offsetZ < 0) {
    // full copy
    if (meta_info.sizeY != sizeY || meta_info.sizeZ != sizeZ ||
        meta_info.sizeX != sizeX) {
      throw std::runtime_error(fmt::format(
          "shape of header: {} x {} x {}, but shape of source: {} x {} x {}",
          meta_info.sizeZ, meta_info.sizeY, meta_info.sizeX, sizeZ, sizeY,
          sizeX));
    }

    uint64_t need_size =
        kTextualHeaderSize + kBinaryHeaderSize +
        trace_count * (sizeX * sizeof(float) + kTraceHeaderSize);
    std::ofstream ffst(segy_name, std::ios::binary);
    if (!ffst) {
      throw std::runtime_error("create file failed");
    }
    for (int i = 0; i < int(need_size / kMaxLSeekSize) + 1; i++) {
      uint64_t move_point =
          need_size > kMaxLSeekSize ? kMaxLSeekSize : need_size;
      ffst.seekp(move_point - 1, std::ios_base::cur);
      ffst.put(0);
      need_size -= move_point;
    }
    if (need_size != 0) {
      throw std::runtime_error("create file failed");
    }
    ffst.close();

    mio::mmap_sink rw_mmap = mio::make_mmap_sink(segy_name, error);
    if (error) {
      throw std::runtime_error("mmap fail when write data");
    }

    char *outptr = rw_mmap.data();
    memcpy(outptr, m_source.data(), kTextualHeaderSize + kBinaryHeaderSize);
    outptr += (kTextualHeaderSize + kBinaryHeaderSize);

    copy_traces(m_source, line_info, meta_info, src, outptr, sizeX, sizeY,
                sizeZ, 0, 0, 0, -1);

    rw_mmap.unmap();

  } else if (offsetX >= 0 && offsetY >= 0 && offsetZ >= 0) {
    // subset data
    if (sizeX + offsetX > meta_info.sizeX ||
        sizeY + offsetY > meta_info.sizeY ||
        sizeZ + offsetZ > meta_info.sizeZ) {
      throw std::runtime_error("size_{src} + offset > size_{header}");
    }

    int max_size = kTextualHeaderSize + kBinaryHeaderSize +
                   (kTraceHeaderSize + sizeX * sizeof(float)) * (sizeY * sizeZ);
    std::vector<char> outsegy(max_size, 0);
    char *outptr = outsegy.data();
    int start_time =
        meta_info.start_time + offsetX * meta_info.sample_interval / 1000;
    MetaInfo msub = meta_info;
    msub.sizeX = sizeX;
    msub.sizeY = sizeY;
    msub.sizeZ = sizeZ;
    msub.start_time = start_time;
    msub.min_inline = meta_info.min_inline + offsetZ * meta_info.inline_step;
    msub.min_crossline =
        meta_info.min_crossline + offsetY * meta_info.crossline_step;
    msub.max_inline = msub.min_inline + (sizeZ - 1) * meta_info.inline_step;
    msub.max_crossline =
        msub.min_crossline + (sizeY - 1) * meta_info.crossline_step;

    std::string textual_header =
        create_textual_header(msub, segy_name, custom_info);
    memcpy(outptr, textual_header.c_str(), kTextualHeaderSize);
    outptr += kTextualHeaderSize;

    // copy binary header & modify sample count
    memcpy(outptr, m_source.data() + kTextualHeaderSize, kBinaryHeaderSize);
    int16_t *int16ptr = reinterpret_cast<int16_t *>(outptr);
    int16ptr[(kBSampleCountField - 1) / 2] =
        swap_endian(static_cast<int16_t>(sizeX));
    int16ptr = nullptr;
    outptr += kBinaryHeaderSize;

    int64_t size =
        copy_traces(m_source, line_info, meta_info, src, outptr, sizeX, sizeY,
                    sizeZ, offsetX, offsetY, offsetZ, start_time);

    outsegy.resize(kTextualHeaderSize + kBinaryHeaderSize + size);

    std::ofstream file(segy_name, std::ios::binary | std::ios::out);
    if (file.is_open()) {
      file.write(outsegy.data(), outsegy.size());
      file.close();
    } else {
      throw std::runtime_error("create file failed");
    }
  } else {
    throw std::runtime_error(
        fmt::format("offset (X, Y, Z) must be both positive or negetive, "
                    "but now is: {}, {}, {}",
                    offsetX, offsetY, offsetZ));
  }
  m_source.unmap();
}

// if src is very huge, memmap it
void create_by_sharing_header(const std::string &segy_name,
                              const std::string &header_segy,
                              const std::string &src_file, int sizeX, int sizeY,
                              int sizeZ, int iline, int xline, int istep,
                              int xstep, int offsetX, int offsetY, int offsetZ,
                              const std::vector<std::string> &custom_info) {
  mio::mmap_source m_src;
  std::error_code error;
  m_src.map(src_file, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }

  const float *src = reinterpret_cast<const float *>(m_src.data());
  create_by_sharing_header(segy_name, header_segy, src, sizeX, sizeY, sizeZ,
                           iline, xline, istep, xstep, offsetX, offsetY,
                           offsetZ, custom_info);

  m_src.unmap();
}

void load_prestack3D(float *dst, const std::string &segy_name, int sizeX,
                     int min_inline, int max_inline, int min_crossline,
                     int max_crossline, int min_offset, int max_offset,
                     int inline_step, int crossline_step, int offset_step,
                     int inline_field, int crossline_field, int offset_field,
                     float fill) {
  float *start = dst;

  std::error_code error;
  mio::mmap_source m_source;
  m_source.map(segy_name, error);
  if (error) {
    throw std::runtime_error("Cannot mmap segy file");
  }
  const char *src = m_source.data() + kTextualHeaderSize + kBinaryHeaderSize;
  int trace_size = kTraceHeaderSize + sizeX * sizeof(float);
  int data_format =
      swap_endian<int16_t>(m_source.data() + kTextualHeaderSize + 24);

  int sizeZ = (max_inline - min_inline) / inline_step + 1;
  int step = sizeZ >= 100 ? 1 : 100 / sizeZ + 1;
  int nbar = sizeZ >= 100 ? sizeZ : step * sizeZ;
  progressbar bar(nbar);

  int iline, xline, offset;
  uint64_t trace_idx = 0;
  for (int i = min_inline; i <= max_inline; i += inline_step) {
    if (showpbar) {
      updatebar(bar, step);
    }

    for (int j = min_crossline; j <= max_crossline; j += crossline_step) {
      for (int k = min_offset; k <= max_offset; k += offset_step) {
        iline = swap_endian<int32_t>(src + trace_idx * trace_size +
                                     inline_field - 1);
        xline = swap_endian<int32_t>(src + trace_idx * trace_size +
                                     crossline_field - 1);
        offset = swap_endian<int32_t>(src + trace_idx * trace_size +
                                      offset_field - 1);

        if (i == iline && j == xline && k == offset) {
          memcpy(dst, src + trace_idx * trace_size + kTraceHeaderSize,
                 sizeX * sizeof(float));

          for (int iX = 0; iX < sizeX; iX++) {
            if (data_format == 1) {
              dst[iX] = ibm_to_ieee(dst[iX], true);
            } else if (data_format == 5) {
              dst[iX] = swap_endian(dst[iX]);
            } else {
              throw std::runtime_error("Unsuport sample format");
            }
          }

          trace_idx++;
        } else {
          std::fill(dst, dst + sizeX, fill);
        }

        dst += sizeX;
      }
    }
  }
  if (showpbar) {
    fmt::print("\n");
  }
}

} // namespace segy
