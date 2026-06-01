/*********************************************************************
** Copyright (c) 2026 Jintao Li.
** Zhejiang University (ZJU).
** All rights reserved.
*********************************************************************/

#pragma once

#include "utils.hpp"
#include <cstdio>
#include <string>

namespace segy {

class SegyBlockWriter {
public:
  SegyBlockWriter(const std::string &outname, const uchar *textual,
                  size_t textual_size, const uchar *binary,
                  size_t binary_size, const uchar *extended_textual,
                  size_t extended_textual_size, int sample_format,
                  size_t sample_count, bool overwrite);
  ~SegyBlockWriter();

  void close();
  void finalize(const uchar *data_trailer, size_t data_trailer_size);
  void write_trace_block(const uchar *trace_headers, size_t ntrace,
                         const float *samples, size_t nt);
  void write_raw_trace_block(const uchar *trace_headers, size_t ntrace,
                             const uchar *sample_bytes,
                             size_t bytes_per_trace);

  size_t trace_count() const { return m_trace_count; }
  size_t sample_count() const { return m_nt; }
  int sample_format() const { return m_dformat; }
  bool closed() const { return m_closed; }

private:
  std::string m_outname;
  FILE *m_fp = nullptr;
  int m_dformat = 0;
  size_t m_nt = 0;
  size_t m_esize = 0;
  size_t m_trace_count = 0;
  bool m_closed = false;
  bool m_finalized = false;
  WriteFunc m_wfunc;

  void set_sample_count_from_block(size_t sample_count);
  void ensure_open() const;
};

} // namespace segy
