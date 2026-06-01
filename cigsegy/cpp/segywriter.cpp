/*********************************************************************
** Copyright (c) 2026 Jintao Li.
** Zhejiang University (ZJU).
** All rights reserved.
*********************************************************************/

#include "segywriter.h"
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace segy {

namespace {

int16_t read_i2_be(const uchar *header, size_t loc) {
  size_t start = loc - 1;
  uint16_t raw = (static_cast<uint16_t>(header[start]) << 8) |
                 static_cast<uint16_t>(header[start + 1]);
  return static_cast<int16_t>(raw);
}

void write_i2_be(uchar *header, size_t loc, size_t value) {
  if (value > static_cast<size_t>(std::numeric_limits<int16_t>::max())) {
    throw std::runtime_error("SEG-Y 2-byte header value is out of range: " +
                             std::to_string(value));
  }
  size_t start = loc - 1;
  uint16_t raw = static_cast<uint16_t>(value);
  header[start] = static_cast<uchar>((raw >> 8) & 0xffU);
  header[start + 1] = static_cast<uchar>(raw & 0xffU);
}

void write_all(FILE *fp, const void *data, size_t nbytes) {
  const char *ptr = static_cast<const char *>(data);
  size_t written = 0;
  while (written < nbytes) {
    size_t n = std::fwrite(ptr + written, 1, nbytes - written, fp);
    if (n == 0) {
      throw std::runtime_error("failed to write SEG-Y file");
    }
    written += n;
  }
}

bool file_exists(const std::string &path) {
  FILE *fp = std::fopen(path.c_str(), "rb");
  if (fp != nullptr) {
    std::fclose(fp);
    return true;
  }
  return false;
}

} // namespace

SegyBlockWriter::SegyBlockWriter(
    const std::string &outname, const uchar *textual, size_t textual_size,
    const uchar *binary, size_t binary_size, const uchar *extended_textual,
    size_t extended_textual_size, int sample_format, size_t sample_count,
    bool overwrite)
    : m_outname(outname) {
  if (textual_size != kTextualHeaderSize) {
    throw std::runtime_error("textual header size must be 3200, but got " +
                             std::to_string(textual_size));
  }
  if (binary_size != kBinaryHeaderSize) {
    throw std::runtime_error("binary header size must be 400, but got " +
                             std::to_string(binary_size));
  }
  if (extended_textual_size % kTextualHeaderSize != 0) {
    throw std::runtime_error(
        "extended textual header size must be a multiple of 3200");
  }
  if (!overwrite && file_exists(outname)) {
    throw std::runtime_error("output SEG-Y file already exists: " + outname);
  }

  m_dformat = sample_format > 0
                  ? sample_format
                  : static_cast<int>(read_i2_be(binary, kBSampleFormatField));
  m_nt = sample_count > 0
             ? sample_count
             : static_cast<size_t>(read_i2_be(binary, kBSampleCountField));

  auto it = kElementSize.find(m_dformat);
  if (it == kElementSize.end()) {
    throw std::runtime_error("Unknown data format: " +
                             std::to_string(m_dformat));
  }
  m_esize = it->second;
  setWFunc(m_wfunc, m_dformat);

  std::vector<uchar> binary_bytes(binary, binary + binary_size);
  if (sample_format > 0) {
    write_i2_be(binary_bytes.data(), kBSampleFormatField,
                static_cast<size_t>(sample_format));
  }
  if (sample_count > 0) {
    write_i2_be(binary_bytes.data(), kBSampleCountField, sample_count);
  }

  m_fp = std::fopen(outname.c_str(), "wb");
  if (m_fp == nullptr) {
    throw std::runtime_error("failed to open SEG-Y output: " + outname);
  }
  write_all(m_fp, textual, textual_size);
  write_all(m_fp, binary_bytes.data(), binary_bytes.size());
  if (extended_textual_size > 0) {
    write_all(m_fp, extended_textual, extended_textual_size);
  }
}

SegyBlockWriter::~SegyBlockWriter() {
  try {
    close();
  } catch (...) {
  }
}

void SegyBlockWriter::close() {
  if (m_closed) {
    return;
  }
  if (m_fp != nullptr) {
    if (std::fclose(m_fp) != 0) {
      m_fp = nullptr;
      m_closed = true;
      throw std::runtime_error("failed to close SEG-Y output: " + m_outname);
    }
    m_fp = nullptr;
  }
  m_closed = true;
}

void SegyBlockWriter::finalize(const uchar *data_trailer,
                               size_t data_trailer_size) {
  if (m_finalized) {
    return;
  }
  ensure_open();
  if (data_trailer_size > 0) {
    write_all(m_fp, data_trailer, data_trailer_size);
  }
  if (std::fflush(m_fp) != 0) {
    throw std::runtime_error("failed to flush SEG-Y output: " + m_outname);
  }
  m_finalized = true;
  close();
}

void SegyBlockWriter::write_trace_block(const uchar *trace_headers,
                                        size_t ntrace, const float *samples,
                                        size_t nt) {
  ensure_open();
  if (ntrace == 0) {
    return;
  }
  if (m_nt == 0) {
    set_sample_count_from_block(nt);
  }
  if (nt != m_nt) {
    throw std::runtime_error("samples shape[1] does not match sample_count");
  }

  size_t trace_size = kTraceHeaderSize + m_nt * m_esize;
  std::vector<char> buffer(ntrace * trace_size);
  char *dst = buffer.data();
  for (size_t i = 0; i < ntrace; ++i) {
    memcpy(dst, trace_headers + i * kTraceHeaderSize, kTraceHeaderSize);
    dst += kTraceHeaderSize;
    m_wfunc(dst, samples + i * m_nt, m_nt);
    dst += m_nt * m_esize;
  }
  write_all(m_fp, buffer.data(), buffer.size());
  m_trace_count += ntrace;
}

void SegyBlockWriter::write_raw_trace_block(const uchar *trace_headers,
                                            size_t ntrace,
                                            const uchar *sample_bytes,
                                            size_t bytes_per_trace) {
  ensure_open();
  if (ntrace == 0) {
    return;
  }
  if (m_nt == 0) {
    if (bytes_per_trace % m_esize != 0) {
      throw std::runtime_error(
          "sample bytes per trace is not divisible by sample element size");
    }
    set_sample_count_from_block(bytes_per_trace / m_esize);
  }
  size_t expected = m_nt * m_esize;
  if (bytes_per_trace != expected) {
    throw std::runtime_error("sample bytes per trace does not match "
                             "sample_count * sample element size");
  }

  size_t trace_size = kTraceHeaderSize + expected;
  std::vector<char> buffer(ntrace * trace_size);
  char *dst = buffer.data();
  for (size_t i = 0; i < ntrace; ++i) {
    memcpy(dst, trace_headers + i * kTraceHeaderSize, kTraceHeaderSize);
    dst += kTraceHeaderSize;
    memcpy(dst, sample_bytes + i * expected, expected);
    dst += expected;
  }
  write_all(m_fp, buffer.data(), buffer.size());
  m_trace_count += ntrace;
}

void SegyBlockWriter::set_sample_count_from_block(size_t sample_count) {
  if (sample_count == 0) {
    throw std::runtime_error("sample_count must be positive");
  }
  m_nt = sample_count;
  uchar bytes[2] = {0, 0};
  write_i2_be(bytes, 1, sample_count);
  long current = std::ftell(m_fp);
  if (current < 0) {
    throw std::runtime_error("failed to query SEG-Y output position");
  }
  long offset = static_cast<long>(kTextualHeaderSize + kBSampleCountField - 1);
  if (std::fseek(m_fp, offset, SEEK_SET) != 0) {
    throw std::runtime_error("failed to seek SEG-Y binary header");
  }
  write_all(m_fp, bytes, sizeof(bytes));
  if (std::fseek(m_fp, current, SEEK_SET) != 0) {
    throw std::runtime_error("failed to restore SEG-Y output position");
  }
}

void SegyBlockWriter::ensure_open() const {
  if (m_closed || m_fp == nullptr) {
    throw std::runtime_error("SEG-Y block writer is closed");
  }
  if (m_finalized) {
    throw std::runtime_error("SEG-Y block writer is finalized");
  }
}

} // namespace segy
