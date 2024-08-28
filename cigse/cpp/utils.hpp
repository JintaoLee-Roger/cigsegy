/*********************************************************************
** Copyright (c) 2024 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
*********************************************************************/

#ifndef CIG_UTILS_H
#define CIG_UTILS_H

#include <algorithm>
#include <cctype>
#include <climits>
#include <cstring>
#include <limits>
#include <map>
#include <string.h>
#include <utility>

#ifdef _WIN32
#define NOMINMAX
#undef max
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

#ifdef USE_PYBIND11
#include <pybind11/pybind11.h>

inline void checkSignals() {
  if (PyErr_CheckSignals() != 0) {
    throw pybind11::error_already_set();
  }
}
#define CHECK_SIGNALS() checkSignals()

#else
#define CHECK_SIGNALS()                                                        \
  do {                                                                         \
  } while (0)
#endif

using uchar = unsigned char;

namespace segy {

// const size
const int kTextualHeaderSize = 3200;
const int kBinaryHeaderSize = 400;
const int kTraceHeaderStart = 3600;
const int kTraceHeaderSize = 240;
const int kTextualColumns = 80;
const int kTextualRows = 40;

const int kMaxSizeOneDimemsion = 100000;

// const binary header field
const int kBSampleIntervalField = 17;
const int kBSampleCountField = 21;
const int kBSampleFormatField = 25;
const int kBTraceSortingCodeField = 29;

// const trace header field
const int kTStartTimeField = 105; // in ms
const int kTDelayTimeField = 109; // in ms
const int kTScalarField = 71;
const int kTSampleCountField = 115;
const int kTSampleIntervalField = 117;

const int kDefaultInlineField = 189;
const int kDefaultCrosslineField = 193;
const int kDefaultXField = 73;
const int kDefaultYField = 77;

// A key map that convert EBCDIC to ASCII format
const std::map<unsigned char, char> kEBCDICtoASCIImap = {
    {64, ' '},   {75, '.'},  {76, '<'},   {77, '('},  {78, '+'},  {79, '|'},
    {80, '&'},   {90, '!'},  {91, '$'},   {92, '*'},  {93, ')'},  {94, ';'},
    {96, '-'},   {97, '/'},  {106, '|'},  {107, ','}, {108, '%'}, {109, '_'},
    {110, '>'},  {111, '?'}, {121, '`'},  {122, ':'}, {123, '#'}, {124, '@'},
    {125, '\''}, {126, '='}, {127, '"'},  {129, 'a'}, {130, 'b'}, {131, 'c'},
    {132, 'd'},  {133, 'e'}, {134, 'f'},  {135, 'g'}, {136, 'h'}, {137, 'i'},
    {145, 'j'},  {146, 'k'}, {147, 'l'},  {148, 'm'}, {149, 'n'}, {150, 'o'},
    {151, 'p'},  {152, 'q'}, {153, 'r'},  {161, '~'}, {162, 's'}, {163, 't'},
    {164, 'u'},  {165, 'v'}, {166, 'w'},  {167, 'x'}, {168, 'y'}, {169, 'z'},
    {192, '{'},  {193, 'A'}, {194, 'B'},  {195, 'C'}, {196, 'D'}, {197, 'E'},
    {198, 'F'},  {199, 'G'}, {200, 'H'},  {201, 'I'}, {208, '}'}, {209, 'J'},
    {210, 'K'},  {211, 'L'}, {212, 'M'},  {213, 'N'}, {214, 'O'}, {215, 'P'},
    {216, 'Q'},  {217, 'R'}, {224, '\\'}, {226, 'S'}, {227, 'T'}, {228, 'U'},
    {229, 'V'},  {230, 'W'}, {231, 'X'},  {232, 'Y'}, {233, 'Z'}, {240, '0'},
    {241, '1'},  {242, '2'}, {243, '3'},  {244, '4'}, {245, '5'}, {246, '6'},
    {247, '7'},  {248, '8'}, {249, '9'}};

const std::map<char, unsigned char> kASCIItoEBCDICmap = {
    {' ', 64},   {'.', 75},  {'<', 76},   {'(', 77},  {'+', 78},  {'|', 79},
    {'&', 80},   {'!', 90},  {'$', 91},   {'*', 92},  {')', 93},  {';', 94},
    {'-', 96},   {'/', 97},  {'|', 106},  {',', 107}, {'%', 108}, {'_', 109},
    {'>', 110},  {'?', 111}, {'`', 121},  {':', 122}, {'#', 123}, {'@', 124},
    {'\'', 125}, {'=', 126}, {'"', 127},  {'a', 129}, {'b', 130}, {'c', 131},
    {'d', 132},  {'e', 133}, {'f', 134},  {'g', 135}, {'h', 136}, {'i', 137},
    {'j', 145},  {'k', 146}, {'l', 147},  {'m', 148}, {'n', 149}, {'o', 150},
    {'p', 151},  {'q', 152}, {'r', 153},  {'~', 161}, {'s', 162}, {'t', 163},
    {'u', 164},  {'v', 165}, {'w', 166},  {'x', 167}, {'y', 168}, {'z', 169},
    {'{', 192},  {'A', 193}, {'B', 194},  {'C', 195}, {'D', 196}, {'E', 197},
    {'F', 198},  {'G', 199}, {'H', 200},  {'I', 201}, {'}', 208}, {'J', 209},
    {'K', 210},  {'L', 211}, {'M', 212},  {'N', 213}, {'O', 214}, {'P', 215},
    {'Q', 216},  {'R', 217}, {'\\', 224}, {'S', 226}, {'T', 227}, {'U', 228},
    {'V', 229},  {'W', 230}, {'X', 231},  {'Y', 232}, {'Z', 233}, {'0', 240},
    {'1', 241},  {'2', 242}, {'3', 243},  {'4', 244}, {'5', 245}, {'6', 246},
    {'7', 247},  {'8', 248}, {'9', 249}};

// NOTE: only support 1, 2, 4 bytes
const std::map<int, int> kElementSize = {{1, 4}, {2, 4},  {3, 2},  {5, 4},
                                         {8, 1}, {10, 4}, {11, 2}, {16, 1}};

inline void swap_endian_inplace(void *dst, const void *src, int n) {
  uchar *_dst = static_cast<uchar *>(dst);
  const uchar *_src = static_cast<const uchar *>(src);
  static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");
  memcpy(_dst, _src, n);
  if (n == 1) {
    return;
  }
  if (n <= 8) {
    // don't do anything if n > 16
    // uchar u8[17] = {0}; // for src == dst
    // for (int i = 0; i < n; i++) {
    //   u8[n - i - 1] = _src[i];
    // }
    // memcpy(_dst, u8, n);
    std::reverse(_dst, _dst + n);
  }
}

template <typename T> T swap_endian(const void *src) {
  T u = 0;
  swap_endian_inplace(&u, src, sizeof(T));
  return u;
}

template <typename T> T swap_endian(T u) {
  static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");

  union {
    T u;
    unsigned char u8[sizeof(T)];
  } source, dest;

  source.u = u;

  for (size_t k = 0; k < sizeof(T); k++)
    dest.u8[k] = source.u8[sizeof(T) - k - 1];

  return dest.u;
}

inline float ieee_to_ibm(float value, bool is_litte_endian_input) {
  if (!is_litte_endian_input)
    value = swap_endian<float>(value);

  int32_t *addr = reinterpret_cast<int32_t *>(&value);
  int32_t int_val = *addr;

  int32_t sign = (int_val >> 31) & 1;
  int32_t exponent = ((int_val & 0x7f800000) >> 23) - 127;
  int32_t fraction = int_val & 0x007fffff;

  if ((int_val & 0x7fffffff) == 0) {
    return sign ? -0.0f : 0.0f;
  }

  fraction <<= 1; // 24 bits

  fraction |= 0x01000000; // add 1, 25 bits

  // convert 2-base to 16-base
  fraction <<= (exponent & 3); // 28 bits
  exponent >>= 2;

  if (fraction & 0x0f000000) { // 24 bits
    fraction >>= 4;
    exponent += 1;
  }

  exponent += 64;

  int32_t ibm_value;
  if (exponent > 127) {
    return (sign ? -std::numeric_limits<float>::max()
                 : std::numeric_limits<float>::max());
  } else if (exponent <= 0) {
    ibm_value = (sign << 31) | fraction;
  } else {
    ibm_value = (sign << 31) | (exponent << 24) | fraction;
  }

  float *float_addr = reinterpret_cast<float *>(&ibm_value);

  return *float_addr;
}

inline float ibm_to_ieee(float value, bool is_big_endian_input) {
  if (is_big_endian_input) {
    value = swap_endian<float>(value);
  }

  int32_t *int_addr = reinterpret_cast<int32_t *>(&value);
  int32_t int_val = *int_addr;

  int32_t sign = (int_val >> 31) & 1;
  int32_t fraction = int_val & 0x00ffffff;

  if (fraction == 0) {
    return sign ? -0.0f : 0.0f;
  }

  // Convert exponent to be of base 2 and remove IBM exponent bias.
  int32_t exponent = ((int_val & 0x7f000000) >> 22) - 256;

  // Drop the last bit since we can store only 23 bits in IEEE.
  fraction >>= 1;

  // Normalize such that the implicit leading bit of the fraction is 1.
  while (fraction && (fraction & 0x00800000) == 0) {
    fraction <<= 1;
    --exponent;
  }

  // Drop the implicit leading bit.
  fraction &= 0x007fffff;

  // Add IEEE bias to the exponent.
  exponent += 127;

  // Handle overflow.
  if (exponent >= 255) {
    return (sign ? -std::numeric_limits<float>::max()
                 : std::numeric_limits<float>::max());
  }

  int32_t ieee_value;

  // Handle underflow.
  if (exponent <= 0)
    ieee_value = (sign << 31) | fraction;
  else
    ieee_value = (sign << 31) | (exponent << 23) | fraction;

  float *float_addr = reinterpret_cast<float *>(&ieee_value);
  return *float_addr;
}

inline char getASCIIfromEBCDIC(char c) {
  if (kEBCDICtoASCIImap.find(c) != kEBCDICtoASCIImap.end())
    return kEBCDICtoASCIImap.at(c);
  return ' ';
}

inline char getEBCIDfromASCII(char c) {
  if (kASCIItoEBCDICmap.find(c) != kASCIItoEBCDICmap.end()) {
    return kASCIItoEBCDICmap.at(c);
  }
  return ' ';
}

inline bool isTextInEBCDICFormat(const char *text, size_t size) {
  int alnumASCII = 0;
  for (size_t i = 0; i < size; i++) {
    if (std::isalnum(text[i]))
      alnumASCII++;
  }

  int alnumEBCDIC = 0;
  for (size_t i = 0; i < size; i++) {
    if (std::isalnum(getASCIIfromEBCDIC(text[i])))
      alnumEBCDIC++;
  }

  if (alnumASCII > alnumEBCDIC)
    return false;
  return true;
}

template <typename T> void convert2npT(float *dst, const char *src, int size) {
  const T *_src = reinterpret_cast<const T *>(src);
  for (int i = 0; i < size; ++i) {
    dst[i] = ibm_to_ieee(_src[i], true);
  }
}

inline void convert2npibm(float *dst, const char *src, int size) {
  const float *_src = reinterpret_cast<const float *>(src);
  for (int i = 0; i < size; ++i) {
    dst[i] = ibm_to_ieee(_src[i], true);
  }
}

template <typename T> void float2sgyT(char *dst, const float *src, int size) {
  float *_dst = reinterpret_cast<float *>(dst);

  for (int i = 0; i < size; ++i) {
    _dst[i] = swap_endian<T>(T(src[i]));
  }
}

inline void float2sgyibm(char *dst, const float *src, int size) {
  float *_dst = reinterpret_cast<float *>(dst);

  for (int i = 0; i < size; ++i) {
    _dst[i] = ieee_to_ibm(src[i], true);
  }
}

using ReadFunc = std::function<void(float *, const char *, int)>;
using WriteFunc = std::function<void(char *, const float *, int)>;

inline void setRFunc(ReadFunc &m_readfunc, int dformat) {
  switch (dformat) {
  case 1:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npibm(dst, src, size);
    };
    break;
  case 2:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<int32_t>(dst, src, size);
    };
    break;
  case 3:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<int16_t>(dst, src, size);
    };
    break;
  case 5:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<float>(dst, src, size);
    };
    break;
  case 8:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<int8_t>(dst, src, size);
    };
    break;
  case 10:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<uint32_t>(dst, src, size);
    };
    break;
  case 11:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<uint16_t>(dst, src, size);
    };
    break;
  case 16:
    m_readfunc = [](float *dst, const char *src, int size) {
      convert2npT<uint8_t>(dst, src, size);
    };
    break;
  default:
    throw std::invalid_argument("Unsupported dformat value: " +
                                std::to_string(dformat));
  }
}

inline void setWFunc(WriteFunc &m_wfunc, int dformat) {
  switch (dformat) {
  case 1:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyibm(dst, src, size);
    };
    break;
  case 2:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<int32_t>(dst, src, size);
    };
    break;
  case 3:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<int16_t>(dst, src, size);
    };
    break;
  case 5:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<float>(dst, src, size);
    };
    break;
  case 8:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<int8_t>(dst, src, size);
    };
    break;
  case 10:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<uint32_t>(dst, src, size);
    };
    break;
  case 11:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<uint16_t>(dst, src, size);
    };
    break;
  case 16:
    m_wfunc = [](char *dst, const float *src, int size) {
      float2sgyT<uint8_t>(dst, src, size);
    };
    break;
  default:
    throw std::invalid_argument("Unsupported dformat value: " +
                                std::to_string(dformat));
  }
}

// ofstream is very slow on out disk, so we use low-level I/O functions
inline void create_file(const std::string &file_name, uint64_t file_size) {
#ifdef _WIN32
  HANDLE file = CreateFile(file_name.c_str(), GENERIC_WRITE, 0, NULL,
                           CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
  if (file == INVALID_HANDLE_VALUE) {
    throw std::runtime_error("Failed to create file");
  }

  LARGE_INTEGER size;
  size.QuadPart = file_size;
  if (SetFilePointerEx(file, size, NULL, FILE_BEGIN) == 0 ||
      SetEndOfFile(file) == 0) {
    CloseHandle(file);
    throw std::runtime_error("Failed to set file size");
  }

  CloseHandle(file);
#else
  int fd = open(file_name.c_str(), O_CREAT | O_RDWR, 0666);
  if (fd < 0) {
    throw std::runtime_error("Failed to create file");
  }

  if (ftruncate(fd, file_size) != 0) {
    close(fd);
    throw std::runtime_error("Failed to set file size");
  }

  close(fd);
#endif
}

inline void append_to_file(const std::string &file_name, uint64_t size_to_add) {
#ifdef _WIN32
  HANDLE file = CreateFile(file_name.c_str(), GENERIC_WRITE, 0, NULL,
                           OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (file == INVALID_HANDLE_VALUE) {
    throw std::runtime_error("Failed to open file");
  }

  LARGE_INTEGER current_size;
  if (!GetFileSizeEx(file, &current_size)) {
    CloseHandle(file);
    throw std::runtime_error("Failed to get file size");
  }

  LARGE_INTEGER new_size;
  new_size.QuadPart = current_size.QuadPart + size_to_add;

  if (SetFilePointerEx(file, new_size, NULL, FILE_BEGIN) == 0 ||
      SetEndOfFile(file) == 0) {
    CloseHandle(file);
    throw std::runtime_error("Failed to append to file");
  }

  CloseHandle(file);
#else
  int fd = open(file_name.c_str(), O_RDWR);
  if (fd < 0) {
    throw std::runtime_error("Failed to open file");
  }

  off_t current_size = lseek(fd, 0, SEEK_END);
  if (current_size == (off_t)-1) {
    close(fd);
    throw std::runtime_error("Failed to get file size");
  }

  off_t new_size = current_size + size_to_add;
  if (ftruncate(fd, new_size) != 0) {
    close(fd);
    throw std::runtime_error("Failed to append to file");
  }

  close(fd);
#endif
}

inline void truncate_file(const std::string &file_name,
                          uint64_t size_to_remove) {
#ifdef _WIN32
  HANDLE file = CreateFile(file_name.c_str(), GENERIC_WRITE, 0, NULL,
                           OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (file == INVALID_HANDLE_VALUE) {
    throw std::runtime_error("Failed to open file");
  }

  LARGE_INTEGER current_size;
  if (!GetFileSizeEx(file, &current_size)) {
    CloseHandle(file);
    throw std::runtime_error("Failed to get file size");
  }

  LARGE_INTEGER new_size;
  new_size.QuadPart = current_size.QuadPart > size_to_remove
                          ? current_size.QuadPart - size_to_remove
                          : 0;

  if (SetFilePointerEx(file, new_size, NULL, FILE_BEGIN) == 0 ||
      SetEndOfFile(file) == 0) {
    CloseHandle(file);
    throw std::runtime_error("Failed to truncate file");
  }

  CloseHandle(file);
#else
  int fd = open(file_name.c_str(), O_RDWR);
  if (fd < 0) {
    throw std::runtime_error("Failed to open file");
  }

  off_t current_size = lseek(fd, 0, SEEK_END);
  if (current_size == (off_t)-1) {
    close(fd);
    throw std::runtime_error("Failed to get file size");
  }

  off_t new_size =
      current_size > size_to_remove ? current_size - size_to_remove : 0;

  if (ftruncate(fd, new_size) != 0) {
    close(fd);
    throw std::runtime_error("Failed to truncate file");
  }

  close(fd);
#endif
}

} // namespace segy

#endif
