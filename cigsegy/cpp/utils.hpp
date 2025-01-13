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
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <iostream>
#include <sstream>

#ifdef _WIN32
#define NOMINMAX
#undef max
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

using uchar = unsigned char;

namespace segy {

extern bool g_show_progress;
extern std::function<void(int current, int total)> g_progress_callback;
extern std::function<void()> g_check_signals_callback;

// const size
constexpr size_t kTextualHeaderSize = 3200;
constexpr size_t kBinaryHeaderSize = 400;
constexpr size_t kTraceHeaderStart = 3600;
constexpr size_t kTraceHeaderSize = 240;
constexpr size_t kTextualColumns = 80;
constexpr size_t kTextualRows = 40;

constexpr size_t kMaxSizeOneDimemsion = 100000;

// const binary header field
constexpr size_t kBSampleIntervalField = 17;
constexpr size_t kBSampleCountField = 21;
constexpr size_t kBSampleFormatField = 25;
constexpr size_t kBTraceSortingCodeField = 29;

// const trace header field
constexpr size_t kTStartTimeField = 105; // in ms
constexpr size_t kTDelayTimeField = 109; // in ms
constexpr size_t kTScalarField = 71;
constexpr size_t kTSampleCountField = 115;
constexpr size_t kTSampleIntervalField = 117;

constexpr size_t kDefaultInlineField = 189;
constexpr size_t kDefaultCrosslineField = 193;
constexpr size_t kDefaultXField = 73;
constexpr size_t kDefaultYField = 77;
constexpr size_t kInvalid = std::numeric_limits<size_t>::max();

// A key map that convert EBCDIC to ASCII format
const std::unordered_map<uchar, char> kEBCDICtoASCIImap = {
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

const std::unordered_map<char, uchar> kASCIItoEBCDICmap = {
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
const std::unordered_map<size_t, size_t> kElementSize = {
    {1, 4}, {2, 4}, {3, 2}, {5, 4}, {8, 1}, {10, 4}, {11, 2}, {16, 1}};

inline void swap_endian_inplace(void *dst, const void *src, size_t n) {
  uchar *_dst = static_cast<uchar *>(dst);
  const uchar *_src = static_cast<const uchar *>(src);
  static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");
  memcpy(_dst, _src, n);
  if (n == 1) {
    return;
  }
  if (n <= 8) {
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

inline float ieee_to_ibm(float value, bool is_litte_endian_input, bool is_big_endian_output=true) {
  if (!is_litte_endian_input) {
    value = swap_endian<float>(value);
  }

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

  if (is_big_endian_output) {
    return swap_endian<float>(*float_addr);
  }
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

inline float ibm_to_ieee(const void *src) {
  return ibm_to_ieee(*reinterpret_cast<const float *>(src), true);
}

inline float ieee_to_ibm(const void *src) {
  return ieee_to_ibm(*reinterpret_cast<const float *>(src), true);
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
  size_t alnumASCII = 0;
  for (size_t i = 0; i < size; i++) {
    if (std::isalnum(text[i]))
      alnumASCII++;
  }

  size_t alnumEBCDIC = 0;
  for (size_t i = 0; i < size; i++) {
    if (std::isalnum(getASCIIfromEBCDIC(text[i])))
      alnumEBCDIC++;
  }

  if (alnumASCII > alnumEBCDIC)
    return false;
  return true;
}

template <typename T>
void convert2npT(float *dst, const char *src, size_t size) {
  const T *_src = reinterpret_cast<const T *>(src);
  for (size_t i = 0; i < size; ++i) {
    dst[i] = static_cast<float>(swap_endian<T>(_src[i]));
  }
}

inline void convert2npibm(float *dst, const char *src, size_t size) {
  const float *_src = reinterpret_cast<const float *>(src);
  for (size_t i = 0; i < size; ++i) {
    dst[i] = ibm_to_ieee(_src[i], true);
  }
}

template <typename T>
void float2sgyT(char *dst, const float *src, size_t size) {
  T *_dst = reinterpret_cast<T *>(dst);

  for (size_t i = 0; i < size; ++i) {
    _dst[i] = swap_endian<T>(static_cast<T>(src[i]));
  }
}

inline void float2sgyibm(char *dst, const float *src, size_t size) {
  float *_dst = reinterpret_cast<float *>(dst);

  for (size_t i = 0; i < size; ++i) {
    _dst[i] = ieee_to_ibm(src[i], true, true);
  }
}

using ReadFunc = std::function<void(float *, const char *, size_t)>;
using WriteFunc = std::function<void(char *, const float *, size_t)>;
using ReadFuncOne = std::function<float(const char *)>;

inline void setRFunc(ReadFunc &m_readfunc, int dformat) {
  switch (dformat) {
  case 1:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      convert2npibm(dst, src, size);
    };
    break;
  case 2:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      convert2npT<int32_t>(dst, src, size);
    };
    break;
  case 3:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      convert2npT<int16_t>(dst, src, size);
    };
    break;
  case 5:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      convert2npT<float>(dst, src, size);
    };
    break;
  case 8:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      const int8_t *_src = reinterpret_cast<const int8_t *>(src);
      for (size_t i = 0; i < size; ++i) {
        dst[i] = static_cast<float>(_src[i]);
      }
    };
    break;
  case 10:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      convert2npT<uint32_t>(dst, src, size);
    };
    break;
  case 11:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      convert2npT<uint16_t>(dst, src, size);
    };
    break;
  case 16:
    m_readfunc = [](float *dst, const char *src, size_t size) {
      const uint8_t *_src = reinterpret_cast<const uint8_t *>(src);
      for (size_t i = 0; i < size; ++i) {
        dst[i] = static_cast<float>(_src[i]);
      }
    };
    break;
  default:
    throw std::invalid_argument("Unsupported dformat value: " +
                                std::to_string(dformat));
  }
}

inline void setRFuncOne(ReadFuncOne &m_readfunc, int dformat) {
  switch (dformat) {
  case 1:
    m_readfunc = [](const char *src) -> float { return ibm_to_ieee(src); };
    break;
  case 2:
    m_readfunc = [](const char *src) -> float {
      return static_cast<float>(swap_endian<int32_t>(src));
    };
    break;
  case 3:
    m_readfunc = [](const char *src) -> float {
      return static_cast<float>(swap_endian<int16_t>(src));
    };
    break;
  case 5:
    m_readfunc = [](const char *src) -> float {
      return swap_endian<float>(src);
    };
    break;
  case 8:
    m_readfunc = [](const char *src) -> float {
      return static_cast<float>(*reinterpret_cast<const int8_t *>(src));
    };
    break;
  case 10:
    m_readfunc = [](const char *src) -> float {
      return static_cast<float>(swap_endian<uint32_t>(src));
    };
    break;
  case 11:
    m_readfunc = [](const char *src) -> float {
      return static_cast<float>(swap_endian<uint16_t>(src));
    };
    break;
  case 16:
    m_readfunc = [](const char *src) -> float {
      return static_cast<float>(*reinterpret_cast<const uint8_t *>(src));
    };
    break;
  default:
    throw std::runtime_error("Unsupported dformat value: " + std::to_string(dformat));
  }
}

inline void setWFunc(WriteFunc &m_wfunc, int dformat) {
  switch (dformat) {
  case 1:
    m_wfunc = [](char *dst, const float *src, size_t size) {
      float2sgyibm(dst, src, size);
    };
    break;
  case 2:
    m_wfunc = [](char *dst, const float *src, size_t size) {
      float2sgyT<int32_t>(dst, src, size);
    };
    break;
  case 3:
    m_wfunc = [](char *dst, const float *src, size_t size) {
      float2sgyT<int16_t>(dst, src, size);
    };
    break;
  case 5:
    m_wfunc = [](char *dst, const float *src, size_t size) {
      float2sgyT<float>(dst, src, size);
    };
    break;
  case 8:
    m_wfunc = [](char *dst, const float *src, size_t size) {
      float2sgyT<int8_t>(dst, src, size);
    };
    break;
  case 10:
    m_wfunc = [](char *dst, const float *src, size_t size) {
      float2sgyT<uint32_t>(dst, src, size);
    };
    break;
  case 11:
    m_wfunc = [](char *dst, const float *src, size_t size) {
      float2sgyT<uint16_t>(dst, src, size);
    };
    break;
  case 16:
    m_wfunc = [](char *dst, const float *src, size_t size) {
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

inline void append_to_file(const std::string &file_name, int64_t size_to_add) {
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
                          int64_t size_to_remove) {
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

inline void default_progress_callback(int current, int total) {
    const int barWidth = 50;

    std::ostringstream oss; // 使用字符串流存储输出

    if (current == 0) {
        oss << "Progress: [>" << std::string(barWidth, ' ') << "] 0% (0/" << total << ")";
    } else if (current > 0 && current <= total) {
        float progress = static_cast<float>(current) / total;
        int pos = static_cast<int>(barWidth * progress);

        oss << "\rProgress: [";
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) oss << "=";
            else if (i == pos) oss << ">";
            else oss << " ";
        }
        oss << "] " << static_cast<int>(progress * 100) << "% (" << current << "/" << total << ")";
    } else if (current == -1) {
        oss << "\rProgress: [" << std::string(barWidth, '=') << ">] 100% (" << total << "/" << total << ")";
    }

    // 一次性输出
    std::cout << oss.str() << std::flush;
}

inline int cal_progress_steps(int total, bool show_progress, int min_iters = 100, int max_iters = 50) {
    if (!show_progress || !g_show_progress || total < min_iters) {
        return 0;
    }

    return total / max_iters;

}

} // namespace segy

#endif
