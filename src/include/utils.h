/*********************************************************************
** Copyright (c) 2023 Roger Lee.
** Computational and Interpretation Group (CIG),
** University of Science and Technology of China (USTC).
**
** @File: utils.h
** @Description :
*********************************************************************/
#ifndef CIG_UTILS_H
#define CIG_UTILS_H

#include <algorithm>
#include <cctype>
#include <climits>
#include <limits>
#include <map>
#include <string.h>
#include <utility>

#ifdef _WIN32
#undef max
#endif

using uchar = unsigned char;

namespace segy {

// const size
const int kTextualHeaderSize = 3200;
const int kBinaryHeaderSize = 400;
const int kTraceHeaderSize = 240;
const int kTextualColumns = 80;
const int kTextualRows = 40;

const int kMaxSizeOneDimemsion = 100000;

// const binary header field
const int kBSampleIntervalField = 17;
const int kBSampleCountField = 21;
const int kBSampleFormatField = 25;

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

// const int kMaxTempSize = 512 * 512 * 512 * 4;

// const int kMaxThreadsNum = 8;

#pragma pack(push, 1)
struct BinaryHeader {
  int32_t jobID;                       // 1-4
  int32_t line_number;                 // 5-8
  int32_t reel_number;                 // 9-12
  int16_t num_traces_per_ensemble;     // 13-14
  int16_t num_aux_traces_per_ensemble; // 15-16
  int16_t sample_interval;             // 17-18
  int16_t sample_interval_orig;        // 19-20
  int16_t trace_length;                // 21-22
  int16_t trace_length_orig;           // 23-24
  int16_t data_format;                 // 25-26
  int16_t ensemble_fold;               // 27-28
  int16_t trace_sorting_code;          // 29-30
  int16_t v_sum_code;                  // 31-32
  int16_t sweep_freq_start;            // 33-34
  int16_t sweep_freq_end;              // 35-36
  int16_t sweep_length;                // 37-38
  int16_t sweep_type_code;             // 39-40
  int16_t trace_num_sweep_channel;     // 41-42
  int16_t sweep_trace_taper_start;     // 43-44
  int16_t sweep_trace_taper_end;       // 45-46
  int16_t taper_type;                  // 47-48
  int16_t correlated_data_trace;       // 49-50
  int16_t bin_gain_recover;            // 51-52
  int16_t amplitude_recover_method;    // 53-54
  int16_t measurement_system;          // 55-56
  int16_t impulse_signal_polarity;     // 57-58
  int16_t vibratory_polarity_code;     // 59-60
  int32_t extend_num_data_tarces;      // 61-64
  int32_t extend_num_aux_data_tarces;  // 65-68
  int32_t extend_trace_length;         // 69-72
  double extend_sample_intervel;       // 73-80
  double extend_sample_intervel_orig;  // 81-88
  int32_t extend_trace_length_orig;    // 89-92
  int32_t extend_ensement_fold;        // 93-96
  int32_t const_1;                     // 97-100
  uchar dummy[200];                    // 101-300
  unsigned char major_version;         // 301
  unsigned char minor_version;         // 302
  int16_t fixed_length_trace;          // 303-304
  int16_t extend_textual_header;       // 305-306
  int32_t max_extend_trace_header;     // 307-310
  int16_t time_bias_code;              // 311-312
  uint64_t num_traces;                 // 313-320
  uint64_t byte_offset;                // 321-328
  int32_t num_trailer_stanza;          // 329-332
  uchar dummy2[68];                    // 333-400
};

struct TraceHeader {
  int32_t trace_sequence_number_in_line; // 1-4
  int32_t trace_sequence_number_in_file; // 5-8
  int32_t orig_field_num;                // 9-12
  int32_t trace_num_in_orig;             // 13-16
  int32_t source_point_num;              // 17-20
  int32_t ensemble_num;                  // 21-24
  int32_t trace_num_in_ensemble;         // 25-28
  int16_t trace_ID_code;                 // 29-30
  int16_t num_v_summed_traces;           // 31-32
  int16_t num_h_stacked_tarces;          // 33-34
  int16_t data_used_for;                 // 35-36
  int32_t distance_from_center;          // 37-40
  int32_t elevation_rev;                 // 41-44
  int32_t surface_elevation_source;      // 45-48
  int32_t source_depth;                  // 49-52
  int32_t seis_datum_elevation_rev;      // 53-56
  int32_t seis_datum_elevation_source;   // 57-60
  int32_t water_col_height_source;       // 61-64
  int32_t water_col_height_rev;          // 65-68
  int16_t scalar_for_elev_and_depth;     // 69-70
  int16_t scalar_for_coord;              // 71-72
  int32_t source_coord_X;                // 73-76
  int32_t source_coord_Y;                // 77-80
  int32_t group_coord_X;                 // 81-84
  int32_t group_coord_Y;                 // 85-88
  int16_t coord_units;                   // 89-90
  int16_t weather_vel;                   // 91-92
  int16_t subweather_vel;                // 93-94
  int16_t uphole_time_source;            // 95-96
  int16_t uphole_time_rev;               // 97-98
  int16_t source_static_corr;            // 99-100
  int16_t group_static_corr;             // 101-102
  int16_t total_static;                  // 103-104
  int16_t lag_time_A;                    // 105-106
  int16_t lag_time_B;                    // 107-108
  int16_t delay_record_time;             // 109-110
  int16_t mute_time_start;               // 111-112
  int16_t mute_time_end;                 // 113-114
  int16_t num_sample;                    // 115-116
  int16_t sample_interval;               // 117-118
  int16_t gain_type;                     // 119-120
  int16_t instrument_gain_constant;      // 121-122
  int16_t instrument_early;              // 123-124
  int16_t correlated;                    // 125-126
  int16_t sweep_freq_start;              // 127-128
  int16_t sweep_freq_end;                // 129-130
  int16_t sweep_length;                  // 131-132
  int16_t sweep_type_code;               // 133-134
  int16_t sweep_trace_taper_start;       // 135-136
  int16_t sweep_trace_taper_end;         // 137-138
  int16_t taper_type;                    // 139-140
  int16_t alias_filter_freq;             // 141-142
  int16_t alias_filter_slope;            // 143-144
  int16_t notch_filter_freq;             // 145-146
  int16_t notch_filter_slope;            // 147-148
  int16_t lowcut_freq;                   // 149-150
  int16_t highcut_freq;                  // 151-152
  int16_t lowcut_scope;                  // 153-154
  int16_t highcut_scope;                 // 155-156
  int16_t years;                         // 157-158
  int16_t day;                           // 159-160
  int16_t hour;                          // 161-162
  int16_t minute;                        // 163-164
  int16_t secend;                        // 165-166
  int16_t time_basis;                    // 167-168
  int16_t trace_weight_factor;           // 169-170
  int16_t geophone[3];                   // 171-176
  int16_t gap_size;                      // 177-178
  int16_t down_or_up;                    // 179-180
  int32_t X;                             // 181-184
  int32_t Y;                             // 185-188
  int32_t inline_num;                    // 189-192
  int32_t crossline_num;                 // 193-196
  int32_t shotpoint_num;                 // 197-200
  int16_t scalar_for_shotpoint;          // 201-202
  int16_t trace_value_measurement_uint;  // 203-204
  uchar transduction[6];                 // 205-210
  int16_t transduction_unit;             // 211-212
  int16_t traceID;                       // 213-214
  int16_t scalar_for_95;                 // 215-216
  int16_t source_type;                   // 217-218
  uchar source_energy_direction[6];      // 219-224
  uchar source_measurement[6];           // 225-230
  int16_t source_measu_unit;             // 231-232
  int64_t dummy;                         // 233-240
};
#pragma pack(pop)

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

const std::map<int, std::pair<const char *, int>> kBinaryHeaderHelp = {
    {1, {"Job ID", 4}},
    {5, {"Line number", 4}},
    {9, {"Reel Number", 4}},
    {13, {"N traces per ensembel", 2}},
    {15, {"N auxiliary traces per ensembel", 2}},
    {17, {"Sample interval(dt)", 2}},
    {19, {"dt of original", 2}},
    {21, {"N samples per traces(ns)", 2}},
    {23, {"ns of orignal", 2}},
    {25, {"Data sample format code (1-IBM, 5-IEEE)", 2}},
    {27, {"Ensemble fold", 2}},
    {29, {"Trace sorting code", 2}},
    {31, {"vertical sum code", 2}},
    {33, {"Sweep freq at start(Hz)", 2}},
    {35, {"Sweep freq at end(HZ)", 2}},
    {37, {"Sweep length(ms)", 2}},
    {39, {"Sweep type code", 2}},
    {41, {"Trace number of sweep channel", 2}},
    {43, {"Sweep trace taper length in ms at strat", 2}},
    {45, {"Sweep trace taper length in ms at end", 2}},
    {47, {"Taper type", 2}},
    {49, {"Correlated data traces", 2}},
    {51, {"Binary gain recovered", 2}},
    {53, {"Amplitude recovery method", 2}},
    {55, {"Measurement system (units)", 2}},
    {57, {"Impulse signal polarity", 2}},
    {59, {"Vibratory polariy code", 2}},
    {61, {"Extended number of data traces per ensemble", 4}},
    {65, {"Extended number of auxiliary traces per ensemble", 4}},
    {69, {"Extended number of samples per data trace", 4}},
    {73, {"Extended sample interval, IEEE double precision (64-bit)", 9}},
    {81,
     {"Extended sample interval of original field recording, IEEE double "
      "precision (64-bit)",
      9}},
    {89,
     {"Extended number of samples per data trace in original recording", 4}},
    {93, {"Extended ensemble fold", 4}},
    {97, {"The integer constant 16909060_10 (01020304_16)", 4}},
    {101, {"Unassigned", 200}},
    {301,
     {"Major SEG-Y Format Revision Number. This is an 8-bit unsigned value",
      1}},
    {302,
     {"Minor SEG-Y Format Revision Number. This is an 8-bit unsigned value "
      "with a radix point between the first and second bytes.",
      1}},
    {303, {"Fixed length trace flag", 2}},
    {305,
     {"Number of 3200-byte, Extended Textual File Header records following "
      "the Binary Header",
      2}},
    {307, {"Max number of additional 240-byte trace header", 4}},
    {311, {"Time basis code", 2}},
    {313,
     {"number of trace header in this file, 64-bit unsigned integer value", 9}},
    {321,
     {"Byte offset of first trace relative to start of file or stream if "
      "known, otherwise zero. (64-bit unsigned integer value) ",
      9}},
    {329,
     {"Number of 3200-byte data trailer stanza records following the last "
      "trace (4 byte signed integer).",
      4}},
    {333, {"Unassigned", 68}}};

const std::map<int, std::pair<const char *, int>> kTraceHeaderHelp = {
    {1, {"Trace sequence number within line", 4}},
    {5, {"Trace sequence number within SEG-Y file", 4}},
    {9, {"Original field record number", 4}},
    {13, {"Trace number within the original field record", 4}},
    {17, {"Energy source point number", 4}},
    {21, {"Ensemble number", 4}},
    {25, {"Trace number within the ensemble", 4}},
    {29, {"Trace identification code", 2}},
    {31, {"Number of vertically summed traces yielding this trace", 2}},
    {33, {"Number of horizontally stacked traces yielding this trace", 2}},
    {35, {"Data use for", 2}},
    {37,
     {"Distance from center of the source point to the center of the receiver "
      "group",
      4}},
    {41, {"Elevation of receiver group", 4}},
    {45, {"Surface elevation at source location", 4}},
    {49, {"Source depth below surface", 4}},
    {53, {"Seismic Datum elevation at receiver group", 4}},
    {57, {"Seismic Datum elevation at source", 4}},
    {61, {"Water column height at source location", 4}},
    {65, {"Water column height at receiver group location", 4}},
    {69, {"Scalar to be applied to all elevations and depths", 2}},
    {71, {"Scalar to be applied to all coordinates", 2}},
    {73, {"Source coordinate - X", 4}},
    {77, {"Source coordinate - Y", 4}},
    {81, {"Group coordinate - X", 4}},
    {85, {"Group coordinate - Y", 4}},
    {89, {"Coordinate units", 2}},
    {91, {"Weathering velocity", 2}},
    {93, {"Subweathering velocity", 2}},
    {95, {"Uphole time at source in ms", 2}},
    {97, {"Uphole time at group in ms", 2}},
    {99, {"Source static correction in ms", 2}},
    {101, {"Group static correction in ms", 2}},
    {103, {"Total static applied in ms", 2}},
    {105, {"Lag time A", 2}},
    {107, {"Lag Time B", 2}},
    {109, {"Delay recording time", 2}},
    {111, {"Mute time - Start time in ms", 2}},
    {113, {"Mute time - End time in ms", 2}},
    {115, {"Number of samples in this trace", 2}},
    {117, {"Sample interval for this trace", 2}},
    {119, {"Gain type of field instruments", 2}},
    {121, {"Instrument gain constant", 2}},
    {123, {"Instrument gain constant", 2}},
    {125, {"Correlated", 2}},
    {127, {"Sweep frequency at start", 2}},
    {129, {"Sweep frequency at end", 2}},
    {131, {"Sweep length in ms", 2}},
    {133, {"Sweep type", 2}},
    {135, {"Sweep trace taper length at start in ms", 2}},
    {137, {"Sweep trace taper length at end in ms", 2}},
    {139, {"Taper type", 2}},
    {141, {"Alias filter frequency (Hz)", 2}},
    {143, {"Alias filter slope (dB/octave)", 2}},
    {145, {"Notch filter frequency (Hz)", 2}},
    {147, {"Notch filter slope (dB/octave)", 2}},
    {149, {"Low-cut frequency (Hz)", 2}},
    {151, {"High-cut frequency (Hz)", 2}},
    {153, {"Low-cut slope (dB/octave)", 2}},
    {155, {"High-cut slope (dB/octave)", 2}},
    {157, {"Year data recorded", 2}},
    {159, {"Day of year", 2}},
    {161, {"Hour of day", 2}},
    {163, {"Minute of hour", 2}},
    {165, {"Second of minute", 2}},
    {167, {"Time basis code", 2}},
    {169, {"Trace weighting factor", 2}},
    {171, {"Geophone group number of roll switch position one", 2}},
    {173,
     {"Geophone group number of trace number one within original field record",
      2}},
    {175,
     {"Geophone group number of last trace within original field record", 2}},
    {177, {"Gap size (total number of groups dropped)", 2}},
    {179, {"Over travel associated with taper at beginning or end of line", 2}},
    {181, {"X coordinate", 4}},
    {185, {"Y coordinate", 4}},
    {189, {"The in-line number", 4}},
    {193, {"The cross-line number", 4}},
    {197, {"Shotpoint number", 4}},
    {201, {"Scalar to be applied to the shotpoint number", 2}},
    {203, {"Trace value measurement unit", 2}},
    {205, {"Transduction Constant", 6}},
    {211, {"Transduction Units", 2}},
    {213, {"Device/Trace Identifier", 2}},
    {215, {"Scalar to be applied to bytes 95-114", 2}},
    {217, {"Source Type/Orientation", 2}},
    {219,
     {"Source Energy Direction with respect to the source orientation", 6}},
    {225,
     {"Source Measurement - Describes the source effort used to generate "
      "the trace",
      6}},
    {231, {"Source Measurement Unit", 2}},
    {233, {"Either binary zeros or chars SEG00000", 8}}};

const std::map<int, const char *> kTraceSortingHelp = {
    {-1, "Other"},
    {0, "Unknown"},
    {1, "As recorded (no sorting)"},
    {2, "CDP ensemble"},
    {3, "Single fold continuous profile"},
    {4, "Horizontally stacked"},
    {5, "Common source point"},
    {6, "Common receiver point"},
    {7, "Common offset point"},
    {8, "Common mid-point"},
    {9, "Common conversion point"}};

const uint64_t kMaxLSeekSize = std::numeric_limits<long>::max();

// NOTE: only support SEGY-v1
const std::map<int, int> kElementSize = {{1, 4}, {2, 4}, {3, 2}, {4, 4}, {5, 4}, {8, 1}};

inline void swap_endian_inplace(void *dst, const void *src, int n) {
  uchar *_dst = static_cast<uchar *>(dst);
  const uchar *_src = static_cast<const uchar *>(src);
  static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");
  memcpy(_dst, _src, n);

  if (n <= 16) {
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

inline void read_binary_header(void *dst, const void *src) {
  uchar *_dst = static_cast<uchar *>(dst);
  const uchar *_src = static_cast<const uchar *>(src);
  for (auto &s : kBinaryHeaderHelp) {
    int loc = s.first - 1;
    int len = s.second.second;
    if (len % 2 != 0 && len > 1 && len <= 16) {
      len -= 1;
    }
    swap_endian_inplace(_dst + loc, _src + loc, len);
  }
}

inline void read_one_trace_header(void *dst, const void *src) {
  uchar *_dst = static_cast<uchar *>(dst);
  const uchar *_src = static_cast<const uchar *>(src);
  for (auto &s : kTraceHeaderHelp) {
    int loc = s.first - 1;
    int len = s.second.second;
    if (len % 2 != 0 && len > 1 && len <= 16) {
      len -= 1;
    }
    swap_endian_inplace(_dst + loc, _src + loc, len);
  }
}


template <typename T> void convert2npT(float *dst, const void *src, int size, int dformat) {
  const T *_src = static_cast<const T *>(src);
  for (int i = 0; i < size; ++i) {
    if (dformat == 1) {
      dst[i] = ibm_to_ieee(_src[i], true);
    } else {
      dst[i] = float(swap_endian<T>(_src[i]));
    }
  }
}


inline void convert2np(float *dst, const char *src, int size, int dformat) {
  if (dformat == 1) convert2npT<float>(dst, src, size, dformat);
  else if (dformat == 2) convert2npT<int32_t>(dst, src, size, dformat);
  else if (dformat == 3) convert2npT<int16_t>(dst, src, size, dformat);
  else if (dformat == 5) convert2npT<float>(dst, src, size, dformat);
  else if (dformat == 8) convert2npT<int8_t>(dst, src, size, dformat);
  else if (dformat == 10) convert2npT<uint32_t>(dst, src, size, dformat);
  else if (dformat == 11) convert2npT<uint16_t>(dst, src, size, dformat);
  else if (dformat == 16) convert2npT<uint8_t>(dst, src, size, dformat);
}

template <typename T> void float2sgyT(void *dst, const float *src, int size, int dformat) {
  T *_dst = static_cast<T *>(dst);

  for (int i = 0; i < size; ++i) {
    if (dformat == 1) {
      _dst[i] = ibm_to_ieee(src[i], true);
    } else {
      _dst[i] = swap_endian<T>(T(src[i]));
    }
  }
}

inline void float2sgy(char *dst, const float *src, int size, int dformat) {
  if (dformat == 1) float2sgyT<float>(dst, src, size, dformat);
  else if (dformat == 2) float2sgyT<int32_t>(dst, src, size, dformat);
  else if (dformat == 3) float2sgyT<int16_t>(dst, src, size, dformat);
  else if (dformat == 5) float2sgyT<float>(dst, src, size, dformat);
  else if (dformat == 8) float2sgyT<int8_t>(dst, src, size, dformat);
  else if (dformat == 10) float2sgyT<uint32_t>(dst, src, size, dformat);
  else if (dformat == 11) float2sgyT<uint16_t>(dst, src, size, dformat);
  else if (dformat == 16) float2sgyT<uint8_t>(dst, src, size, dformat);
}


} // namespace segy

#endif
