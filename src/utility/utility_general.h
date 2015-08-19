/*
 * utility_general.h
 *
 *  Created on: Oct 11, 2014
 *      Author: ivan
 */

#ifndef UTILITY_GENERAL_H_
#define UTILITY_GENERAL_H_

#include <time.h>
#include <string>
#include <stdarg.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric> 


enum SeqOrientation {
  kForward = 0,
  kReverse = 1,
  kUnknown = 2
};

enum ProcessedReadState {
  kAligned = 0,
  kUnaligned = 1,
  kAmbiguous = 2
};

//#define USE_EXTENDED_CIGAR_FORMAT

#define STATE_UNMAPPED  0
#define STATE_MAPPED    1
#define STATE_AMBIGUOUS 2

#define SAM_DEFAULT_QNAME "*"
#define SAM_DEFAULT_FLAG  0
#define SAM_DEFAULT_RNAME "*"
#define SAM_DEFAULT_POS   0
#define SAM_DEFAULT_MAPQ  255
#define SAM_DEFAULT_CIGAR "*"
#define SAM_DEFAULT_RNEXT "*"
#define SAM_DEFAULT_PNEXT 0
#define SAM_DEFAULT_TLEN  0
#define SAM_DEFAULT_SEQ "*"
#define SAM_DEFAULT_QUAL  "*"
#define SAM_DEFAULT_AS  0
#define SAM_DEFAULT_METAGEN_AS  0
#define SAM_DEFAULT_NUM_PERFECT_HITS  1
#define SAM_DEFAULT_NM  -1

#define SAM_MULTIPLE_SEGMENTS             1 << 0
#define SAM_EACH_SEG_PROPERLY_ALIGNED     1 << 1
#define SAM_THIS_SEG_UNMAPPED               1 << 2
#define SAM_NEXT_SEG_UNMAPPED               1 << 3
#define SAM_THIS_SEG_REVERSED             1 << 4
#define SAM_NEXT_SEG_REVERSED             1 << 5
#define SAM_THIS_SEG_FIRST                1 << 6
#define SAM_THIS_SEG_LAST                 1 << 7
#define SAM_SECONDARY_ALIGNMENT           1 << 8
#define SAM_PASSING_QUAL_SCORES           1 << 9
#define SAM_PCR_OR_OPTICAL_DUP            1 << 10
#define SAM_SUPPLEMENTARY_ALIGN           1 << 11

#define ALIGNMENT_GOOD                0
#define ALIGNMENT_NOT_SANE            -123
#define ALIGNMENT_WRONG_CLUSTER_SIZE  -12
#define ALIGNMENT_NOT_CIRCULAR        -15
#define ALIGNMENT_WRONG_DATA          -1
#define ALIGNMENT_CONVERSION_PROBLEM  -2
#define ALIGNMENT_LOCALIZATION_PROBLEM  -3
#define ALIGNMENT_MYERS_INTERNAL_ERROR  -4
#define ALIGNMENT_DISTANCE_BETWEEN_ANCHORS_PROBLEM  -5



int fileExists(const char *fname);
std::string GetUTCTime(std::string fmt="%a, %d %b %y %T %z");
std::string GetLocalTime();
void ProcessMemUsage(double& vm_usage, double& resident_set);
std::string FormatMemoryConsumptionAsString();
size_t getCurrentRSS( );
size_t getPeakRSS( );

// Takes a C-style formatted input string (same as printf), and converts it
// to a C++ std::string.
// Sample usage:
//  std::string formatted_string = FormatString ("This is just a %s! %d + %d = %d\n", "test", 1, 2, 3);
std::string FormatString(const char* additional_message, ...);
std::string FormatStringToLength(std::string original_string, uint32_t length);
void PrintSubstring(char *text, int64_t length, FILE *fp=stdout);
std::string GetSubstring(char *text, int64_t length);
unsigned char* CreateReverseCopy(const unsigned char* seq, uint64_t length);
std::string TrimToFirstSpace(std::string original_string);

// Calculates the sigmoid function with given parameters.
// S(t) = 1.0f / (1.0f + exp(-t));
// t is real, <-inf, inf>, S(t) is real, <0, 1>
// For t = +-6, S(t) is almost equal to 1 (or zero, respectivelly).
// Parameter:
//   @mean  shift on the X-axis
//   @width for (x + width) the value of sigmoid will be equal to S(t + 6)
float sigmoid(float x, float mean, float width);

template <typename T>
std::vector<size_t> ordered_sort_vector(std::vector<T> const& values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

template <typename T>
std::vector<size_t> ordered_sort_array(const T* values, int64_t values_size) {
    std::vector<size_t> indices(values_size);
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

/// User needs to free the memory using 'free'.
template <typename T>
T* reverse_data(const T* data, int64_t data_len) {
  T* ret_data = (T *) malloc(sizeof(T) * data_len);
  for (int64_t i=0; i<data_len; i++)
    ret_data[i] = data[data_len - i - 1];
  return ret_data;
}

int GetClippingOpsFromCigar(const std::string &cigar, char *clip_op_front, int64_t *clip_count_front, char *clip_op_back, int64_t *clip_count_back);



#endif /* UTILITY_GENERAL_H_ */
