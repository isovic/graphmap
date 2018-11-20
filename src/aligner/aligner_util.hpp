#ifndef SRC_ALIGNER_ALIGNER_UTIL_H_
#define SRC_ALIGNER_ALIGNER_UTIL_H_

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include "sam_parser.h"

namespace is {

template<typename T>
std::vector<T> GenerateSimpleMatchMatrix(T match, T mismatch, size_t alphabet_size) {
  std::vector<T> matrix(alphabet_size * alphabet_size, mismatch);  // Set the mismatch score.
  // Goes to "-1" to allow for 'N' bases which should not match to themselves.
  for (size_t i=0; i<(alphabet_size - 1); i++) {
    matrix[i*alphabet_size + i] = match;                  // Set the match score.
    matrix[i*alphabet_size + alphabet_size - 1] = 0;      // Reset the last column to 0.
    matrix[(alphabet_size - 1) * alphabet_size + i] = 0;  // Reset the last row to 0.
  }
  return matrix;
}

std::vector<int8_t> ConvertSeqAlphabet(const int8_t* seq, size_t seqlen, const uint8_t* conv_table);

std::vector<is::CigarOp> ConvertBasicToExtCIGAR(const char* qseq, int64_t qlen,
                                                const char* tseq, int64_t tlen,
                                                const std::vector<is::CigarOp>& basic_cigar);

int64_t EditDistFromExtCIGAR(const std::vector<is::CigarOp>& extended_cigar);

std::vector<is::CigarOp> ExtractCigarBetweenQueryCoords(const std::vector<is::CigarOp>& cigar, int64_t qstart, int64_t qend, int64_t *cigar_length, int64_t *cigar_length_q);

std::string CigarToString(const std::vector<is::CigarOp>& cigar);

}

#endif
