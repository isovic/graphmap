/*
 * cigargen.h
 *
 *  Created on: Aug 28, 2014
 *      Author: ivan
 */

#ifndef CIGARGEN_H_
#define CIGARGEN_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <algorithm>

#include "libs/myers.h"
#include "utility/utility_general.h"

#define EDLIB_M 0
#define EDLIB_EQUAL 0
#define EDLIB_X 3
#define EDLIB_I 1
#define EDLIB_D 2
#define EDLIB_S 4
#define EDLIB_H 5     /// Not used in GraphMap currently (26.01.2016.)
#define EDLIB_NOP 6

std::string AlignmentToCigar(unsigned char *alignment, int alignmentLength, bool extended_format);
int AlignmentToBasicCigar(unsigned char* alignment, int alignmentLength, char** cigar_);
int AlignmentToExtendedCigar(unsigned char* alignment, int alignmentLength, char** cigar_);
std::string AlignmentToMD(std::vector<unsigned char>& alignment, const int8_t *read_data, const int8_t *ref_data, int64_t reference_hit_id, int64_t alignment_position_start);

/// Searches for consecutive EDLIB_I and EDLIB_D (or vice versa) operations, and replaces the overlap with EDLIB_X.
std::vector<unsigned char> FixAlignment(unsigned char* alignment, int alignmentLength);
/// In case an alignment has leading/trailing EDLIB_I operations, they will be replaced with EDLIB_S.
int ConvertInsertionsToClipping(unsigned char* alignment, int alignmentLength);
/// Counts the number of leading and trailing clipped bases (or insertions).
int CountClippedBases(unsigned char* alignment, int alignmentLength, int64_t *ret_num_clipped_front, int64_t *ret_num_clipped_back);
/// Sums up the bases on the reference the alignment spans through (EDLIB_M and EDLIB_D operations).
int64_t CalculateReconstructedLength(unsigned char *alignment, int alignmentLength);
/// Counts each operation type, and calculates the alignment score as well (while rescoring the alignment with the given scores/penalties).
int CountAlignmentOperations(std::vector<unsigned char> &alignment, const int8_t *read_data, const int8_t *ref_data, int64_t reference_hit_id, int64_t alignment_position_start, SeqOrientation orientation,
                             int64_t match, int64_t mismatch, int64_t gap_open, int64_t gap_extend,
                             bool skip_leading_and_trailing_insertions,
                             int64_t *ret_eq, int64_t *ret_x, int64_t *ret_i, int64_t *ret_d, int64_t *ret_alignment_score, int64_t *ret_edit_dist, int64_t *ret_nonclipped_length);
/// Reverses the operations in a CIGAR string.
std::string ReverseCigarString(std::string &cigar);

std::string PrintAlignmentToString(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const int modeCode, int row_width=100);

int GetAlignmentPatterns(const unsigned char* query, const int64_t queryLength,
                         const unsigned char* target, const int64_t targetLength,
                         const unsigned char* alignment, const int64_t alignmentLength,
                         std::string &ret_query, std::string &ret_target, std::string &ret_match_pattern);

#endif /* CIGARGEN_H_ */
