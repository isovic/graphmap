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

std::string AlignmentToCigar(unsigned char *alignment, int alignmentLength, bool extended_format);
int AlignmentToBasicCigar(unsigned char* alignment, int alignmentLength, char** cigar_);
int AlignmentToExtendedCigar(unsigned char* alignment, int alignmentLength, char** cigar_);

int ConvertInsertionsToClipping(unsigned char* alignment, int alignmentLength);
int64_t CalculateReconstructedLength(unsigned char *alignment, int alignmentLength);
int CountAlignmentOperations(std::vector<unsigned char> &alignment, const int8_t *read_data, const int8_t *ref_data, int64_t reference_hit_id, int64_t alignment_position_start, SeqOrientation orientation,
                             int64_t match, int64_t mismatch, int64_t gap_open, int64_t gap_extend,
                             int64_t *ret_eq, int64_t *ret_x, int64_t *ret_i, int64_t *ret_d, int64_t *ret_alignment_score, int64_t *ret_nonclipped_length);

std::string ReverseCigarString(std::string &cigar);

std::string PrintAlignmentToString(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const int modeCode, int row_width=100);

#endif /* CIGARGEN_H_ */
