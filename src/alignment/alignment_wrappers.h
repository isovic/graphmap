/*
 * local_realignment_generic.h
 *
 *  Created on: Jan 16, 2015
 *      Author: isovic
 */

#ifndef LOCAL_REALIGNMENT_GENERIC_H_
#define LOCAL_REALIGNMENT_GENERIC_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <sstream>
#include <vector>

#include "utility/utility_general.h"
#include "program_parameters.h"
#include "libs/edlib.h"
#include "alignment/cigargen.h"
#include "log_system/log_system.h"
#include "containers/region.h"
#include "seqan/basic.h"
#include "seqan/align.h"
#include "seqan/sequence.h"
#include "seqan/stream.h"

#define ALIGNMENT_TYPE_SHW  0     /// Gaps at the end are not penalized.
#define ALIGNMENT_TYPE_HW   1     /// Gaps at the beginning and the end are not penalized.
#define ALIGNMENT_TYPE_NW   2     /// Global alignment (gaps at the beginning and the end are penalized).

#ifndef RELEASE_VERSION
  #include "libs/opal.h"
#endif

typedef int (*AlignmentFunctionType)(const int8_t*, int64_t, const int8_t*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t*, int64_t*, int64_t*, std::vector<unsigned char> &);
typedef int (*EditDistanceFunctionType)(const int8_t*, int64_t, const int8_t*, int64_t, int64_t*, int64_t*, int);

int LocalizeAlignmentPosWithMyers(const int8_t *read_data, int64_t read_length,
                                  const int8_t *reference_data, int64_t reference_length,
                                  int64_t rough_reference_start, int64_t rough_reference_end,
                                  int64_t *ret_alignment_start, int64_t *ret_alignment_end,
                                  int64_t *ret_start_ambiguity, int64_t *ret_end_ambiguity,
                                  int64_t *ret_edit_distance, int64_t *ret_band_width,
                                  bool verbose_debug_output=false);

int MyersSemiglobalWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);

int MyersNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);

int MyersSHWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);

#ifndef RELEASE_VERSION
int OpalNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);

int OpalSHWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);
#endif

int SeqAnSemiglobalWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);
int SeqAnSemiglobalWrapperWithMyersLocalization(const int8_t *read_data, int64_t read_length,
                                                const int8_t *reference_data, int64_t reference_length,
                                                int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                                                int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                                                int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);

int MyersEditDistanceWrapper(const int8_t *read_data, int64_t read_length,
                             const int8_t *reference_data, int64_t reference_length,
                             int64_t *ret_alignment_position_end,
                             int64_t *ret_edit_distance, EdlibAlignMode edlib_mode_code=EDLIB_MODE_HW);

int SeqAnNWWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);
int SeqAnSHWWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);

//int SeqAnAlignmentToEdlibAlignmentNoCigar(seqan::Align<seqan::Dna5String> &align, bool is_global, int64_t *ret_start_offset, int64_t *ret_end_offset, int64_t *edit_distance, std::vector<unsigned char> &ret_alignment);
int SeqAnAlignmentToEdlibAlignmentNoCigar(seqan::Align<seqan::Dna5String> &align, int alignment_type, int64_t *ret_start_offset, int64_t *ret_end_offset, int64_t *edit_distance, std::vector<unsigned char> &ret_alignment);

int CheckAlignmentSaneSimple(std::vector<unsigned char> &alignment);



#endif /* LOCAL_REALIGNMENT_GENERIC_H_ */
