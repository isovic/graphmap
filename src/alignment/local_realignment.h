/*
 * local_realignment.h
 *
 *  Created on: Sep 3, 2014
 *      Author: Ivan Sovic
 */

#ifndef MYERS_WRAPPER_H_
#define MYERS_WRAPPER_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <sstream>
#include <vector>

#include "index/index_sa.h"
#include "sequences/single_sequence.h"
#include "utility/utility_general.h"
#include "utility/program_parameters.h"
#include "containers/path_graph_entry.h"
#include "alignment/myers.h"
#include "alignment/cigargen.h"
#include "log_system/log_system.h"
#include "containers/region.h"
#include "seqan/basic.h"
#include "seqan/align.h"
#include "seqan/sequence.h"
#include "seqan/stream.h"

#include "alignment/local_realignment_wrappers.h"



//void VerboseAlignment(SingleSequence *read, ProgramParameters &parameters, VertexData &vertex, std::string &cigar, int current_score,
//                      const unsigned char* query, const int query_length,
//                      const unsigned char* target, const int target_length,
//                      const unsigned char* alignment, const int alignment_length,
//                      const int position, const int mode_code);



//typedef int (*AlignmentFunctionType)(int8_t *read_data, int64_t read_length,
//                 int8_t *reference_data, int64_t reference_length,
//                 int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
//                 int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
//                 int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);
//
//int MyersEditDistanceWrapper(int8_t *read_data, int64_t read_length,
//                 int8_t *reference_data, int64_t reference_length,
//                 int64_t *ret_alignment_position_end,
//                 int64_t *ret_edit_distance, int myers_mode_code=MYERS_MODE_HW);

typedef int (*AlignmentFunctionType)(const int8_t*, int64_t, const int8_t*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t*, int64_t*, int64_t*, std::vector<unsigned char> &);
typedef int (*EditDistanceFunctionType)(const int8_t*, int64_t, const int8_t*, int64_t, int64_t*, int64_t*, int);

//int ClipCircularAlignment(int64_t alignment_start, int64_t alignment_end, unsigned char *alignment, int64_t alignment_length,
//                   int64_t read_length, int64_t reference_start, int64_t reference_length, int64_t split_region_start_offset, int64_t position_of_ref_end,
//                   unsigned char **ret_clipped_alignment, int64_t *ret_clipped_alignment_length,
//                   int64_t *ret_clipped_alignment_pos_start, int64_t *ret_clipped_alignment_pos_end);
int ClipCircularAlignment(int64_t alignment_start, int64_t alignment_end, unsigned char *alignment, int64_t alignment_length,
                          int64_t read_length, int64_t reference_start, int64_t reference_length, int64_t split_region_start_offset, int64_t position_of_ref_end,
                          unsigned char **ret_left_alignment, int64_t *ret_left_alignment_length,
                          int64_t *ret_left_alignment_pos_start, int64_t *ret_left_alignment_pos_end,
                          unsigned char **ret_right_alignment, int64_t *ret_right_alignment_length,
                          int64_t *ret_right_alignment_pos_start, int64_t *ret_right_alignment_pos_end);

//int LocalRealignmentLinear(AlignmentFunctionType AlignmentFunction, SingleSequence *read, Index *index, ProgramParameters &parameters, PathGraphEntry *best_path,
//                           int64_t unused_reference_length, int64_t *ret_alignment_position, std::string *ret_cigar,
//                           SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
//                           int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement=true);
//int LocalRealignmentCircular(AlignmentFunctionType AlignmentFunction, SingleSequence *read, Index *index, ProgramParameters &parameters, PathGraphEntry *best_path,
//                           int64_t unused_reference_length, int64_t *ret_alignment_position, std::string *ret_cigar,
//                           SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
//                           int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement=true);
//int HybridRealignment(SingleSequence *read, Index *index, ProgramParameters &parameters, PathGraphEntry *best_path, int64_t reference_length,  int64_t *ret_alignment_position, std::string *ret_cigar,
//                      SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
//                      int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement=true);

int LocalRealignmentLinear(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                           int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part,
                           int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part,
                           SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                           int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement=true);

int LocalRealignmentCircular(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                             int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part,
                             int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part,
                             SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                             int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement=true);
int HybridRealignment(const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                      int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part,
                      int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part,
                      SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                      int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement=true);

int CalcEditDistanceLinear(EditDistanceFunctionType EditDistanceFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                           int64_t *ret_alignment_position, int64_t *ret_edit_distance);
int CalcEditDistanceCircular(EditDistanceFunctionType EditDistanceFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                             int64_t *ret_alignment_position, int64_t *ret_edit_distance);



int CheckAlignmentSane(std::vector<unsigned char> &alignment, const SingleSequence* read=NULL, const Index* index=NULL, int64_t reference_hit_id=-1, int64_t reference_hit_pos=-1);
int CountAlignmentOperations(std::vector<unsigned char> &alignment, const int8_t *read_data, const int8_t *ref_data, int64_t reference_hit_id, int64_t alignment_position_start, SeqOrientation orientation,
                             int64_t match, int64_t mismatch, int64_t gap_open, int64_t gap_extend,
                             int64_t *ret_eq, int64_t *ret_x, int64_t *ret_i, int64_t *ret_d, int64_t *ret_alignment_score);

//int64_t RescoreAlignment(unsigned char *alignment, int64_t alignment_length, int64_t match, int64_t mismatch, int64_t gap_open, int64_t gap_extend);



//#ifndef RELEASE_VERSION
//  int LocalRealignmentLinearExperimental(AlignmentFunctionType AlignmentFunction, SingleSequence *read, Index *index, ProgramParameters &parameters, PathGraphEntry *best_path,
//                             int64_t unused_reference_length, int64_t *ret_alignment_position, std::string *ret_cigar,
//                             SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity);
//  int LocalRealignmentCircularExperimental(AlignmentFunctionType AlignmentFunction, SingleSequence *read, Index *index, ProgramParameters &parameters, PathGraphEntry *best_path,
//                             int64_t unused_reference_length, int64_t *ret_alignment_position, std::string *ret_cigar,
//                             SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity);
//#endif

#endif /* MYERS_WRAPPER_H_ */
