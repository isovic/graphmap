/*
 * alignment.cc
 *
 *  Created on: Jan 17, 2016
 *      Author: isovic
 */

#include "alignment/alignment.h"

//int Align(const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *region_data, bool extend_to_end,
//                      int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part, int64_t *ret_nonclipped_left_part,
//                      int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part, int64_t *ret_nonclipped_right_part,
//                      SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
//                      int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op) {
//
//
//    if (parameters.alignment_algorithm == "gotoh") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh for alignment!\n", "Alignment");
//      if (region_data->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//      else
//        return LocalRealignmentCircular(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//
//    } else if (parameters.alignment_algorithm == "myers") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
//      if (region_data->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//      else
//        return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//
//    } else if (parameters.alignment_algorithm == "anchor") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//      bool is_linear = region_data->get_region_data().is_split == false || parameters.is_reference_circular == false;
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignment(is_linear, extend_to_end, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//
//    } else if (parameters.alignment_algorithm == "anchorgotoh") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//      bool is_linear = region_data->get_region_data().is_split == false || parameters.is_reference_circular == false;
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh's algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignment(is_linear, extend_to_end, SeqAnNWWrapper, SeqAnSHWWrapper, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//
//#ifndef RELEASE_VERSION
//    } else if (parameters.alignment_algorithm == "anchormex") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//      /// Extension to read ends is not currently supported. SHW alignment needs to be implemented in Opal first.
//      extend_to_end = false;
//      bool is_linear = region_data->get_region_data().is_split == false || parameters.is_reference_circular == false;
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Match Extend algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignmentMex(is_linear, extend_to_end, OpalNWWrapper, OpalSHWWrapper, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//#endif
//
//    } else {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Warning: Unknown alignment algorithm selected. Using Myers' bit-vector alignment instead.\n", "Alignment");
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
//      if (region_data->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//      else
//        return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, region_data, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//    }
//
//  return -1;
//}
