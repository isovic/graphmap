/*
 * myers_wrapper.cc
 *
 *  Created on: Sep 3, 2014
 *      Author: ivan
 */

#include "alignment/local_realignment.h"
#include "index/index_spaced_hash_fast.h"



int HybridRealignment(const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                      int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part, int64_t *ret_nonclipped_left_part,
                      int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part, int64_t *ret_nonclipped_right_part,
                      SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                      int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement) {


    if (parameters.alignment_algorithm == "gotoh") {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh for alignment!\n", "Alignment");
      if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
        return LocalRealignmentLinear(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
      else
        return LocalRealignmentCircular(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);

    } else if (parameters.alignment_algorithm == "myers") {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
      if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
        return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
      else
        return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);

    } else if (parameters.alignment_algorithm == "anchor") {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
      bool extend_to_end = true;
//      bool extend_to_end = false;
      bool is_linear = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
      return AnchoredAlignment(is_linear, extend_to_end, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//      return AnchoredAlignment(is_linear, extend_to_end, MyersNWWrapper, SeqAnSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//      return AnchoredAlignment(is_linear, extend_to_end, SeqAnNWWrapper, SeqAnSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);

    } else if (parameters.alignment_algorithm == "anchorgotoh") {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
      bool extend_to_end = true;
//      bool extend_to_end = false;
      bool is_linear = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh's algorithm for alignment!\n", "Alignment");
      return AnchoredAlignment(is_linear, extend_to_end, SeqAnNWWrapper, SeqAnSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//      return AnchoredAlignment(is_linear, extend_to_end, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);

    } else if (parameters.alignment_algorithm == "anchormex") {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//      bool extend_to_end = true;
      bool extend_to_end = false;
      bool is_linear = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Match Extend algorithm for alignment!\n", "Alignment");
      return AnchoredAlignmentMex(is_linear, extend_to_end, OpalNWWrapper, OpalSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);

    } else {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Warning: Unknown alignment algorithm selected. Using Myers' bit-vector alignment instead.\n", "Alignment");
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
      if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
        return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
      else
        return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
    }



//  if (parameters.alignment_approach == "sg") {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
//
//    if (parameters.alignment_algorithm == "gotoh") {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh for alignment!\n", "Alignment");
//      if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//      else
//        return LocalRealignmentCircular(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//
//    } else if (parameters.alignment_algorithm == "myers") {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
//      if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//      else
//        return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//
//    } else {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Warning: Unknown alignment algorithm selected. Using Myers' bit-vector alignment instead.\n", "Alignment");
//      if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//      else
//        return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//    }
//
//  } else if (parameters.alignment_approach == "anchor" || parameters.alignment_approach == "overlap") {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//    /// TODO: Make this an update in the next release. Enable both non-overlap anchored mode and overlap. Also, enable Gotoh for anchored alignment.
////    bool extend_to_end = (parameters.alignment_approach == "overlap");
//    bool extend_to_end = true;
//    bool is_linear = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
//
////    if (parameters.alignment_algorithm == "gotoh") {
////      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh for alignment!\n", "Alignment");
////      return AnchoredAlignment(is_linear, extend_to_end, SeqAnNWWrapper, SeqAnSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
////
////    } else if (parameters.alignment_algorithm == "myers") {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignment(is_linear, extend_to_end, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//
////    } else {
////      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Warning: Unknown alignment algorithm selected. Using Myers' bit-vector alignment instead.\n", "Alignment");
////      return AnchoredAlignment(is_linear, extend_to_end, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
////    }
//
//  }







//  } else if (parameters.alignment_algorithm == "anchor") {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment for local realignment!\n\n", "HybridRealignment");
//    bool is_linear = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
//    return AnchoredAlignment(is_linear, true, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//
////    if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
////      return AnchoredAlignmentLinear(MyersNWWrapper, MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
////    else
////      return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
////      return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
////      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Circular anchored alignment not implemented yet!!\n\n", "HybridRealignment");
//
//  } else if (parameters.alignment_algorithm == "anchorsq") {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment with SeqAn for local realignment!\n\n", "HybridRealignment");
//    bool is_linear = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
//    return AnchoredAlignment(is_linear, true, SeqAnNWWrapper, SeqAnSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//
//  } else if (parameters.alignment_algorithm == "anchornoext") {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment for local realignment!\n\n", "HybridRealignment");
//    bool is_circular = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
//    return AnchoredAlignment(is_circular, false, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//
//  } else if (parameters.alignment_algorithm == "anchorsqnoext") {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment with SeqAn for local realignment!\n\n", "HybridRealignment");
//    bool is_linear = best_path->get_region_data().is_split == false || parameters.is_reference_circular == false;
//    return AnchoredAlignment(is_linear, false, SeqAnNWWrapper, SeqAnSHWWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//
//  } else {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using EDlib for local realignment!\n\n", "HybridRealignment");
//    if (best_path->get_region_data().is_split == false || parameters.is_reference_circular == false)
//      return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//    else
//      return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, best_path, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op, perform_reverse_complement);
//  }

  return -1;
}

//void VerboseAlignment(SingleSequence *read, ProgramParameters &parameters, VertexData &vertex, std::string &cigar, int current_score,
//                      const unsigned char* query, const int query_length,
//                      const unsigned char* target, const int target_length,
//                      const unsigned char* alignment, const int alignment_length,
//                      const int position, const int mode_code) {
//  std::string alignment_as_string = "";
//  alignment_as_string = PrintAlignmentToString(query, query_length,
//                                               target, target_length,
//                                               alignment, alignment_length,
//                                               position, mode_code);
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
//                                           FormatString("%s\nRef:   %s\nRead:  %s\nCIGAR: %s\nScore: %d\nAlignment:\n%s\n\n",
//                                                        vertex.VerboseToString().c_str(),
//                                                        GetSubstring((char *) query, query_length).c_str(),
//                                                        GetSubstring((char *) target, target_length).c_str(),
//                                                        cigar.c_str(), current_score, alignment_as_string.c_str()), "[]");
//}



// Example usage:
// ClipCircularAlignment(best_aligning_position_start, best_aligning_position_end, alignment, alignment_length,
// read->get_sequence_length(), index->get_reference_starting_pos()[absolute_reference_id], index->get_reference_lengths()[absolute_reference_id], start_offset, position_of_ref_end,
// &new_alignment, &new_alignment_length, &new_alignment_start, &new_alignment_end);
// if (alignment)
//  free(alignment);
// alignment = new_alignment;
// alignment_length = new_alignment_length;
// best_aligning_position_start = new_alignment_start;
// best_aligning_position_end = new_alignment_end;
//
// User should free the memory allocated under ret_left_alignment and ret_right_alignment parameters. Use free() for freeing the memory.
int ClipCircularAlignment(int64_t alignment_start, int64_t alignment_end, unsigned char *alignment, int64_t alignment_length,
                          int64_t read_length, int64_t reference_start, int64_t reference_length, int64_t split_region_start_offset, int64_t position_of_ref_end,
                          unsigned char **ret_left_alignment, int64_t *ret_left_alignment_length,
                          int64_t *ret_left_alignment_pos_start, int64_t *ret_left_alignment_pos_end,
                          unsigned char **ret_right_alignment, int64_t *ret_right_alignment_length,
                          int64_t *ret_right_alignment_pos_start, int64_t *ret_right_alignment_pos_end) {
  // Check for clipping and report the final result.
  // In this case, we need to clip a part of the alignment.
  if (position_of_ref_end >= alignment_start && position_of_ref_end < alignment_end) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Alignment needs to be clipped."), "[]");

    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Alignment needs to be clipped."), "[]");

    // Find the base on the read which is on the position of the 0th base of the reference
    int64_t pos_on_read = 0, pos_on_ref = alignment_start, pos_on_alignment = 0;
    int64_t last_m_pos_on_read = -1, last_m_pos_on_ref = -1, last_m_pos_on_alignment = -1;
    int64_t next_m_pos_on_read = -1, next_m_pos_on_ref = -1, next_m_pos_on_alignment = -1;
    for (int64_t i=0; i<alignment_length; i++) {
      if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X) {
        if (pos_on_ref <= position_of_ref_end) {
          last_m_pos_on_read = pos_on_read;
          last_m_pos_on_ref = pos_on_ref;
          last_m_pos_on_alignment = pos_on_alignment;
        }
        if (pos_on_ref >= (position_of_ref_end + 1) && next_m_pos_on_read == -1) {
          next_m_pos_on_read = pos_on_read;
          next_m_pos_on_ref = pos_on_ref;
          next_m_pos_on_alignment = pos_on_alignment;
          break;
        }
      }

      if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X || alignment[i] == EDLIB_S || alignment[i] == EDLIB_I)
        pos_on_read += 1;
      if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X || alignment[i] == EDLIB_D)
        pos_on_ref += 1;
      pos_on_alignment += 1;
    }

  //  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("read_length = %ld\n", read->get_sequence_length()), "MyersLocalRealignment3");
  //  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("pos_on_read = %ld, pos_on_ref = %ld, pos_on_alignment = %ld\n", pos_on_read, pos_on_ref, pos_on_alignment), "MyersLocalRealignment3");
  //  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("last_m_pos_on_read = %ld, last_m_pos_on_ref = %ld, last_m_pos_on_alignment = %ld\n", last_m_pos_on_read, last_m_pos_on_ref, last_m_pos_on_alignment), "MyersLocalRealignment3");
  //  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("next_m_pos_on_read = %ld, next_m_pos_on_ref = %ld, next_m_pos_on_alignment = %ld\n\n", next_m_pos_on_read, next_m_pos_on_ref, next_m_pos_on_alignment), "MyersLocalRealignment3");

    // For the left part:
    // last_m_pos_on_read is the last base on read that is located left of the reference end.
    // That doesn't mean that it is exactly on the last base of the reference, as there could have been a deletion.
    // We need to check if last_m_pos_on_ref is < position_of_ref_end . If yes, the rest needs to be filled with S operations.
    //
    // split_region_start_offset - the absolute offset of the starting position of the split region, used in the data_ variable of the reference index. In other words, it's the
    //                             position of the region start before the reference end, in absolute coordinates. So it has the information about the reference start already.
    *ret_left_alignment_pos_start = ((alignment_start + split_region_start_offset - reference_start) % (reference_length)) + reference_start;
    /// The last part compensates if there was actually an indel in between the last M operation and the required clip length.
    *ret_left_alignment_pos_end = reference_start + reference_length - 1 - (position_of_ref_end - last_m_pos_on_ref);
//    printf ("last_m_pos_on_ref = %ld\n", last_m_pos_on_ref);
//    printf ("position_of_ref_end = %ld\n", position_of_ref_end);
//    fflush(stdout);

//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Left part:\n"), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("ret_left_alignment_pos_start = %ld\n", (*ret_left_alignment_pos_start)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("ret_left_alignment_pos_end = %ld\n", (*ret_left_alignment_pos_end)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("alignment_start = %ld\n", (alignment_start)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("split_region_start_offset = %ld\n", (split_region_start_offset)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("reference_start = %ld\n", (reference_start)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("reference_length = %ld\n", (reference_length)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, "", "[]");

    int64_t num_clipped_bases_left = read_length - last_m_pos_on_read - 1;
    int64_t new_alignment_length_left = last_m_pos_on_alignment + 1 + num_clipped_bases_left;
    unsigned char *new_alignment_left = (unsigned char *) malloc(sizeof(char) * (new_alignment_length_left + 1));
    // Copy the useful operations into a new array.
    memmove(new_alignment_left, alignment, (last_m_pos_on_alignment + 1));
    // Add soft clipping after useful operations.
    memset(&new_alignment_left[last_m_pos_on_alignment + 1], EDLIB_S, num_clipped_bases_left);
    new_alignment_left[new_alignment_length_left] = '\0';
    *ret_left_alignment = new_alignment_left;
    *ret_left_alignment_length = new_alignment_length_left;

    // For the right part:
    // If next_m_pos_on_ref > (position_of_ref_end + 1), then the bases in between need to be filled with S operations
    // If there were I or D bases between the start of the reference and the first base of the read that was mapped,
    // compensate the starting position of the alignment.
    //
    // There is a bug in this code. I think it's because of the "- reference_start" part before the modulo operation. Need to further test it.
//    *ret_right_alignment_pos_start = (((next_m_pos_on_ref - (position_of_ref_end + 1)) - reference_start) % (reference_length)) + reference_start;
    *ret_right_alignment_pos_start = (((next_m_pos_on_ref - (position_of_ref_end + 1))) % (reference_length)) + reference_start;
    *ret_right_alignment_pos_end = ((alignment_end + split_region_start_offset - reference_start) % (reference_length)) + reference_start;
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Right part:\n"), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("ret_right_alignment_pos_start = %ld\n", (*ret_right_alignment_pos_start)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("ret_right_alignment_pos_end = %ld\n", (*ret_right_alignment_pos_end)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("next_m_pos_on_ref = %ld\n", (next_m_pos_on_ref)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("position_of_ref_end = %ld\n", (position_of_ref_end)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("alignment_start = %ld\n", (alignment_start)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("alignment_end = %ld\n", (alignment_end)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("reference_start = %ld\n", (reference_start)), "ClipCircularAlignment");
//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("reference_length = %ld\n", (reference_length)), "ClipCircularAlignment");
//Problem je izgleda ovdje! Iz SAM filea izgleda da je stvar u right dijelu, ali mi koordinata nema smisla.
    int64_t num_clipped_bases_right = next_m_pos_on_read;
    int64_t new_alignment_length_right = (alignment_length - next_m_pos_on_alignment) + num_clipped_bases_right;
    unsigned char *new_alignment_right = (unsigned char *) malloc(sizeof(char) * (new_alignment_length_right + 1));

//    printf ("\n");
//    printf ("new_alignment_length_right = %ld\n", new_alignment_length_right);
//    printf ("num_clipped_bases_right = %ld\n", num_clipped_bases_right);
//    printf ("position_of_ref_end = %ld\n", position_of_ref_end);
//    printf ("alignment_length = %ld\n", alignment_length);
//    printf ("next_m_pos_on_alignment = %ld\n", next_m_pos_on_alignment);
//    printf ("position_of_ref_end = %ld\n", position_of_ref_end);
//    printf ("alignment_start = %ld\n", alignment_start);
//    printf ("alignment_end = %ld\n", alignment_end);
//    printf ("\n");
//    fflush(stdout);

    memset(&new_alignment_right[0], EDLIB_S, num_clipped_bases_right);
    memmove(&new_alignment_right[num_clipped_bases_right], &alignment[next_m_pos_on_alignment], (alignment_length - next_m_pos_on_alignment));
    new_alignment_right[new_alignment_length_right] = '\0';
    *ret_right_alignment = new_alignment_right;
    *ret_right_alignment_length = new_alignment_length_right;

    return 1;

  // In this case, the entire alignment is located at one end of the reference. No clipping is necessary.
  } else {
    // Only positions are returned, as the alignment does not need to be modified.
    *ret_left_alignment_pos_start = ((alignment_start + split_region_start_offset - reference_start) % (reference_length)) + reference_start;
    *ret_left_alignment_pos_end = ((alignment_end + split_region_start_offset - reference_start) % (reference_length)) + reference_start;
    *ret_left_alignment = NULL;
    *ret_left_alignment_length = 0;
    *ret_right_alignment = NULL;
    *ret_right_alignment_length = 0;

    return 0;
  }

  return -1;
}






//int64_t RescoreAlignment(unsigned char *alignment, int64_t alignment_length, int64_t match, int64_t mismatch, int64_t gap_open, int64_t gap_extend) {
//  int64_t alignment_score = 0;
//
//  for (int64_t i=0; i<alignment_length; i++) {
//    if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL) {
//      alignment_score += match;
//    }
//  }
//
//  return alignment_score;
//}

int LocalRealignmentLinear(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                           int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part, int64_t *ret_nonclipped_left_part,
                           int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part, int64_t *ret_nonclipped_right_part,
                           SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                           int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement) {

//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, "More generic implementation of the alignment step.\n", "LocalRealignmentLinear");

  int64_t absolute_reference_id = best_path->get_region_data().reference_id;
  int64_t reference_id = best_path->get_region_data().reference_id;
  int64_t reference_start = index->get_reference_starting_pos()[absolute_reference_id];
  int64_t reference_length = index->get_reference_lengths()[absolute_reference_id];
  SeqOrientation orientation = (best_path->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t best_aligning_position = 0;

  int64_t l1_reference_start = best_path->get_l1_data().l1_lmin;
  int64_t l1_reference_end = ((int64_t) (best_path->get_l1_data().l1_k * read->get_sequence_length())) + best_path->get_l1_data().l1_lmax;
  if (l1_reference_start < reference_start)
    l1_reference_start = reference_start;
  if (l1_reference_end >= (reference_start + reference_length))
    l1_reference_end = reference_start + reference_length - 1;
  int64_t reference_data_length = l1_reference_end - l1_reference_start + 1;

#ifndef RELEASE_VERSION
  if (parameters.verbose_level > 5 && read->get_sequence_id() == parameters.debug_read) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                             FormatString("\nl1_reference_start = %ld\nl1_reference_end = %ld\nreference_data_length = %ld\nreference_start = %ld\nreference_length = %ld\nabsolute_reference_id = %ld\norientation = %ld\n\n",
                                                          l1_reference_start, l1_reference_end, reference_data_length, reference_start, reference_length, absolute_reference_id, ((orientation == kForward) ? "forward" : "reverse")), "LocalRealignmentLinear");
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                             FormatString("\n%s\n", best_path->VerboseInfoToString().c_str()), "LocalRealignmentLinear");

  }
#endif

  int64_t alignment_position_start = 0, alignment_position_end = 0, edit_distance = 0;
  std::vector<unsigned char> alignment;
  int ret_code = AlignmentFunction(read->get_data(), read->get_sequence_length(),
                                   (int8_t *) (index->get_data() + l1_reference_start), reference_data_length,
                                   -1, parameters.match_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                   &alignment_position_start, &alignment_position_end,
                                   &edit_distance, alignment);
  alignment_position_start += l1_reference_start;
  alignment_position_end += l1_reference_start;

  if (ret_code != 0 || alignment.size() == 0)
    return ret_code;

  ConvertInsertionsToClipping((unsigned char *) &(alignment[0]), alignment.size());
  *ret_cigar_left_part = AlignmentToCigar((unsigned char *) &(alignment[0]), alignment.size());
//  *ret_AS_left_part = RescoreAlignment((unsigned char *) &(alignment[0]), alignment.size(), parameters.match_score, parameters.mismatch_penalty, parameters.gap_open_penalty, parameters.gap_extend_penalty);
  *ret_cigar_right_part = "";


#ifndef RELEASE_VERSION
  if (parameters.verbose_level > 5 && read->get_sequence_id() == parameters.debug_read) {
    std::string alignment_as_string = "";
    alignment_as_string = PrintAlignmentToString((const unsigned char *) (read->get_data()), read->get_sequence_length(),
                                               (const unsigned char *) (index->get_data() + alignment_position_start), (alignment_position_end - alignment_position_start + 1),
                                               (unsigned char *) &(alignment[0]), alignment.size(),
                                               (0), MYERS_MODE_NW);
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                             FormatString("Alignment:\n%s\n\nalignment_position_start = %ld\n\n", alignment_as_string.c_str(), alignment_position_start), "LocalRealignmentLinear");
  }
#endif



  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
    std::string alignment_as_string = "";
    alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data(), read->get_sequence_length(),
                                                 (const unsigned char *) (index->get_data() + l1_reference_start), reference_data_length,
                                                 (unsigned char *) &(alignment[0]), alignment.size(),
                                                 (alignment_position_end - l1_reference_start), MYERS_MODE_HW);

    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                             FormatString("alignment_position_start = %ld\nalignment_position_end = %ld\nl1_reference_start = %ld\n",
                                                          alignment_position_start, alignment_position_end, l1_reference_start), "LocalRealignmentLinear");

//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
//                                             FormatString("Ref:   %s\nRead:  %s\nCIGAR: %s\nEdit distance: %ld\nAlignment:\n%s\n\n",
//                                                          GetSubstring((char *) (index->get_data() + alignment_position_start), (alignment_position_end - alignment_position_start + 1)).c_str(),
//                                                          GetSubstring((char *) read->get_data(), read->get_sequence_length()).c_str(),
//                                                          (*ret_cigar).c_str(), edit_distance, alignment_as_string.c_str()), "[]");
  }

  if (orientation == kForward) {
    (index)->RawPositionConverterWithRefId(alignment_position_start, reference_id, 0, NULL, &best_aligning_position, NULL);
//    printf ("Tu sam 1!\n");
//    printf ("best_aligning_position = %ld\n", best_aligning_position);
//    fflush(stdout);
  } else {
//    printf ("Tu sam 2!\n");
//    fflush(stdout);
    index->RawPositionConverterWithRefId(alignment_position_end, reference_id, 0, NULL, &best_aligning_position, NULL);
    if (perform_reverse_complement == true) {
//      *ret_cigar_left_part = ReverseCigarString((*ret_cigar_left_part));
//      read->ReverseComplement();
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                               FormatString("ERROR! Tried to explicitly reverse complement a read! This functionality is no longer allowed!"), "LocalRealignmentLinear");
      exit(1);
    }
    reference_id -= index->get_num_sequences_forward();
  }

  *ret_alignment_position_left_part = best_aligning_position;
  *ret_alignment_position_right_part = 0;
  *ret_orientation = orientation;
  *ret_reference_id = reference_id;
  *ret_position_ambiguity = 0;

  if (CheckAlignmentSane(alignment, read, index, reference_id, best_aligning_position))
    return -1;

  CountAlignmentOperations(alignment, read->get_data(), index->get_data(), reference_id, alignment_position_start, orientation,
                           parameters.evalue_match, parameters.evalue_mismatch, parameters.evalue_gap_open, parameters.evalue_gap_extend,
                           ret_eq_op, ret_x_op, ret_i_op, ret_d_op, ret_AS_left_part, ret_nonclipped_left_part);

  alignment.clear();

  return ((int) edit_distance);
}

int CalcEditDistanceLinear(EditDistanceFunctionType EditDistanceFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                           int64_t *ret_alignment_position, int64_t *ret_edit_distance) {

  int64_t absolute_reference_id = best_path->get_region_data().reference_id;
  int64_t reference_id = best_path->get_region_data().reference_id;
  int64_t reference_start = index->get_reference_starting_pos()[absolute_reference_id];
  int64_t reference_length = index->get_reference_lengths()[absolute_reference_id];
  SeqOrientation orientation = (best_path->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t best_aligning_position = 0;

  int64_t l1_reference_start = best_path->get_l1_data().l1_lmin;
  int64_t l1_reference_end = ((int64_t) (best_path->get_l1_data().l1_k * read->get_sequence_length())) + best_path->get_l1_data().l1_lmax;
  if (l1_reference_start < reference_start)
    l1_reference_start = reference_start;
  if (l1_reference_end >= (reference_start + reference_length))
    l1_reference_end = reference_start + reference_length - 1;
  int64_t reference_data_length = l1_reference_end - l1_reference_start + 1;

  int64_t alignment_position_end = 0, edit_distance = 0;
  int ret_code = EditDistanceFunction(read->get_data(), read->get_sequence_length(),
                                   (int8_t *) (index->get_data() + l1_reference_start), reference_data_length,
                                   &alignment_position_end, &edit_distance, MYERS_MODE_HW);
  alignment_position_end += l1_reference_start;

  if (ret_code != 0)
    return ret_code;

  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                             FormatString("alignment_position_end = %ld\nl1_reference_start = %ld\n",
                                                          alignment_position_end, l1_reference_start), "LocalRealignmentLinear");
  }

  *ret_alignment_position = alignment_position_end;
  *ret_edit_distance = edit_distance;

  return 0;
}





int LocalRealignmentCircular(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                             int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part, int64_t *ret_nonclipped_left_part,
                             int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part, int64_t *ret_nonclipped_right_part,
                             SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                             int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op, bool perform_reverse_complement) {

  if (best_path->get_region_data().is_split == false) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Called the function for handling the circular part of the genome, but alignment is not split. best_path->region.is_split == false.\n\n"), "LocalRealignmentCircular");
    return -1;
  }

  int8_t *data_copy = NULL;
  int64_t region_length_joined = 0, start_offset = 0, position_of_ref_end = 0;
  ConcatenateSplitRegion(index, (Region &) best_path->get_region_data(), &data_copy, &region_length_joined, &start_offset, &position_of_ref_end);

  int64_t absolute_reference_id = best_path->get_region_data().reference_id;
  int64_t reference_id = best_path->get_region_data().reference_id;
  SeqOrientation orientation = kForward;
  int64_t l1_reference_start = best_path->get_l1_data().l1_lmin;
  int64_t l1_reference_end = ((int64_t) (best_path->get_l1_data().l1_k * read->get_sequence_length())) + best_path->get_l1_data().l1_lmax;
  if (l1_reference_start < 0)
    l1_reference_start = 0;
  if (l1_reference_end >= region_length_joined)
    l1_reference_end = region_length_joined - 1;



  int64_t best_aligning_position_start = 0, best_aligning_position_end = 0, edit_distance = 0;
  std::vector<unsigned char> alignment_left_part;
  int ret_code = AlignmentFunction(read->get_data(), read->get_sequence_length(),
                   (int8_t *) &(data_copy[l1_reference_start]), (l1_reference_end - l1_reference_start + 1),
                   -1, parameters.match_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                   &best_aligning_position_start, &best_aligning_position_end,
                   &edit_distance, alignment_left_part);
  best_aligning_position_start += l1_reference_start;
  best_aligning_position_end += l1_reference_start;

  if (best_path->get_region_data().reference_id >= index->get_num_sequences_forward()) {
    orientation = kReverse;
    if (perform_reverse_complement == true) {
//      read->ReverseComplement();
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                               FormatString("ERROR! Tried to explicitly reverse complement a read! This functionality is no longer allowed!"), "LocalRealignmentCircular");
      exit(1);
    }
    reference_id -= index->get_num_sequences_forward();
  }


  ConvertInsertionsToClipping((unsigned char *) &(alignment_left_part[0]), alignment_left_part.size());
  CountAlignmentOperations(alignment_left_part, read->get_data(), data_copy, reference_id, best_aligning_position_start, orientation,
                           parameters.evalue_match, parameters.evalue_mismatch, parameters.evalue_gap_open, parameters.evalue_gap_extend,
                           ret_eq_op, ret_x_op, ret_i_op, ret_d_op, ret_AS_left_part, ret_nonclipped_left_part);
  
  *ret_AS_right_part = *ret_AS_left_part;
  *ret_nonclipped_right_part = *ret_nonclipped_left_part;



  int64_t best_aligning_position = 0;

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("best_aligning_position_start = %ld\n", best_aligning_position_start), "LocalRealignmentCircular");
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("best_aligning_position_end = %ld\n", best_aligning_position_end), "LocalRealignmentCircular");

  int64_t left_alignment_length = 0, left_alignment_start = 0, left_alignment_end = 0;
  int64_t right_alignment_length = 0, right_alignment_start = 0, right_alignment_end = 0;
  unsigned char *left_alignment = NULL;
  unsigned char *right_alignment = NULL;
  std::vector<unsigned char> alignment_right_part;
  if (ClipCircularAlignment(best_aligning_position_start, best_aligning_position_end, (unsigned char *) &(alignment_left_part[0]), alignment_left_part.size(),
                        (int64_t) (read->get_sequence_length()), (int64_t) (index->get_reference_starting_pos()[absolute_reference_id]),
                        (int64_t) (index->get_reference_lengths()[absolute_reference_id]),
                        start_offset, position_of_ref_end,
                        &left_alignment, &left_alignment_length, &left_alignment_start, &left_alignment_end,
                        &right_alignment, &right_alignment_length, &right_alignment_start, &right_alignment_end) != 0) {
    alignment_left_part.clear();
    alignment_left_part.assign(left_alignment, (left_alignment + left_alignment_length));
    if (left_alignment)
      free(left_alignment);

    alignment_right_part.clear();
    alignment_right_part.assign(right_alignment, (right_alignment + right_alignment_length));
    if (right_alignment)
      free(right_alignment);
  }

  int64_t best_aligning_position_left_part = 0;
  if (alignment_left_part.size() > 0) {
    *ret_cigar_left_part = AlignmentToCigar((unsigned char *) &(alignment_left_part[0]), alignment_left_part.size());
//    *ret_AS_left_part = RescoreAlignment((unsigned char *) &(alignment_left_part[0]), alignment_left_part.size(), parameters.match_score, parameters.mismatch_penalty, parameters.gap_open_penalty, parameters.gap_extend_penalty);

    if (orientation == kForward) {
      index->RawPositionConverterWithRefId(left_alignment_start, reference_id, 0, NULL, &best_aligning_position_left_part, NULL);
    } else {
      index->RawPositionConverterWithRefId(left_alignment_end, reference_id, 0, NULL, &best_aligning_position_left_part, NULL);
      if (perform_reverse_complement == true) {
        *ret_cigar_left_part = ReverseCigarString((*ret_cigar_left_part));
      }
    }
  } else {
    *ret_cigar_left_part = "";
  }

  int64_t best_aligning_position_right_part = 0;
  if (alignment_right_part.size() > 0) {
    *ret_cigar_right_part = AlignmentToCigar((unsigned char *) &(alignment_right_part[0]), alignment_right_part.size());
//    *ret_AS_right_part = RescoreAlignment((unsigned char *) &(alignment_right_part[0]), alignment_right_part.size(), parameters.match_score, parameters.mismatch_penalty, parameters.gap_open_penalty, parameters.gap_extend_penalty);

    if (orientation == kForward) {
      index->RawPositionConverterWithRefId(right_alignment_start, reference_id, 0, NULL, &best_aligning_position_right_part, NULL);
    } else {
      index->RawPositionConverterWithRefId(right_alignment_end, reference_id, 0, NULL, &best_aligning_position_right_part, NULL);
      if (perform_reverse_complement == true) {
        *ret_cigar_right_part = ReverseCigarString((*ret_cigar_right_part));
      }
    }
  } else {
    *ret_cigar_right_part = "";
  }

//  CountAlignmentOperations(alignment, read->get_data(), index->get_data(), reference_id, best_aligning_position_start, orientation, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);

  if (CheckAlignmentSane(alignment_left_part, read, index, reference_id, best_aligning_position_left_part))
    return -1;

  if (CheckAlignmentSane(alignment_right_part, read, index, reference_id, best_aligning_position_right_part))
    return -1;



//#ifndef RELEASE_VERSION
//  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
//    printf ("Reference:\t");
//    PrintSubstring((char *) &(index->get_data()[best_aligning_position + index->get_reference_starting_pos()[reference_id]]), 100, stdout);
//    printf ("\n");
//    printf ("Read:\t\t");
//    PrintSubstring((char *) read->get_data(), 100, stdout);
//    printf ("\n");
//    fflush(stdout);
//  }
//#endif
//
//#ifndef RELEASE_VERSION
//  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read)
//  {
//    printf ("\nqname: %s\n", read->get_header());
//    printf ("CIGAR begining:\t");
//    PrintSubstring((char *) (*ret_cigar).c_str(), 25, stdout);
//    printf ("\n");
//    printf ("CIGAR end:\t");
//    PrintSubstring((char *) &((*ret_cigar).c_str()[ret_cigar->size() - 25]), 25, stdout);
//    printf ("\n");
//    fflush(stdout);
//  }
//#endif

  *ret_alignment_position_left_part = best_aligning_position_left_part;
  *ret_alignment_position_right_part = best_aligning_position_right_part;
  *ret_orientation = orientation;
  *ret_reference_id = reference_id;
  *ret_position_ambiguity = 0;

  if (data_copy)
    delete[] data_copy;

  return ((int) edit_distance);

}



int CalcEditDistanceCircular(EditDistanceFunctionType EditDistanceFunction, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                           int64_t *ret_alignment_position, int64_t *ret_edit_distance) {

  if (best_path->get_region_data().is_split == false) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Called the function for handling the circular part of the genome, but alignment is not split. best_path->region.is_split == false.\n\n"), "LocalRealignmentCircular");
    return -1;
  }

  int8_t *data_copy = NULL;
  int64_t region_length_joined = 0, start_offset = 0, position_of_ref_end = 0;
  ConcatenateSplitRegion(index, best_path->get_region_data(), &data_copy, &region_length_joined, &start_offset, &position_of_ref_end);


  int64_t absolute_reference_id = best_path->get_region_data().reference_id;
  int64_t reference_id = best_path->get_region_data().reference_id;
  SeqOrientation orientation = kForward;
  int64_t l1_reference_start = best_path->get_l1_data().l1_lmin;
  int64_t l1_reference_end = ((int64_t) (best_path->get_l1_data().l1_k * read->get_sequence_length())) + best_path->get_l1_data().l1_lmax;
  if (l1_reference_start < 0)
    l1_reference_start = 0;
  if (l1_reference_end >= region_length_joined)
    l1_reference_end = region_length_joined - 1;



  int64_t best_aligning_position_end = 0, edit_distance = 0;
  int ret_code = EditDistanceFunction(read->get_data(), read->get_sequence_length(),
                                      (int8_t *) &(data_copy[l1_reference_start]), (l1_reference_end - l1_reference_start + 1),
                                      &best_aligning_position_end, &edit_distance, MYERS_MODE_HW);
  best_aligning_position_end += l1_reference_start;

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("best_aligning_position_end = %ld\n", best_aligning_position_end), "LocalRealignmentCircular");

  *ret_alignment_position = best_aligning_position_end;
  *ret_edit_distance = edit_distance;

  if (data_copy)
    delete[] data_copy;

  return 0;
}



int CountAlignmentOperations(std::vector<unsigned char>& alignment, const int8_t *read_data, const int8_t *ref_data, int64_t reference_hit_id, int64_t alignment_position_start, SeqOrientation orientation,
                             int64_t match, int64_t mismatch, int64_t gap_open, int64_t gap_extend,
                             int64_t* ret_eq, int64_t* ret_x, int64_t* ret_i, int64_t* ret_d, int64_t *ret_alignment_score, int64_t *ret_nonclipped_length) {
  unsigned char last_move = -1;  // Code of last move.
  int64_t num_same_moves = 0;
  int64_t read_position = 0;
  int64_t ref_position = 0;

  int64_t num_eq = 0;
  int64_t num_x = 0;
  int64_t num_i = 0;
  int64_t num_d = 0;
  int64_t alignment_score = 0;

  int64_t nonclipped_length = 0;

  for (int i = 0; i < alignment.size(); i++) {
    char align_op = 255;
    align_op = alignment[i];

//    if (align_op == EDLIB_S) {
//      read_position += 1;
//
//    } else
    if (align_op == EDLIB_M || align_op == EDLIB_EQUAL || align_op == EDLIB_X) {
      if (read_data[read_position] == ref_data[alignment_position_start + ref_position]) {
        num_eq += 1;
        alignment_score += match;
//        printf ("=");
//        fflush(stdout);

      } else {
        num_x += 1;
        alignment_score -= mismatch;
//        printf ("X");
//        fflush(stdout);
      }

      nonclipped_length += 1;

    } else if (align_op == EDLIB_I) {
      num_i += 1;
      // This is for the (gap_open + (N) * gap_extend).
//      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open + gap_extend) : gap_extend);
      // This is for the (gap_open + (N - 1) * gap_extend).
      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open) : gap_extend);
      nonclipped_length += 1;
//      printf ("I");
//      fflush(stdout);

    } else if (align_op == EDLIB_D) {
      num_d += 1;
      // This is for the (gap_open + (N) * gap_extend).
//      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open + gap_extend) : gap_extend);
      // This is for the (gap_open + (N - 1) * gap_extend).
      alignment_score -= ((i == 0 || (i > 0 && alignment[i-1] != align_op)) ? (gap_open) : gap_extend);
//      printf ("D");
//      fflush(stdout);
    }

    // Increase coordinates.
    if (align_op == EDLIB_M || align_op == EDLIB_EQUAL || align_op == EDLIB_X || align_op == EDLIB_I || align_op == EDLIB_S)
      read_position += 1;
    if (align_op == EDLIB_M || align_op == EDLIB_EQUAL || align_op == EDLIB_X || align_op == EDLIB_D)
      ref_position += 1;
  }

//  printf ("\n\n");
//  fflush(stdout);

  *ret_eq = num_eq;
  *ret_x = num_x;
  *ret_i = num_i;
  *ret_d = num_d;
  *ret_alignment_score = alignment_score;
  *ret_nonclipped_length = nonclipped_length;

  return 0;
}



// Checks if there is a strange (large) number of insertions and deletions, or consecutive insertion/deletion operations.
// SeqAn likes to make such alignments.
// Returns 0 if everything went ok.
int CheckAlignmentSane(std::vector<unsigned char> &alignment, const SingleSequence* read, const Index* index, int64_t reference_hit_id, int64_t reference_hit_pos) {
  unsigned char last_move = -1;  // Code of last move.
  int64_t num_same_moves = 0;
  int64_t read_length = 0;
  int64_t ref_length = 0;

  for (int i = 0; i <= alignment.size(); i++) {
    char alignment_char = 255;
    if (i < alignment.size()) {
      alignment_char = alignment[i];
      if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X || alignment[i] == EDLIB_I || alignment[i] == EDLIB_S)
        read_length += 1;
      if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X || alignment[i] == EDLIB_D)
        ref_length += 1;
    }

    // If new sequence of same moves started
    if (i == alignment.size() || alignment_char != last_move) {
        if (i > 0) {  // if previous sequence of same moves ended
          // num_same_moves
          // last_move == EDLIB_M
          // If the number of consecutive insertions or deletions is very high, something is wrong.
//          if ((last_move == EDLIB_I || last_move == EDLIB_D) && num_same_moves >= 99) {
//          if ((last_move == EDLIB_I || last_move == EDLIB_D)) {
//            printf ("CheckAlignmentSane returned false! return 1.\n");
//            fflush(stdout);
//            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("CheckAlignmentSane returned false! return 1.\n"), "CheckAlignmentSane");
//            return 1;
//          }
          // If there are insertions following deletions (or other way around), something is wrong again.
          if ((last_move == EDLIB_I && alignment_char == EDLIB_D) || (last_move == EDLIB_D && alignment_char == EDLIB_I)) {
//            printf ("CheckAlignmentSane returned false! return 2.\n");
//            fflush(stdout);
            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("CheckAlignmentSane returned false! return 2. num_same_moves = %ld, qname: '%s', read_length: %ld, ref_length: %ld\n", num_same_moves, read->get_header(), read_length, ref_length), "CheckAlignmentSane");
            return 2;
          }
        }
        if (i < alignment.size()) {
            num_same_moves = 0;
            last_move = alignment_char;
        }
    }
    num_same_moves++;
  }

  if (read != NULL && read_length != read->get_sequence_length()) {
//    printf ("CheckAlignmentSane returned false! return 3.\n");
//    printf ("read_length = %ld\n", read_length);
//    printf ("read->get_sequence_length() = %ld\n", read->get_sequence_length());
//    fflush(stdout);
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("CheckAlignmentSane returned false! return 3. read_length = %ld, read->get_sequence_length() = %ld\\n", read_length, read->get_sequence_length()), "CheckAlignmentSane");
    return 3;
  }
  if ((index != NULL && reference_hit_id >= 0 && reference_hit_pos >= 0) && ref_length > index->get_reference_lengths()[reference_hit_id]) {
//    printf ("CheckAlignmentSane returned false! return 4.\n");
//    fflush(stdout);
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("CheckAlignmentSane returned false! return 4.\n"), "CheckAlignmentSane");
    return 4;
  }
  if ((index != NULL && reference_hit_id >= 0 && reference_hit_pos >= 0) &&
      (reference_hit_pos + ref_length) > (index->get_reference_starting_pos()[reference_hit_id] + index->get_reference_lengths()[reference_hit_id])) {
//    printf ("CheckAlignmentSane returned false! return 5.\n");
//    fflush(stdout);
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("CheckAlignmentSane returned false! return 5.\n"), "CheckAlignmentSane");
    return 5;
  }

  return 0;
}
