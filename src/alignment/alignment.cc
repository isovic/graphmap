/*
 * alignment.cc
 *
 *  Created on: Jan 17, 2016
 *      Author: isovic
 */

#include "alignment/alignment.h"


int AlignRegion(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, bool extend_to_end, PathGraphEntry *region_results) {

    if (parameters->alignment_algorithm == "gotoh") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh for alignment!\n", "Alignment");
//      if (region_results->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, region_results, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//      else
//        return LocalRealignmentCircular(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, region_results, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);

    } else if (parameters->alignment_algorithm == "myers") {

      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using semiglobal alignment approach.\n", "Alignment");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");

      return SemiglobalAlignment(MyersSemiglobalWrapper, read, index, parameters, evalue_params, region_results);

//    } else if (parameters.alignment_algorithm == "anchor") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//      bool is_linear = region_results->get_region_data().is_split == false || parameters.is_reference_circular == false;
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignment(is_linear, extend_to_end, MyersNWWrapper, MyersSHWWrapper, read, index, parameters, region_results, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//
//    } else if (parameters.alignment_algorithm == "anchorgotoh") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//      bool is_linear = region_results->get_region_data().is_split == false || parameters.is_reference_circular == false;
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Gotoh's algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignment(is_linear, extend_to_end, SeqAnNWWrapper, SeqAnSHWWrapper, read, index, parameters, region_results, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//
//#ifndef RELEASE_VERSION
//    } else if (parameters.alignment_algorithm == "anchormex") {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using anchored alignment approach.\n", "Alignment");
//      /// Extension to read ends is not currently supported. SHW alignment needs to be implemented in Opal first.
//      extend_to_end = false;
//      bool is_linear = region_results->get_region_data().is_split == false || parameters.is_reference_circular == false;
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using Match Extend algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignmentMex(is_linear, extend_to_end, OpalNWWrapper, OpalSHWWrapper, read, index, parameters, region_results, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//#endif

    } else {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Warning: Unknown alignment algorithm selected. Using Myers' bit-vector alignment instead.\n", "Alignment");
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Using semiglobal alignment approach.\n", "Alignment");
//      if (region_results->get_region_data().is_split == false || parameters.is_reference_circular == false)
//        return LocalRealignmentLinear(MyersSemiglobalWrapper, read, index, parameters, region_results, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
//      else
//        return LocalRealignmentCircular(MyersSemiglobalWrapper, read, index, parameters, region_results, ret_alignment_position_left_part, ret_cigar_left_part, ret_AS_left_part, ret_nonclipped_left_part, ret_alignment_position_right_part, ret_cigar_right_part, ret_AS_right_part, ret_nonclipped_right_part, ret_orientation, ret_reference_id, ret_position_ambiguity, ret_eq_op, ret_x_op, ret_i_op, ret_d_op);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using semiglobal alignment approach.\n", "Alignment");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");

      return SemiglobalAlignment(MyersSemiglobalWrapper, read, index, parameters, evalue_params, region_results);
    }

  return -10;
}

void VerboseAlignment(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const AlignmentResults *aln) {
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    std::string alignment_as_string = "";

    alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data(), read->get_sequence_length(),
       (const unsigned char *) (index->get_data() + aln->aln_window_start), (aln->aln_window_end - aln->aln_window_start + 1),
       (unsigned char *) &(aln->raw_alignment[0]), aln->raw_alignment.size(),
       (aln->raw_pos_end - aln->aln_window_start), aln->aln_mode_code);

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read,
       FormatString("alignment_position_start = %ld\nalignment_position_end = %ld\nl1_reference_start = %ld\nedit_distance = %ld\n%s\n",
       aln->raw_pos_start, aln->raw_pos_end, aln->aln_window_start, aln->edit_distance, alignment_as_string.c_str()), std::string(__FUNCTION__));
  }
}

/// Determines the start and end locations for semiglobal alignment, keeping in mind the boundaries of the reference being aligned to. Works with circular alignment as well.
int CalcL1AlignmentBounds(const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const PathGraphEntry *region_results, int64_t *ret_l1_start, int64_t *ret_l1_end, int64_t *ret_l1_len) {
  SeqOrientation orientation = (region_results->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t abs_ref_id = region_results->get_region_data().reference_id;
  int64_t ref_id = (orientation == kForward) ? (abs_ref_id) : (abs_ref_id - index->get_num_sequences_forward());
  int64_t ref_start = index->get_reference_starting_pos()[abs_ref_id];
  int64_t ref_length = index->get_reference_lengths()[abs_ref_id];
  int64_t best_aligning_position = 0;

  int64_t l1_ref_start = region_results->get_l1_data().l1_lmin;
  int64_t l1_ref_end = ((int64_t) (region_results->get_l1_data().l1_k * read->get_sequence_length())) + region_results->get_l1_data().l1_lmax;
  if (l1_ref_start < ref_start)
    l1_ref_start = ref_start;
  if (l1_ref_end >= (ref_start + ref_length))
    l1_ref_end = ref_start + ref_length - 1;
  int64_t reference_data_length = l1_ref_end - l1_ref_start + 1;

  *ret_l1_start = l1_ref_start;
  *ret_l1_end = l1_ref_end;
  *ret_l1_len = reference_data_length;

  /// Debug output.
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read,
              FormatString("\nl1_reference_start = %ld\nl1_reference_end = %ld\nreference_data_length = %ld\nreference_start = %ld\nreference_length = %ld\nabsolute_reference_id = %ld\norientation = %s\n\n",
              l1_ref_start, l1_ref_end, reference_data_length, ref_start, ref_length, abs_ref_id, ((orientation == kForward) ? "forward" : "reverse")), std::string(__FUNCTION__));
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, FormatString("\n%s\n", region_results->VerboseInfoToString().c_str()), std::string(__FUNCTION__));
  }

  return 0;
}

int SemiglobalAlignment(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, PathGraphEntry *region_results) {
  SeqOrientation orientation = (region_results->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t abs_ref_id = region_results->get_region_data().reference_id;
  int64_t ref_id = (orientation == kForward) ? (abs_ref_id) : (abs_ref_id - index->get_num_sequences_forward());

  int64_t l1_ref_start = 0, l1_ref_end = 0, reference_data_length = 0;
  CalcL1AlignmentBounds(read, index, parameters, region_results, &l1_ref_start, &l1_ref_end, &reference_data_length);

  /// Perform alignment and store results.
  AlignmentResults aln;
  std::vector<unsigned char> alignment1;
  int64_t aln_pos_start = 0, aln_pos_end = 0;
  int ret_code = AlignmentFunction(read->get_data(), read->get_sequence_length(),
                                   (int8_t *) (index->get_data() + l1_ref_start), reference_data_length,
                                   -1, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
                                   &aln_pos_start, &aln_pos_end,
                                   &aln.edit_distance,
                                   aln.raw_alignment);
  /// Sanity check.
  if (ret_code != 0 || aln.raw_alignment.size() == 0) { return ret_code; }

  aln_pos_start += l1_ref_start;
  aln_pos_end += l1_ref_start;
  aln.raw_alignment = FixAlignment((unsigned char *) &(aln.raw_alignment), aln.raw_alignment.size());
  ConvertInsertionsToClipping((unsigned char *) &(aln.raw_alignment[0]), aln.raw_alignment.size());

  /// This part converts the global alignment coordinates to local. 'Global' meaning the coordinates on the entire data array from the index.
  /// 'Local' meaning the coordinates on the reference that was hit (in range [0, ref_len>).
  /// final_aln_pos_start is the alignment position within the reference (local to the reference), and is also reversed (subtracted from reference_length) in case orientation == kReverse.
  int64_t final_aln_pos_start = 0, final_aln_pos_end = 0;
  index->RawPositionConverterWithRefId((orientation == kForward) ? aln_pos_start : aln_pos_end, abs_ref_id, 0, NULL, &final_aln_pos_start, NULL);
  index->RawPositionConverterWithRefId((orientation == kForward) ? aln_pos_end : aln_pos_start, abs_ref_id, 0, NULL, &final_aln_pos_end, NULL);

  /// Assign resulting values.
  // aln.raw_alignment is assigned directly in the AlignmentFunction call to avoid copying of the data twice.
  // aln.edit_distance is assigned directly in the AlignmentFunction call.
  aln.orientation = orientation;
  aln.is_reverse = (orientation == kReverse) ? true : false;
  aln.pos_start = final_aln_pos_start;
  aln.pos_end = final_aln_pos_end; // aln.pos_start + (aln_pos_end - aln_pos_start);
  aln.raw_pos_start = aln_pos_start;
  aln.raw_pos_end = aln_pos_end;
  aln.ref_id = ref_id;
  aln.ref_header = index->get_headers()[ref_id];
  aln.ref_len = index->get_reference_lengths()[ref_id];
  aln.query_id = read->get_sequence_absolute_id();
  aln.query_header = read->get_header();
  aln.query_len = read->get_sequence_length();
  aln.is_aligned = true;
  aln.aln_window_start = l1_ref_start;
  aln.aln_window_end = l1_ref_end;
  aln.aln_mode_code = MYERS_MODE_HW;

  region_results->AddAlignmentData(aln);

  /// Verbose debug.
  VerboseAlignment(read, index, parameters, &aln);

  /// Fill out statistics for each alignment (E-value calculation, couting of CIGAR operations, etc.) and check if the alignments are sane.
  for (int32_t i=0; i<region_results->get_alignments().size(); i++) {
    AlignmentResults *curr_aln = &region_results->get_alignments()[i];

    /// In case alignment needs to be reversed now. This is because previously we mapped to the reverse complement of the reference, but in the output the read is reversed, so the alignment is backwards.
    if (curr_aln->orientation == kForward) {
      curr_aln->cigar = AlignmentToCigar((unsigned char *) &(curr_aln->raw_alignment[0]), curr_aln->raw_alignment.size(), parameters->use_extended_cigar);
      curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), index->get_data(), ref_id, aln_pos_start);

    } else {
      std::vector<int8_t> alignment_reverse (curr_aln->raw_alignment.rbegin(), curr_aln->raw_alignment.rend());
      curr_aln->cigar = AlignmentToCigar((unsigned char *) &(alignment_reverse[0]), alignment_reverse.size(), parameters->use_extended_cigar);
      curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) alignment_reverse, read->get_data(), index->get_data(), ref_id, aln_pos_start);
    }

    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), index->get_data(), ref_id, aln_pos_start, orientation,
                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend,
                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->nonclipped_length);
    CalculateEValueDNA(curr_aln->alignment_score, curr_aln->nonclipped_length, index->get_data_length_forward(), evalue_params, &curr_aln->evalue);
    if (CheckAlignmentSane((std::vector<unsigned char> &) curr_aln->raw_alignment, read, index, curr_aln->ref_id, curr_aln->pos_start) != 0) {
      curr_aln->is_aligned = false;
    }
  }

  /// Count the number of 'unaligned' alignments.
  int32_t num_unaligned = 0;
  for (int32_t i=0; i<region_results->get_alignments().size(); i++) { if (region_results->get_alignments()[i].is_aligned == false) num_unaligned += 1; }
  if (num_unaligned == region_results->get_alignments().size()) { return ALIGNMENT_NOT_SANE; }

  return 0;
}

//Ovo je iz GenerateAlignments_:     primary_alignment.pos_start = relative_position_left_part + 1;
//Sto znaci da cu to trebati handleati drukcije drugdje, prilikom zapisa SAM filea!


