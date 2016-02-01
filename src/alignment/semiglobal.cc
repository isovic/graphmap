/*
 * semiglobal.cc
 *
 *  Created on: Jan 27, 2016
 *      Author: isovic
 */

#include "alignment/alignment.h"

// pos_of_ref_end is the distance from the beginning of the region to the last base of the reference, before the beginning of the reference is concatenated.
int SplitCircularAlignment(const AlignmentResults *aln, int64_t pos_of_ref_end, AlignmentResults *aln_l, AlignmentResults *aln_r) {
  *aln_l = *aln;
  *aln_r = *aln;

  if (aln->reg_pos_start > pos_of_ref_end) {
    // The entire alignment is on the right part of the circular region.
    aln_l->is_aligned = false;
    return 0;

  } else if (aln->reg_pos_end <= pos_of_ref_end) {
    // The entire alignment is on the left part of the circular region.
    aln_r->is_aligned = false;
    return 0;
  }

  int64_t end_on_aln = 0, end_on_read = 0, end_on_ref = 0;
  int64_t start_on_aln = 0, start_on_read = 0, start_on_ref = 0;

  FindCircularEnd(aln->raw_alignment, pos_of_ref_end - aln->reg_pos_start,
                  &end_on_aln, &end_on_read, &end_on_ref,
                  &start_on_aln, &start_on_read, &start_on_ref);

  // Trim the left part of the alignment.
  aln_l->raw_pos_end = aln_l->raw_pos_start + end_on_ref - 1;
  LOG_DEBUG("pos_of_ref_end = %ld\n", pos_of_ref_end);
  LOG_DEBUG("pos_of_ref_end - aln->reg_pos_start = %ld\n", pos_of_ref_end - aln->reg_pos_start);
  LOG_DEBUG("end_on_aln = %ld\n", end_on_aln);
  LOG_DEBUG("end_on_read = %ld\n", end_on_read);
  LOG_DEBUG("end_on_ref = %ld\n", end_on_ref);
  LOG_DEBUG("start_on_aln = %ld\n", start_on_aln);
  LOG_DEBUG("start_on_read = %ld\n", start_on_read);
  LOG_DEBUG("start_on_ref = %ld\n", start_on_ref);
  LOG_DEBUG("aln_l->query_len - (end_on_read + 1) = %ld\n", aln_l->query_len - (end_on_read + 1));

  int64_t clip_l = aln_l->query_len - (end_on_read + 1) + 1;
  uint8_t *alignment_l = (uint8_t *) malloc(sizeof(uint8_t) * ((end_on_aln + 1) + clip_l + 1));
  memmove(&alignment_l[0], &aln->raw_alignment[0], end_on_aln + 1);
  memset(&alignment_l[end_on_aln+1], EDLIB_S, clip_l);
  aln_l->raw_alignment.clear();
  aln_l->raw_alignment.assign(alignment_l, alignment_l + ((end_on_aln + 1) + clip_l));

//  // Trim the right part of the alignment.
//  aln_r->raw_pos_start = aln_r->raw_pos_start + end_on_ref - 1;
//  int64_t clip_l = aln_r->query_len - (end_on_read + 1) + 1;
//  uint8_t *alignment_l = (uint8_t *) malloc(sizeof(uint8_t) * ((end_on_aln + 1) + clip_l + 1));
//  memmove(&alignment_l[0], &aln->raw_alignment[0], end_on_aln + 1);
//  memset(&alignment_l[end_on_aln+1], EDLIB_S, clip_l);
//  aln_r->raw_alignment.clear();
//  aln_r->raw_alignment.assign(alignment_l, alignment_l + ((end_on_aln + 1) + clip_l));



  printf ("Before:\n");
  for (int64_t i=0; i<aln->raw_alignment.size(); i++) {
    printf ("%d", aln->raw_alignment[i]);
  }
  printf ("\n");
  printf ("\n");
  fflush(stdout);

  printf ("After:\n");
  for (int64_t i=0; i<aln_l->raw_alignment.size(); i++) {
    printf ("%d", aln_l->raw_alignment[i]);
  }
  printf ("\n");
  printf ("\n");
  fflush(stdout);

//  LOG_DEBUG_SPEC("end_on_read = %ld\n", end_on_read);
//  LOG_DEBUG_SPEC("Read trim:\n%s\n", GetSubstring((char *) (read->get_data() + ((end_on_read + 1))), read->get_data_length() - (end_on_read + 1)).c_str());

  return 0;
}

int SemiglobalAlignment(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, PathGraphEntry *region_results) {
  /// General useful things.
  const Region &region = region_results->get_region_data();
  SeqOrientation orientation = (region_results->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t abs_ref_id = region_results->get_region_data().reference_id;
  int64_t ref_id = (orientation == kForward) ? (abs_ref_id) : (abs_ref_id - index->get_num_sequences_forward());
  int64_t ref_start = index->get_reference_starting_pos()[region.reference_id];
  int64_t ref_len = index->get_reference_lengths()[region.reference_id];

  int8_t *reg_data = NULL;       // The data of the region.
  int64_t reg_data_len = 0;         // Length of the region_data.
  int64_t index_pos = 0;            // Position of the start of the region on the original Index data. E.g. reg_data[0] is the same base as index->get_data()[index_pos].
  int64_t pos_of_ref_end = 0; // If the region was circular, it crosses the boundary between the end and the start of the data. ref_data[index_pos_of_ref_end] is the last base of the reference before the split part is concatenated.
  bool is_cleanup_required = NULL;  // If true, the region_data will have to be freed manually.
  int64_t l1_start = 0;             // The start position of the L1 determined boundaries, but within the region (relative to region start).
  int64_t l1_end = 0;               // The end position of the L1 determined boundaries, but within the region (relative to region start).
  // Here, a pointer to the beginning of the region is obtained. If the region was linear, reg_data just points to a location in the Index.
  // If the region was circular, reg_data points to a newly allocated memory for the concatenated region, and needs to be cleared.
  int retval_region = GetRegionData(index, &region, &reg_data, &reg_data_len, &index_pos, &pos_of_ref_end, &is_cleanup_required);
  GetL1PosInRegion(read, index, parameters, region_results, &l1_start, &l1_end);

  if (retval_region) {
    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "GetRegionData returned with error!"));
    if (is_cleanup_required) { free(reg_data); }
    return 1;
  }

  if (is_cleanup_required) { LOG_DEBUG_SPEC("Alignment is circular and manual cleanup will be required.\n"); }

  /// Perform alignment and store results.
  AlignmentResults aln;
  std::vector<unsigned char> alignment1;
  int64_t aln_pos_start = 0, aln_pos_end = 0;
  int ret_code = AlignmentFunction(read->get_data(), read->get_sequence_length(),
                                   reg_data + l1_start, (int64_t) (l1_end - l1_start),
                                   -1, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
                                   &aln_pos_start, &aln_pos_end,
                                   &aln.edit_distance,
                                   aln.raw_alignment);

  /// Sanity check.
  if (ret_code != 0 || aln.raw_alignment.size() == 0) {
    if (is_cleanup_required) { free(reg_data); }
    LOG_DEBUG_SPEC("Something went wrong with AlignmentFunction.\n");
    return ret_code;
  }

  int64_t aln_start_in_reg = l1_start + aln_pos_start;
  int64_t aln_end_in_reg = l1_start + aln_pos_end;

  int64_t aln_abs_start = aln_start_in_reg + index_pos;
  int64_t aln_abs_end = aln_end_in_reg + index_pos;

  /// Assign resulting values.
  // aln.raw_alignment is assigned directly in the AlignmentFunction call to avoid copying of the data twice.
  // aln.edit_distance is assigned directly in the AlignmentFunction call.
//  aln.ref_data = ref_data;
//  aln.is_cleanup_required = is_cleanup_required;
  aln.orientation = orientation;
  aln.is_reverse = (orientation == kReverse) ? true : false;
  aln.ref_start = -1;
  aln.ref_end = -1;
  aln.query_start = -1;
  aln.query_end = -1;
  aln.raw_pos_start = aln_abs_start;
  aln.raw_pos_end = aln_abs_end;
  aln.reg_pos_start = aln_start_in_reg;
  aln.reg_pos_end = aln_end_in_reg;
  aln.ref_id = ref_id;
  aln.ref_header = index->get_headers()[ref_id];
  aln.ref_len = index->get_reference_lengths()[ref_id];
  aln.query_id = read->get_sequence_absolute_id();
  aln.query_header = read->get_header();
  aln.query_len = read->get_sequence_length();
  aln.is_aligned = true;
  aln.aln_mode_code = MYERS_MODE_HW;
  aln.alignment = aln.raw_alignment;

  /// If alignment is linear, just add it. Otherwise, it needs to be split first.
  if (region.is_split == false || parameters->is_reference_circular == false) {
    if (orientation == kReverse) ReverseArray(aln.alignment);
    region_results->AddAlignmentData(aln);

  } else {
    LOG_DEBUG_SPEC("Alignment is circular and will be split.\n");

    AlignmentResults aln_l;
    AlignmentResults aln_r;
    SplitCircularAlignment(&aln, pos_of_ref_end, &aln_l, &aln_r);
    region_results->AddAlignmentData(aln_l);
//    region_results->AddAlignmentData(aln_r);

//
//    int64_t end_on_aln = 0, end_on_read = 0, end_on_ref = 0;
//    int64_t start_on_aln = 0, start_on_read = 0, start_on_ref = 0;
//    FindCircularEnd(aln.raw_alignment, pos_of_ref_end,
//                    &end_on_aln, &end_on_read, &end_on_ref,
//                    &start_on_aln, &start_on_read, &start_on_ref);
//
////    aln.ref_end = end_on_ref;
////    aln.query_end = end_on_read;
//    aln.ref_end = aln.ref_start + end_on_ref;
//    aln.raw_pos_end = aln.raw_pos_start + end_on_ref;
//    uint8_t *alignment_l = (uint8_t *) malloc(sizeof(uint8_t) * (end_on_aln + (read->get_sequence_length() - (end_on_read + 1)) + 1));
//    memmove(&alignment_l[0], &aln.raw_alignment[0], end_on_aln + 1);
//    memset(&alignment_l[end_on_aln+1], EDLIB_S, (read->get_sequence_length() - (end_on_read + 1)  ));
//    LOG_DEBUG_SPEC("end_on_read = %ld\n", end_on_read);
//    LOG_DEBUG_SPEC("Read trim:\n%s\n", GetSubstring((char *) (read->get_data() + ((end_on_read + 1))), read->get_data_length() - (end_on_read + 1)).c_str());
//
//    aln_r.ref_start = start_on_ref;
//    aln_r.query_start = start_on_read;
//    uint8_t *alignment_r = (uint8_t *) malloc(sizeof(uint8_t) * ((aln.raw_alignment.size() - start_on_aln) + (start_on_read) + 1));
//    memset(&alignment_r[0], EDLIB_S, (start_on_read));
//    memmove(&alignment_r[start_on_read], &aln.raw_alignment[start_on_aln], (aln.raw_alignment.size() - start_on_aln));
//
////    printf ("Before:\n");
////    for (int64_t i=0; i<aln.raw_alignment.size(); i++) {
////      printf ("%d", aln.raw_alignment[i]);
////    }
////    printf ("\n");
////    printf ("\n");
////    fflush(stdout);
//
//    aln.raw_alignment.assign(alignment_l, alignment_l + (end_on_aln + (read->get_sequence_length() - (end_on_read + 1))));
//    aln.alignment = aln.raw_alignment;
//    if (orientation == kReverse) ReverseArray(aln.alignment);
//    // ref_data points to the start of the region (either on the index->data or on the newly malloced array). win_start points to the l1_start location which is located in ref_data.
//    // offset_from_ref_start is the distance from ref_data to the start of the reference. ref_start is the position of the beginning of the reference within index->data
//    aln.aln_window_start = (win_start - ref_data) + offset_from_ref_start + ref_start;
//    aln.aln_window_end = aln.aln_window_start + ((int64_t) (end_on_ref));
//    LOG_DEBUG_SPEC("aln.aln_window_start = %ld\n", aln.aln_window_start);
//    LOG_DEBUG_SPEC("aln.aln_window_end = %ld\n", aln.aln_window_end);
//
//    aln_r.raw_alignment.assign(alignment_r, alignment_r + ((aln.raw_alignment.size() - start_on_aln) + (start_on_read)));
//    aln_r.alignment = aln_r.raw_alignment;
//    if (orientation == kReverse) ReverseArray(aln_r.alignment);
//    aln_r.aln_window_start = start_on_ref + ref_start;
//    aln_r.aln_window_end = aln_r.aln_window_start + ((int64_t) (win_end - win_start - start_on_ref)) - 1;
//    LOG_DEBUG_SPEC("aln_r.aln_window_start = %ld\n", aln_r.aln_window_start);
//    LOG_DEBUG_SPEC("aln_r.aln_window_end = %ld\n", aln_r.aln_window_end);
//
//    printf ("After:\n");
//    for (int64_t i=0; i<aln.raw_alignment.size(); i++) {
//      printf ("%d", aln.raw_alignment[i]);
//    }
//    printf ("\n");
//    printf ("\n");
//    fflush(stdout);
//    LOG_DEBUG_SPEC("read:\n%s\n", GetSubstring((char *) (read->get_data()), 15).c_str());
//    LOG_DEBUG_SPEC("Ref:\n%s\n", GetSubstring((char *) (index->get_data() + aln.aln_window_start), 15).c_str());
//    LOG_DEBUG_SPEC("ref_data:\n%s\n", GetSubstring((char *) (win_start), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start:\n%s\n", GetSubstring((char *) (index->get_data() + aln_pos_start), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start2:\n%s\n", GetSubstring((char *) (win_start + aln_pos_start - (ref_start + offset_from_ref_start)), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start3:\n%s\n", GetSubstring((char *) (index->get_data() + ref_start + offset_from_ref_start + (win_start - ref_data) + aln_pos_start - (ref_start + offset_from_ref_start)), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start4:\n%s\n", GetSubstring((char *) (index->get_data() + aln_pos_start + (win_start - ref_data)), 15).c_str());
////    exit(1);
//
//    region_results->AddAlignmentData(aln);
//    region_results->AddAlignmentData(aln_r);
  }

  /// Fill out statistics for each alignment (E-value calculation, couting of CIGAR operations, etc.) and check if the alignments are sane.
  for (int32_t i=0; i<region_results->get_alignments().size(); i++) {
    AlignmentResults *curr_aln = &region_results->get_alignments()[i];

    /// Verbose debug.
    LOG_DEBUG_SPEC("Alignment part %d / %d:\n", i, region_results->get_alignments().size());

    ConvertInsertionsToClipping((unsigned char *) &(curr_aln->raw_alignment[0]), curr_aln->raw_alignment.size());
    int64_t num_clipped_front = 0, num_clipped_back = 0;
    CountClippedBases((unsigned char *) &(curr_aln->raw_alignment[0]), curr_aln->raw_alignment.size(), &num_clipped_front, &num_clipped_back);
    curr_aln->query_start = (curr_aln->orientation == kForward) ? num_clipped_front : num_clipped_back;
    curr_aln->query_end = read->get_sequence_length() - ((curr_aln->orientation == kForward) ? num_clipped_back : num_clipped_front) - 1;

    /// This part converts the global alignment coordinates to local. 'Global' meaning the coordinates on the entire data array from the index.
    /// 'Local' meaning the coordinates on the reference that was hit (in range [0, ref_len>).
    /// final_aln_pos_start is the alignment position within the reference (local to the reference), and is also reversed (subtracted from reference_length) in case orientation == kReverse.
    int64_t final_aln_pos_start = 0, final_aln_pos_end = 0;
    index->RawPositionConverterWithRefId((curr_aln->orientation == kForward) ? curr_aln->raw_pos_start : curr_aln->raw_pos_end, abs_ref_id, 0, NULL, &final_aln_pos_start, NULL);
    index->RawPositionConverterWithRefId((curr_aln->orientation == kForward) ? curr_aln->raw_pos_end : curr_aln->raw_pos_start, abs_ref_id, 0, NULL, &final_aln_pos_end, NULL);
    curr_aln->ref_start = final_aln_pos_start % ref_len;
    curr_aln->ref_end = final_aln_pos_end % ref_len;

    LOG_DEBUG_SPEC("Converting alignment to CIGAR string.\n");
    curr_aln->cigar = AlignmentToCigar((unsigned char *) &(curr_aln->alignment[0]), curr_aln->alignment.size(), parameters->use_extended_cigar);

    LOG_DEBUG_SPEC("Converting alignment to MD string.\n");
    curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) curr_aln->alignment, read->get_data(), index->get_data(), ref_id, curr_aln->ref_start + index->get_reference_starting_pos()[curr_aln->ref_id]);

    LOG_DEBUG_SPEC("Counting alignment operations.\n");
    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), reg_data, ref_id, curr_aln->reg_pos_start, orientation,
                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend,
                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->nonclipped_length);

    LOG_DEBUG_SPEC("Calculating the E-value.\n");
    CalculateEValueDNA(curr_aln->alignment_score, curr_aln->nonclipped_length, index->get_data_length_forward(), evalue_params, &curr_aln->evalue);

    LOG_DEBUG_SPEC("Checking if the alignment is sane.\n");
    if (CheckAlignmentSane((std::vector<unsigned char> &) curr_aln->raw_alignment, read, index, curr_aln->ref_id, curr_aln->ref_start) != 0) {
      curr_aln->is_aligned = false;
      LOG_DEBUG_SPEC("Alignment is insane!\n");
    }

    VerboseAlignment(read, index, parameters, curr_aln);
  }

  if (is_cleanup_required) { free(reg_data); }

  /// Count the number of 'unaligned' alignments.
  int32_t num_unaligned = 0;
  for (int32_t i=0; i<region_results->get_alignments().size(); i++) { if (region_results->get_alignments()[i].is_aligned == false) num_unaligned += 1; }
  if (num_unaligned == region_results->get_alignments().size()) { return ALIGNMENT_NOT_SANE; }

  return 0;
}

//int SemiglobalAlignment2(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, PathGraphEntry *region_results) {
//  /// General useful things.
//  const Region &region = region_results->get_region_data();
//  SeqOrientation orientation = (region_results->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
//  int64_t abs_ref_id = region_results->get_region_data().reference_id;
//  int64_t ref_id = (orientation == kForward) ? (abs_ref_id) : (abs_ref_id - index->get_num_sequences_forward());
//  int64_t ref_start = index->get_reference_starting_pos()[region.reference_id];
//
//  /// Get the pointers to the data for alignment.
//  int8_t *ref_data = NULL, *win_start = NULL, *win_end = NULL;
//  int64_t ref_data_len = 0, offset_from_ref_start = 0, pos_of_ref_end = 0;
//  bool is_cleanup_required = NULL;
//  int ret_windata = GetAlignmentWindowData(read, index, parameters, region_results,
//                                           &ref_data, &ref_data_len, &win_start, &win_end,
//                                           &offset_from_ref_start, &pos_of_ref_end, &is_cleanup_required);
//  if (ret_windata) {
//    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "GetAlignmentWindowData returned with error!"));
//    if (is_cleanup_required) { free(ref_data); }
//    return 1;
//  }
//
//  if (is_cleanup_required) {
//    LOG_DEBUG_SPEC("Alignment is circular and manual cleanup will be required.\n");
//  }
//
////  LOG_DEBUG_SPEC("\n%s\n", (char *) GetSubstring((char *) win_start, (win_end - win_start - 1)).c_str());
////  LOG_DEBUG_SPEC("\nRef data:\n%s\n", (char *) GetSubstring((char *) ref_data, ref_data_len).c_str());
//
//  /// Perform alignment and store results.
//  AlignmentResults aln;
//  std::vector<unsigned char> alignment1;
//  int64_t aln_pos_start = 0, aln_pos_end = 0;
//  int ret_code = AlignmentFunction(read->get_data(), read->get_sequence_length(),
//                                   win_start, (int64_t) (win_end - win_start),
//                                   -1, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
//                                   &aln_pos_start, &aln_pos_end,
//                                   &aln.edit_distance,
//                                   aln.raw_alignment);
//
//  /// Sanity check.
//  if (ret_code != 0 || aln.raw_alignment.size() == 0) {
//    if (is_cleanup_required) { free(ref_data); }
//    return ret_code;
//  }
//
//  aln_pos_start += (ref_start + offset_from_ref_start) + (win_start - ref_data - 1);
//  aln_pos_end += (ref_start + offset_from_ref_start) + (win_start - ref_data - 1);
//
////  aln.raw_alignment = FixAlignment((unsigned char *) &(aln.raw_alignment), aln.raw_alignment.size());
//  ConvertInsertionsToClipping((unsigned char *) &(aln.raw_alignment[0]), aln.raw_alignment.size());
//  int64_t num_clipped_front = 0, num_clipped_back = 0;
//  CountClippedBases((unsigned char *) &(aln.raw_alignment[0]), aln.raw_alignment.size(), &num_clipped_front, &num_clipped_back);
//
//  /// This part converts the global alignment coordinates to local. 'Global' meaning the coordinates on the entire data array from the index.
//  /// 'Local' meaning the coordinates on the reference that was hit (in range [0, ref_len>).
//  /// final_aln_pos_start is the alignment position within the reference (local to the reference), and is also reversed (subtracted from reference_length) in case orientation == kReverse.
//  int64_t final_aln_pos_start = 0, final_aln_pos_end = 0;
//  index->RawPositionConverterWithRefId((orientation == kForward) ? aln_pos_start : aln_pos_end, abs_ref_id, 0, NULL, &final_aln_pos_start, NULL);
//  index->RawPositionConverterWithRefId((orientation == kForward) ? aln_pos_end : aln_pos_start, abs_ref_id, 0, NULL, &final_aln_pos_end, NULL);
//
//  /// Assign resulting values.
//  // aln.raw_alignment is assigned directly in the AlignmentFunction call to avoid copying of the data twice.
//  // aln.edit_distance is assigned directly in the AlignmentFunction call.
////  aln.ref_data = ref_data;
////  aln.is_cleanup_required = is_cleanup_required;
//  aln.orientation = orientation;
//  aln.is_reverse = (orientation == kReverse) ? true : false;
//  aln.ref_start = final_aln_pos_start;
//  aln.ref_end = final_aln_pos_end; // aln.pos_start + (aln_pos_end - aln_pos_start);
//  aln.query_start = (orientation == kForward) ? num_clipped_front : num_clipped_back;
//  aln.query_end = read->get_sequence_length() - ((orientation == kForward) ? num_clipped_back : num_clipped_front) - 1;
//  aln.raw_pos_start = aln_pos_start;
//  aln.raw_pos_end = aln_pos_end;
//  aln.ref_id = ref_id;
//  aln.ref_header = index->get_headers()[ref_id];
//  aln.ref_len = index->get_reference_lengths()[ref_id];
//  aln.query_id = read->get_sequence_absolute_id();
//  aln.query_header = read->get_header();
//  aln.query_len = read->get_sequence_length();
//  aln.is_aligned = true;
//  aln.aln_window_start = (ref_start + offset_from_ref_start);
//  aln.aln_window_end = aln.aln_window_start + ((int64_t) (win_end - win_start));
//  aln.aln_mode_code = MYERS_MODE_HW;
//  aln.alignment = aln.raw_alignment;
//
//  /// If alignment is linear, just add it. Otherwise, it needs to be split first.
//  if (region.is_split == false || parameters->is_reference_circular == false) {
//    if (orientation == kReverse) ReverseArray(aln.alignment);
//    region_results->AddAlignmentData(aln);
//  } else {
//    LOG_DEBUG_SPEC("Alignment is circular and will be split.\n");
////    printf ("Circular!\n");
////    exit(1);
//
//    AlignmentResults aln_r;
//    aln_r = aln;
//
//    int64_t end_on_aln = 0, end_on_read = 0, end_on_ref = 0;
//    int64_t start_on_aln = 0, start_on_read = 0, start_on_ref = 0;
//    FindCircularEnd(aln.raw_alignment, pos_of_ref_end,
//                    &end_on_aln, &end_on_read, &end_on_ref,
//                    &start_on_aln, &start_on_read, &start_on_ref);
//
////    aln.ref_end = end_on_ref;
////    aln.query_end = end_on_read;
//    aln.ref_end = aln.ref_start + end_on_ref;
//    aln.raw_pos_end = aln.raw_pos_start + end_on_ref;
//    uint8_t *alignment_l = (uint8_t *) malloc(sizeof(uint8_t) * (end_on_aln + (read->get_sequence_length() - (end_on_read + 1)) + 1));
//    memmove(&alignment_l[0], &aln.raw_alignment[0], end_on_aln + 1);
//    memset(&alignment_l[end_on_aln+1], EDLIB_S, (read->get_sequence_length() - (end_on_read + 1)  ));
//    LOG_DEBUG_SPEC("end_on_read = %ld\n", end_on_read);
//    LOG_DEBUG_SPEC("Read trim:\n%s\n", GetSubstring((char *) (read->get_data() + ((end_on_read + 1))), read->get_data_length() - (end_on_read + 1)).c_str());
//
//    aln_r.ref_start = start_on_ref;
//    aln_r.query_start = start_on_read;
//    uint8_t *alignment_r = (uint8_t *) malloc(sizeof(uint8_t) * ((aln.raw_alignment.size() - start_on_aln) + (start_on_read) + 1));
//    memset(&alignment_r[0], EDLIB_S, (start_on_read));
//    memmove(&alignment_r[start_on_read], &aln.raw_alignment[start_on_aln], (aln.raw_alignment.size() - start_on_aln));
//
////    printf ("Before:\n");
////    for (int64_t i=0; i<aln.raw_alignment.size(); i++) {
////      printf ("%d", aln.raw_alignment[i]);
////    }
////    printf ("\n");
////    printf ("\n");
////    fflush(stdout);
//
//    aln.raw_alignment.assign(alignment_l, alignment_l + (end_on_aln + (read->get_sequence_length() - (end_on_read + 1))));
//    aln.alignment = aln.raw_alignment;
//    if (orientation == kReverse) ReverseArray(aln.alignment);
//    // ref_data points to the start of the region (either on the index->data or on the newly malloced array). win_start points to the l1_start location which is located in ref_data.
//    // offset_from_ref_start is the distance from ref_data to the start of the reference. ref_start is the position of the beginning of the reference within index->data
//    aln.aln_window_start = (win_start - ref_data) + offset_from_ref_start + ref_start;
//    aln.aln_window_end = aln.aln_window_start + ((int64_t) (end_on_ref));
//    LOG_DEBUG_SPEC("aln.aln_window_start = %ld\n", aln.aln_window_start);
//    LOG_DEBUG_SPEC("aln.aln_window_end = %ld\n", aln.aln_window_end);
//
//    aln_r.raw_alignment.assign(alignment_r, alignment_r + ((aln.raw_alignment.size() - start_on_aln) + (start_on_read)));
//    aln_r.alignment = aln_r.raw_alignment;
//    if (orientation == kReverse) ReverseArray(aln_r.alignment);
//    aln_r.aln_window_start = start_on_ref + ref_start;
//    aln_r.aln_window_end = aln_r.aln_window_start + ((int64_t) (win_end - win_start - start_on_ref)) - 1;
//    LOG_DEBUG_SPEC("aln_r.aln_window_start = %ld\n", aln_r.aln_window_start);
//    LOG_DEBUG_SPEC("aln_r.aln_window_end = %ld\n", aln_r.aln_window_end);
//
//    printf ("After:\n");
//    for (int64_t i=0; i<aln.raw_alignment.size(); i++) {
//      printf ("%d", aln.raw_alignment[i]);
//    }
//    printf ("\n");
//    printf ("\n");
//    fflush(stdout);
//    LOG_DEBUG_SPEC("read:\n%s\n", GetSubstring((char *) (read->get_data()), 15).c_str());
//    LOG_DEBUG_SPEC("Ref:\n%s\n", GetSubstring((char *) (index->get_data() + aln.aln_window_start), 15).c_str());
//    LOG_DEBUG_SPEC("ref_data:\n%s\n", GetSubstring((char *) (win_start), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start:\n%s\n", GetSubstring((char *) (index->get_data() + aln_pos_start), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start2:\n%s\n", GetSubstring((char *) (win_start + aln_pos_start - (ref_start + offset_from_ref_start)), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start3:\n%s\n", GetSubstring((char *) (index->get_data() + ref_start + offset_from_ref_start + (win_start - ref_data) + aln_pos_start - (ref_start + offset_from_ref_start)), 15).c_str());
//    LOG_DEBUG_SPEC("aln_pos_start4:\n%s\n", GetSubstring((char *) (index->get_data() + aln_pos_start + (win_start - ref_data)), 15).c_str());
////    exit(1);
//
//    region_results->AddAlignmentData(aln);
//    region_results->AddAlignmentData(aln_r);
//  }
//
//  /// Fill out statistics for each alignment (E-value calculation, couting of CIGAR operations, etc.) and check if the alignments are sane.
//  for (int32_t i=0; i<region_results->get_alignments().size(); i++) {
//    AlignmentResults *curr_aln = &region_results->get_alignments()[i];
//
//    /// Verbose debug.
//    LOG_DEBUG_SPEC("Alignment part %d / %d:\n", i, region_results->get_alignments().size());
//    VerboseAlignment(read, index, parameters, curr_aln);
//
//    LOG_DEBUG_SPEC("Converting alignment to CIGAR string.\n");
//    curr_aln->cigar = AlignmentToCigar((unsigned char *) &(curr_aln->alignment[0]), curr_aln->alignment.size(), parameters->use_extended_cigar);
//
//    LOG_DEBUG_SPEC("Converting alignment to MD string.\n");
//    curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) curr_aln->alignment, read->get_data(), index->get_data(), ref_id, curr_aln->ref_start + index->get_reference_starting_pos()[curr_aln->ref_id]);
//
//    LOG_DEBUG_SPEC("Counting alignment operations.\n");
//    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), index->get_data(), ref_id, aln_pos_start, orientation,
//                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend,
//                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->nonclipped_length);
//
//    LOG_DEBUG_SPEC("Calculating the E-value.\n");
//    CalculateEValueDNA(curr_aln->alignment_score, curr_aln->nonclipped_length, index->get_data_length_forward(), evalue_params, &curr_aln->evalue);
//
//    LOG_DEBUG_SPEC("Checking if the alignment is sane.\n");
//    if (CheckAlignmentSane((std::vector<unsigned char> &) curr_aln->raw_alignment, read, index, curr_aln->ref_id, curr_aln->ref_start) != 0) {
//      curr_aln->is_aligned = false;
//    }
//  }
//
//  if (is_cleanup_required) { free(ref_data); }
//
//  /// Count the number of 'unaligned' alignments.
//  int32_t num_unaligned = 0;
//  for (int32_t i=0; i<region_results->get_alignments().size(); i++) { if (region_results->get_alignments()[i].is_aligned == false) num_unaligned += 1; }
//  if (num_unaligned == region_results->get_alignments().size()) { return ALIGNMENT_NOT_SANE; }
//
//  return 0;
//}

//int SemiglobalAlignment(AlignmentFunctionType AlignmentFunction, const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, PathGraphEntry *region_results) {
//  Region &region = region_results->get_region_data();
//
//  SeqOrientation orientation = (region_results->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
//  int64_t abs_ref_id = region_results->get_region_data().reference_id;
//  int64_t ref_id = (orientation == kForward) ? (abs_ref_id) : (abs_ref_id - index->get_num_sequences_forward());
//
//  int8_t *ref_data = NULL;
//  int64_t region_length_joined = 0, start_offset = 0, position_of_ref_end = 0;
//  if (parameters->is_reference_circular == false || region.is_split == false) {
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Linear alignment will be applied.\n", std::string(__FUNCTION__));
//    ref_data = (int8_t *) (index->get_data());
//  } else {
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Concatenating regions for circular alignment.\n", std::string(__FUNCTION__));
//    ConcatenateSplitRegion(index, region, &ref_data, &region_length_joined, &start_offset, &position_of_ref_end);
//  }
//
//  /// Calculate te offset for the L1 region filtering. This function considers whether the region was linear or circular.
//  int64_t l1_ref_start = 0, l1_ref_end = 0, reference_data_length = 0;
//  GetAlignmentWindowFromRegion(read, index, parameters, region_results, &l1_ref_start, &l1_ref_end, &reference_data_length);
//
//  /// Perform alignment and store results.
//  AlignmentResults aln;
//  std::vector<unsigned char> alignment1;
//  int64_t aln_pos_start = 0, aln_pos_end = 0;
//  int ret_code = AlignmentFunction(read->get_data(), read->get_sequence_length(),
//                                   ref_data + l1_ref_start, reference_data_length,
//                                   -1, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
//                                   &aln_pos_start, &aln_pos_end,
//                                   &aln.edit_distance,
//                                   aln.raw_alignment);
//  /// Sanity check.
//  if (ret_code != 0 || aln.raw_alignment.size() == 0) { return ret_code; }
//
//  aln_pos_start += l1_ref_start;
//  aln_pos_end += l1_ref_start;
//
////  aln.raw_alignment = FixAlignment((unsigned char *) &(aln.raw_alignment), aln.raw_alignment.size());
//  ConvertInsertionsToClipping((unsigned char *) &(aln.raw_alignment[0]), aln.raw_alignment.size());
//  int64_t num_clipped_front = 0, num_clipped_back = 0;
//  CountClippedBases((unsigned char *) &(aln.raw_alignment[0]), aln.raw_alignment.size(), &num_clipped_front, &num_clipped_back);
//
//  /// This part converts the global alignment coordinates to local. 'Global' meaning the coordinates on the entire data array from the index.
//  /// 'Local' meaning the coordinates on the reference that was hit (in range [0, ref_len>).
//  /// final_aln_pos_start is the alignment position within the reference (local to the reference), and is also reversed (subtracted from reference_length) in case orientation == kReverse.
//  int64_t final_aln_pos_start = 0, final_aln_pos_end = 0;
//  index->RawPositionConverterWithRefId((orientation == kForward) ? aln_pos_start : aln_pos_end, abs_ref_id, 0, NULL, &final_aln_pos_start, NULL);
//  index->RawPositionConverterWithRefId((orientation == kForward) ? aln_pos_end : aln_pos_start, abs_ref_id, 0, NULL, &final_aln_pos_end, NULL);
//
//  /// Assign resulting values.
//  // aln.raw_alignment is assigned directly in the AlignmentFunction call to avoid copying of the data twice.
//  // aln.edit_distance is assigned directly in the AlignmentFunction call.
//  aln.orientation = orientation;
//  aln.is_reverse = (orientation == kReverse) ? true : false;
//  aln.ref_start = final_aln_pos_start;
//  aln.ref_end = final_aln_pos_end; // aln.pos_start + (aln_pos_end - aln_pos_start);
//  aln.query_start = (orientation == kForward) ? num_clipped_front : num_clipped_back;
//  aln.query_end = read->get_sequence_length() - ((orientation == kForward) ? num_clipped_back : num_clipped_front) - 1;
//  aln.raw_pos_start = aln_pos_start;
//  aln.raw_pos_end = aln_pos_end;
//  aln.ref_id = ref_id;
//  aln.ref_header = index->get_headers()[ref_id];
//  aln.ref_len = index->get_reference_lengths()[ref_id];
//  aln.query_id = read->get_sequence_absolute_id();
//  aln.query_header = read->get_header();
//  aln.query_len = read->get_sequence_length();
//  aln.is_aligned = true;
//  aln.aln_window_start = l1_ref_start;
//  aln.aln_window_end = l1_ref_end;
//  aln.aln_mode_code = MYERS_MODE_HW;
//  aln.alignment = aln.raw_alignment;
//  if (orientation == kReverse) ReverseArray(aln.alignment);
//
//  region_results->AddAlignmentData(aln);
//
//  /// Verbose debug.
//  VerboseAlignment(read, index, parameters, &aln);
//
//  /// Fill out statistics for each alignment (E-value calculation, couting of CIGAR operations, etc.) and check if the alignments are sane.
//  for (int32_t i=0; i<region_results->get_alignments().size(); i++) {
//    AlignmentResults *curr_aln = &region_results->get_alignments()[i];
//    curr_aln->cigar = AlignmentToCigar((unsigned char *) &(curr_aln->alignment[0]), curr_aln->alignment.size(), parameters->use_extended_cigar);
//    curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) curr_aln->alignment, read->get_data(), index->get_data(), ref_id, curr_aln->ref_start + index->get_reference_starting_pos()[curr_aln->ref_id]);
//    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), index->get_data(), ref_id, aln_pos_start, orientation,
//                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend,
//                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->nonclipped_length);
//    CalculateEValueDNA(curr_aln->alignment_score, curr_aln->nonclipped_length, index->get_data_length_forward(), evalue_params, &curr_aln->evalue);
//    if (CheckAlignmentSane((std::vector<unsigned char> &) curr_aln->raw_alignment, read, index, curr_aln->ref_id, curr_aln->ref_start) != 0) {
//      curr_aln->is_aligned = false;
//    }
//  }
//
//  /// Count the number of 'unaligned' alignments.
//  int32_t num_unaligned = 0;
//  for (int32_t i=0; i<region_results->get_alignments().size(); i++) { if (region_results->get_alignments()[i].is_aligned == false) num_unaligned += 1; }
//  if (num_unaligned == region_results->get_alignments().size()) { return ALIGNMENT_NOT_SANE; }
//
//  return 0;
//}


