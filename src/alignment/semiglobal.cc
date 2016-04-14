/*
 * semiglobal.cc
 *
 *  Created on: Jan 27, 2016
 *      Author: isovic
 */

#include <alignment/alignment_wrappers.h>
#include "alignment/alignment.h"

int SemiglobalAlignment(AlignmentFunctionType AlignmentFunction,
                        const SingleSequence *read, const Index *index, const ProgramParameters *parameters,
                        const EValueParams *evalue_params, PathGraphEntry *region_results) {
  /// General useful things.
  const Region &region = region_results->get_region_data();
  SeqOrientation orientation = (region_results->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t abs_ref_id = region_results->get_region_data().reference_id;
  int64_t ref_id = (orientation == kForward) ? (abs_ref_id) : (abs_ref_id - index->get_num_sequences_forward());
  int64_t ref_start = index->get_reference_starting_pos()[region.reference_id];
  int64_t ref_len = index->get_reference_lengths()[region.reference_id];
  int64_t reference_length = index->get_reference_lengths()[abs_ref_id];

  int8_t *reg_data = NULL;       // The data of the region.
  int64_t reg_data_len = 0;         // Length of the region_data.
  int64_t index_pos = 0;            // Position of the start of the region on the original Index data. E.g. reg_data[0] is the same base as index->get_data()[index_pos].
  int64_t pos_of_ref_end = 0; // If the region was circular, it crosses the boundary between the end and the start of the data. ref_data[index_pos_of_ref_end] is the last base of the reference before the split part is concatenated.
  bool is_cleanup_required = NULL;  // If true, the region_data will have to be freed manually.
  int64_t l1_start = 0;             // The start position of the L1 determined boundaries, but within the region (relative to region start).
  int64_t l1_end = 0;               // The end position of the L1 determined boundaries, but within the region (relative to region start).

  int64_t region_ref_start = region.start;            // Position of the start of the region on the original Index data. E.g. reg_data[0] is the same base as index->get_data()[index_pos].

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
  aln.aln_mode_code = EDLIB_MODE_HW;
  aln.alignment = aln.raw_alignment;

//  VerboseAlignment(read, index, parameters, &aln);
  CountAlignmentOperations((std::vector<unsigned char> &) aln.raw_alignment, read->get_data(), reg_data, ref_id, aln.reg_pos_start, orientation,
                           parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend, true,
                           &aln.num_eq_ops, &aln.num_x_ops, &aln.num_i_ops, &aln.num_d_ops, &aln.alignment_score, &aln.edit_distance, &aln.nonclipped_length);
  LOG_DEBUG_SPEC("Calculating the E-value.\n");
  CalculateEValueDNA(aln.alignment_score, aln.nonclipped_length, index->get_data_length_forward(), evalue_params, &aln.evalue);

  /// If alignment is linear, just add it. Otherwise, it needs to be split first.
  if (region.is_split == false || parameters->is_reference_circular == false) {
    LOG_DEBUG_SPEC("Alignment is linear.\n");
    if (orientation == kReverse) ReverseArray(aln.alignment);
    region_results->AddAlignmentData(aln);

  } else {
    LOG_DEBUG_SPEC("Alignment is circular and will be split.\n");
    AlignmentResults aln_l;
    AlignmentResults aln_r;
    SplitCircularAlignment(&aln, pos_of_ref_end, ref_start, ref_len, &aln_l, &aln_r);

    // These need to be fixed. The rest of the above code thinks that the region is linear and not circular.
    // In case the region is linear, anchors/clusters will have absolute reference coordinates.
    // The absolute coordinates allow the reference data to be accessed seemlesly.
    // For the circular alignment, the data needed to be copied into a new array, which is indexed from zero.
    // That's why the ref coordinates of anchors/clusters do not correspond to their reference positions. The region
    // information still contains still contains the region's start position on the reference, and this one will be
    // subtracted from the ref coordinates of anchors/clusters when calculating aln.reg_pos_start and aln.reg_pos_end.
    // For this reason we need to increase it back to obtain correct coodinates.
//    aln.reg_pos_start += region_ref_start;
//    aln.reg_pos_end += region_ref_start;
//    aln.raw_pos_start += region_ref_start;
//    aln.raw_pos_end += region_ref_start;
//    aln.ref_start += region_ref_start;
//    aln.ref_end += region_ref_start;

//    SplitCircularAlignment(&aln, pos_of_ref_end, 0, reference_length, &aln_l, &aln_r);

    region_results->AddAlignmentData(aln_l);
    region_results->AddAlignmentData(aln_r);
  }

  /// Fill out statistics for each alignment (E-value calculation, couting of CIGAR operations, etc.) and check if the alignments are sane.
  for (int32_t i=0; i<region_results->get_alignments().size(); i++) {
    AlignmentResults *curr_aln = &region_results->get_alignments()[i];

    if (curr_aln->is_aligned == false) { continue; }

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
    curr_aln->ref_start = final_aln_pos_start; // % ref_len;
    curr_aln->ref_end = final_aln_pos_end; // % ref_len;

    LOG_DEBUG_SPEC("Converting alignment to CIGAR string.\n");
    curr_aln->cigar = AlignmentToCigar((unsigned char *) &(curr_aln->alignment[0]), curr_aln->alignment.size(), parameters->use_extended_cigar);

    LOG_DEBUG_SPEC("Converting alignment to MD string.\n");
    curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) curr_aln->alignment, read->get_data(), index->get_data(), ref_id, curr_aln->ref_start + ref_start + index->get_reference_starting_pos()[ref_id]);

    LOG_DEBUG_SPEC("Counting alignment operations.\n");
    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), reg_data, ref_id, curr_aln->reg_pos_start, orientation,
                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend, true,
                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->edit_distance, &curr_aln->nonclipped_length);
//    LOG_DEBUG_SPEC("Calculating the E-value.\n");
//    CalculateEValueDNA(curr_aln->alignment_score, curr_aln->nonclipped_length, index->get_data_length_forward(), evalue_params, &curr_aln->evalue);

//    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->alignment,
//                             read->get_data(), index->get_data() + curr_aln->ref_start + index->get_reference_starting_pos()[ref_id],
//                             ref_id, 0, orientation,
//                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend,
//                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->nonclipped_length);

    LOG_DEBUG_SPEC("Checking if the alignment is sane.\n");
    if (CheckAlignmentSane((std::vector<unsigned char> &) curr_aln->raw_alignment, read, index, curr_aln->ref_id, curr_aln->ref_start) != 0) {
      curr_aln->is_aligned = false;
      LOG_DEBUG_SPEC("Alignment is insane!\n");
    } else {
      LOG_DEBUG_SPEC("Alignment is ok!\n");
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
