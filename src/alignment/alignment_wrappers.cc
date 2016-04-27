/*
 * local_realignment_generic.cc
 *
 *  Created on: Jan 3, 2015
 *      Author: isovic
 */

#include "alignment_wrappers.h"



int LocalizeAlignmentPosWithMyers(const int8_t *read_data, int64_t read_length,
                                  const int8_t *reference_data, int64_t reference_length,
                                  int64_t rough_reference_start, int64_t rough_reference_end,
                                  int64_t *ret_alignment_start, int64_t *ret_alignment_end,
                                  int64_t *ret_start_ambiguity, int64_t *ret_end_ambiguity,
                                  int64_t *ret_edit_distance, int64_t *ret_band_width,
                                  bool verbose_debug_output) {

  // Things needed for the Myers alignment functions.
  int current_score, current_alignment_length = 0, current_num_positions = 0;
  int *current_positions = NULL;
  unsigned char* current_alignment = NULL;
  std::string current_cigar = "";

  int64_t alignment_start = 0;
  int64_t alignment_end = 0;




  const unsigned char *query = (const unsigned char *) (read_data);
  int query_length = read_length;
  const unsigned char *target = (const unsigned char *) (reference_data + rough_reference_start);
  int target_length = rough_reference_end - rough_reference_start + 1;
  // Calculate the edit distance without the traceback using the HW alignment algorithm.
//  int myers_return_code = myersCalcEditDistance(query, query_length, target, target_length,
//                                                128, -1, EDLIB_MODE_HW, &current_score, &current_positions, &current_num_positions,
//                                                false, &current_alignment, &current_alignment_length);
  int *start_locations = NULL;

  int myers_return_code = edlibCalcEditDistance(query, query_length, target, target_length,
                        128, -1, EDLIB_MODE_HW, false, false, &current_score, &current_positions, &start_locations, &current_num_positions,
                        &current_alignment, &current_alignment_length);

  if (current_num_positions == 0 || myers_return_code != EDLIB_STATUS_OK) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, true, "Something went wrong when calculating the ending position using Myers HW. No positions were returned.\n", "CalculateAlignmentStartAndEnd");
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }
  alignment_end = current_positions[0];
  if (ret_end_ambiguity != NULL)
    *ret_end_ambiguity = current_num_positions;

  if (verbose_debug_output) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, "Alignment ending positions:\n", "CalculateAlignmentStartAndEnd");
    for (int64_t i=0; i<current_num_positions; i++) {
      if (i > 0)
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString(" ", current_positions[i]), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString("%ld", current_positions[i]), "[]");
    }
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, "\n", "[]");
  }

  if (current_positions) { free (current_positions); }  current_positions = NULL;
  if (start_locations) { free(start_locations); } start_locations = NULL;
  if (current_alignment) { free(current_alignment); } current_alignment = NULL;





  const unsigned char* reverse_query  = CreateReverseCopy(query, query_length);
  const unsigned char* reverse_target = CreateReverseCopy(target, (alignment_end + 1));

  int current_band_width = 0;

  myers_return_code = edlibCalcEditDistance(reverse_query, query_length, reverse_target, (alignment_end + 1),
                        128, -1, EDLIB_MODE_SHW, false, false, &current_score, &current_positions, &start_locations, &current_num_positions,
                        &current_alignment, &current_alignment_length, &current_band_width);

  if (current_num_positions == 0 || myers_return_code != EDLIB_STATUS_OK) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, "Something went wrong when calculating the starting position using Myers SHW. No positions were returned.\n", "CalculateAlignmentStartAndEnd");
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }
  alignment_start = alignment_end - current_positions[0];

  if (ret_start_ambiguity != NULL)
    *ret_start_ambiguity = current_num_positions;

  if (verbose_debug_output) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, "Alignment starting positions:\n", "CalculateAlignmentStartAndEnd");
    for (int64_t i=0; i<current_num_positions; i++) {
      if (i > 0)
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString(" ", current_positions[i]), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString("%ld", current_positions[i]), "[]");
    }
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, "\n", "[]");
  }

  if (current_positions) { free (current_positions); }  current_positions = NULL;
  if (start_locations) { free(start_locations); } start_locations = NULL;
  if (current_alignment) { free(current_alignment); } current_alignment = NULL;
  if (reverse_target) { delete[] reverse_target; } reverse_target = NULL;
  if (reverse_query) { delete[] reverse_query; } reverse_query = NULL;

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString("Alignment starting position: %ld\n", alignment_start), "CalculateAlignmentStartAndEnd");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString("Alignment ending position: %ld\n", alignment_end), "CalculateAlignmentStartAndEnd");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString("Alignment starting ambiguity: %ld\n", (*ret_start_ambiguity)), "CalculateAlignmentStartAndEnd");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString("Alignment ending ambiguity: %ld\n", (*ret_end_ambiguity)), "CalculateAlignmentStartAndEnd");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, verbose_debug_output, FormatString("Reference length: %ld\n", reference_length), "CalculateAlignmentStartAndEnd");

  *ret_alignment_start = alignment_start + rough_reference_start;
  *ret_alignment_end = alignment_end + rough_reference_start;
  *ret_edit_distance = current_score;
  *ret_band_width = current_band_width;

  return ALIGNMENT_GOOD;
}

int SeqAnAlignmentToEdlibAlignmentNoCigar(seqan::Align<seqan::Dna5String> &align, int alignment_type, int64_t *ret_start_offset, int64_t *ret_end_offset, int64_t *edit_distance, std::vector<unsigned char> &ret_alignment) {
  auto it_seq0 = seqan::begin(row(align, 0));
  auto it_seq0_end = seqan::end(row(align, 0));
  auto it_seq1 = seqan::begin(row(align, 1));
  auto it_seq1_end = seqan::end(row(align, 1));

  int64_t aln_seq0_len = seqan::length(row(align, 0));
  int64_t aln_seq1_len = seqan::length(row(align, 1));

  if (aln_seq0_len != aln_seq1_len) {
    return ALIGNMENT_CONVERSION_PROBLEM;
  }

  ret_alignment.clear();
  ret_alignment.reserve(aln_seq0_len);

  int64_t aln_ops = 0;
  while (it_seq0 != it_seq0_end && it_seq1 != it_seq1_end) {
    if (isGap(it_seq0) && (*it_seq0) != (*it_seq1)) {
      ret_alignment.push_back(EDLIB_I);
    } else if (isGap(it_seq1) && (*it_seq0) != (*it_seq1)) {
      ret_alignment.push_back(EDLIB_D);
    } else if (isGap(it_seq0) && isGap(it_seq1)) {
    } else {
      ret_alignment.push_back(((*it_seq0) == (*it_seq1)) ? EDLIB_EQUAL : EDLIB_X);
    }

    it_seq0++;
    it_seq1++;
  }

  // Find the offset at the beginning/end of the reference sequence through the number of deletion events.
  int64_t aln_len = ret_alignment.size();
  int64_t num_leading_dels = 0;
  for (int64_t i=0; i<aln_len; i++) {
    if (ret_alignment[i] != EDLIB_D) { break; }
    num_leading_dels += 1;
  }
  int64_t num_trailing_dels = 0;
  for (int64_t i=(aln_len-1); i>=0; i--) {
    if (ret_alignment[i] != EDLIB_D) { break; }
    num_trailing_dels += 1;
  }

  *ret_start_offset = num_leading_dels;
  *ret_end_offset = num_trailing_dels;

  if (alignment_type == ALIGNMENT_TYPE_NW) {
    *ret_start_offset = 0;
    *ret_end_offset = 0;
  }
  else if (alignment_type == ALIGNMENT_TYPE_SHW) {
    *ret_start_offset = 0;
    if (num_trailing_dels > 0) { ret_alignment.resize(ret_alignment.size() - num_trailing_dels); }
  }
  else if (alignment_type == ALIGNMENT_TYPE_HW) {
    /// If the alignment was not global, there might be trailing deletions
    /// if the reference has an overhang (which is most likely to happen).
    if (num_leading_dels > 0) { ret_alignment.erase(ret_alignment.begin(), (ret_alignment.begin() + num_leading_dels)); }
    if (num_trailing_dels > 0) { ret_alignment.resize(ret_alignment.size() - num_trailing_dels); }
  }

  int64_t aln_edit_distance = 0;
  for (int64_t i=0; i<ret_alignment.size(); i++) {
    aln_edit_distance += (ret_alignment[i] == EDLIB_EQUAL) ? 0 : 1;
  }

  *edit_distance = aln_edit_distance;

  return ALIGNMENT_GOOD;
}

int SeqAnSemiglobalWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Read length: %ld\n", read->get_sequence_length()), "SeqAnLocalRealignment");

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;



  seqan::Infix<char *>::Type inf_target = seqan::infix((char *) reference_data, 0, reference_length);
  seqan::Dna5String seq_target = inf_target;
  seqan::Infix<char *>::Type inf_query = seqan::infix((char *) read_data, 0, read_length);
  seqan::Dna5String seq_query = inf_query;

  seqan::Align<seqan::Dna5String> align;
  seqan::resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq_target);
  seqan::assignSource(row(align, 1), seq_query);

  int result = -1;
  if (band_width > 0) {
    globalAlignment(align,
                     seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                     seqan::AlignConfig<true, false, false, true>(),  // top, left, right, bottom
                     -band_width, band_width);
  } else {
    globalAlignment(align,
                     seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                     seqan::AlignConfig<true, false, false, true>()); // top, left, right, bottom
  }

  int64_t start_offset = 0, end_offset = 0, seqan_edit_distance = 0;
  if (SeqAnAlignmentToEdlibAlignmentNoCigar(align, ALIGNMENT_TYPE_HW, &start_offset, &end_offset, &seqan_edit_distance, ret_alignment) != 0)
    return ALIGNMENT_CONVERSION_PROBLEM;
  if (CheckAlignmentSaneSimple(ret_alignment))
    return ALIGNMENT_NOT_SANE;


  int64_t reconstructed_length = CalculateReconstructedLength((unsigned char *) &ret_alignment[0], ret_alignment.size());

  *ret_alignment_position_start = start_offset;
  *ret_alignment_position_end = start_offset + (reconstructed_length - 1);
  *ret_edit_distance = (int64_t) seqan_edit_distance;

  return ALIGNMENT_GOOD;
}

int SeqAnNWWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;



  seqan::Infix<char *>::Type inf_target = seqan::infix((char *) reference_data, 0, reference_length);
  seqan::Dna5String seq_target = inf_target;
  seqan::Infix<char *>::Type inf_query = seqan::infix((char *) read_data, 0, read_length);
  seqan::Dna5String seq_query = inf_query;

  seqan::Align<seqan::Dna5String> align;
  seqan::resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq_target);
  seqan::assignSource(row(align, 1), seq_query);

  int result = -1;
  if (band_width > 0) {
    globalAlignment(align,
                     seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                     seqan::AlignConfig<false, false, false, false>(),  // top, left, right, bottom
                     -band_width, band_width);
  } else {
    globalAlignment(align,
                     seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                     seqan::AlignConfig<false, false, false, false>()); // top, left, right, bottom
  }

  int64_t start_offset = 0, end_offset = 0, seqan_edit_distance = 0;
  if (SeqAnAlignmentToEdlibAlignmentNoCigar(align, ALIGNMENT_TYPE_NW, &start_offset, &end_offset, &seqan_edit_distance, ret_alignment) != 0)
    return ALIGNMENT_CONVERSION_PROBLEM;
  if (CheckAlignmentSaneSimple(ret_alignment))
    return ALIGNMENT_NOT_SANE;

  int64_t reconstructed_length = CalculateReconstructedLength((unsigned char *) &ret_alignment[0], ret_alignment.size());

  *ret_alignment_position_start = start_offset;
  *ret_alignment_position_end = start_offset + (reconstructed_length - 1);
  *ret_edit_distance = (int64_t) seqan_edit_distance;

  return ALIGNMENT_GOOD;
}

int SeqAnSHWWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  seqan::Infix<char *>::Type inf_target = seqan::infix((char *) reference_data, 0, reference_length);
  seqan::Dna5String seq_target = inf_target;
  seqan::Infix<char *>::Type inf_query = seqan::infix((char *) read_data, 0, read_length);
  seqan::Dna5String seq_query = inf_query;

  if (read_length <= reference_length) {
    seqan::Align<seqan::Dna5String> align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq_target);
    seqan::assignSource(row(align, 1), seq_query);

    int result = -1;
    if (band_width > 0) {
      globalAlignment(align,
                       seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                       seqan::AlignConfig<false, false, false, true>(),  // top, left, right, bottom
                       -band_width, band_width);
    } else {
      globalAlignment(align,
                       seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                       seqan::AlignConfig<false, false, false, true>()); // top, left, right, bottom
    }

    int64_t start_offset = 0, end_offset = 0, seqan_edit_distance = 0;
    if (SeqAnAlignmentToEdlibAlignmentNoCigar(align, ALIGNMENT_TYPE_SHW, &start_offset, &end_offset, &seqan_edit_distance, ret_alignment) != 0)
      return ALIGNMENT_CONVERSION_PROBLEM;
    if (CheckAlignmentSaneSimple(ret_alignment))
      return ALIGNMENT_NOT_SANE;

    int64_t reconstructed_length = CalculateReconstructedLength((unsigned char *) &ret_alignment[0], ret_alignment.size());

    *ret_alignment_position_start = start_offset;
//    printf("start_offset = %ld\n", start_offset);
//    fflush(stdout);

    *ret_alignment_position_end = start_offset + (reconstructed_length - 1);
    *ret_edit_distance = (int64_t) seqan_edit_distance;

  } else {
    seqan::Align<seqan::Dna5String> align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq_query);
    seqan::assignSource(row(align, 1), seq_target);

    int result = -1;
    if (band_width > 0) {
      globalAlignment(align,
                       seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                       seqan::AlignConfig<false, false, false, true>(),  // top, left, right, bottom
                       -band_width, band_width);
    } else {
      globalAlignment(align,
                       seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                       seqan::AlignConfig<false, false, false, true>()); // top, left, right, bottom
    }

    std::vector<unsigned char> alignment;

    int64_t start_offset = 0, end_offset = 0, seqan_edit_distance = 0;
    if (SeqAnAlignmentToEdlibAlignmentNoCigar(align, ALIGNMENT_TYPE_SHW, &start_offset, &end_offset, &seqan_edit_distance, alignment) != 0)
      return ALIGNMENT_CONVERSION_PROBLEM;
    if (CheckAlignmentSaneSimple(alignment))
      return ALIGNMENT_NOT_SANE;

    int64_t reconstructed_length = CalculateReconstructedLength((unsigned char *) &alignment[0], alignment.size());

    ret_alignment.clear();
    int64_t start_on_read = start_offset; // - (reconstructed_length - 1);
    int64_t end_on_read = start_offset + reconstructed_length - 1;

    // Check if the beginning of the alignment on the read is not at 0. This means we need to insert insertions at the beginning.
    if (start_on_read > 0) {
      std::vector<unsigned char> leading_insertions(start_on_read, EDLIB_I);
      ret_alignment.insert(ret_alignment.end(), leading_insertions.begin(), leading_insertions.end());
    }
    // Convert insertions to deletions and vice versa, becase we are changing contexts back to the original (reference-read).
    for (int64_t i=0; i<alignment.size(); i++) {
      if (alignment[i] == EDLIB_I) { alignment[i] = EDLIB_D; }
      else if  (alignment[i] == EDLIB_D) { alignment[i] = EDLIB_I; }
    }
    // Check if there are leading/trailing deletions. These should be removed.
    int64_t leading_del_offset = 0;
    for (leading_del_offset=0; leading_del_offset<alignment.size(); leading_del_offset++) { if (alignment[leading_del_offset] != 'D') break; }
    int64_t trailing_del_offset = (alignment.size()-1);
    for (trailing_del_offset=(alignment.size()-1); trailing_del_offset>=0; trailing_del_offset--) { if (alignment[trailing_del_offset] != 'D') break; }

    ret_alignment.insert(ret_alignment.end(), alignment.begin()+leading_del_offset, (alignment.begin() + trailing_del_offset + 1));
    // Check if there are any bases left unaligned on the read. If so, insert insertions at the end of the alignment.
    if (end_on_read < read_length) {
      std::vector<unsigned char> trailing_insertions((read_length - end_on_read - 1), EDLIB_I);
      ret_alignment.insert(ret_alignment.end(), trailing_insertions.begin(), trailing_insertions.end());
    }

    *ret_alignment_position_start = 0; // positions[0] - (reconstructed_length - 1);
    *ret_alignment_position_end = reference_length - 1; // positions[0];
    *ret_edit_distance = (int64_t) seqan_edit_distance;
  }

  return ALIGNMENT_GOOD;
}

int SeqAnSemiglobalWrapperWithMyersLocalization(const int8_t *read_data, int64_t read_length,
                                                const int8_t *reference_data, int64_t reference_length,
                                                int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                                                int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                                                int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0) {
    return ALIGNMENT_WRONG_DATA;
  }

  // Find the start and end positions of the optimal alignment with Myers bit-vector algorithm. The Myers' algorithm uses all parameters equal to 1.
  int64_t localized_start = 0, localized_end = 0, ambiguity_start = 0, ambiguity_end = 0, localized_edit_distance = 0, localized_band_width = 0;
  int ret_code = LocalizeAlignmentPosWithMyers(read_data, read_length,
                                               reference_data, reference_length,
                                               0, (reference_length - 1),
                                               &localized_start, &localized_end,
                                               &ambiguity_start, &ambiguity_end,
                                               &localized_edit_distance, &localized_band_width, false);
  if (ret_code != 0)
    return ALIGNMENT_LOCALIZATION_PROBLEM;

  band_width = localized_band_width;

  // Expand the search field a bit to allow for slight jiggling of the optimal alignment positions. The differences are possible because of different alignment parameters.
  localized_start -= 50;
  if (localized_start < 0)
    localized_start = 0;
  localized_end += 50;
  if (localized_end >= reference_length)
    localized_end = reference_length - 1;

  seqan::Infix<char *>::Type inf_target = seqan::infix((char *) (reference_data + localized_start), 0, (localized_end - localized_start));
  seqan::Dna5String seq_target = inf_target;
  seqan::Infix<char *>::Type inf_query = seqan::infix((char *) read_data, 0, read_length);
  seqan::Dna5String seq_query = inf_query;

  seqan::Align<seqan::Dna5String> align;
  seqan::resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq_target);
  seqan::assignSource(row(align, 1), seq_query);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("band_width = %ld\n", band_width), "SeqAnLocalRealignment");

  int result = -1;
//  band_width = 0;
//  band_width *= 2;
  if (band_width > 0) {
    globalAlignment(align,
                     seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                     seqan::AlignConfig<true, false, false, true>(),  // top, left, right, bottom
                     -band_width, band_width);
  } else {
    globalAlignment(align,
                     seqan::Score<int, seqan::Simple>(match_score, mismatch_penalty, gap_extend_penalty, gap_open_penalty),
                     seqan::AlignConfig<true, false, false, true>()); // top, left, right, bottom
  }

  int64_t start_offset = 0, end_offset = 0, seqan_edit_distance = 0;
  if (SeqAnAlignmentToEdlibAlignmentNoCigar(align, ALIGNMENT_TYPE_HW, &start_offset, &end_offset, &seqan_edit_distance, ret_alignment) != 0)
    return ALIGNMENT_CONVERSION_PROBLEM;
  if (CheckAlignmentSaneSimple(ret_alignment))
    return ALIGNMENT_NOT_SANE;

  int64_t reconstructed_length = CalculateReconstructedLength((unsigned char *) &ret_alignment[0], ret_alignment.size());

  *ret_alignment_position_start = start_offset + localized_start;
  *ret_alignment_position_end = start_offset + localized_start + (reconstructed_length - 1);
  *ret_edit_distance = (int64_t) localized_edit_distance;

  return ALIGNMENT_GOOD;
}

int MyersSemiglobalWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  int alphabet_length = 128;
  int score = 0;
  unsigned char* alignment = NULL;
  int alignment_length = 0;

  int *positions = NULL;
  int num_positions = 0;
  int *start_locations = NULL;
  int found_k = 0;

  int myers_return_code = edlibCalcEditDistance((const unsigned char *) read_data, read_length,
                        (const unsigned char *) reference_data, reference_length,
                        alphabet_length, -1, EDLIB_MODE_HW, false, true, &score, &positions, &start_locations, &num_positions,
                        &alignment, &alignment_length, &found_k);

  if (myers_return_code == EDLIB_STATUS_ERROR || num_positions == 0 || alignment_length == 0) {
    if (positions)
      free(positions);
    if (start_locations)
      free(start_locations);
    if (alignment)
        free(alignment);
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }

  int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

  *ret_alignment_position_start = positions[0] - (reconstructed_length - 1);
  *ret_alignment_position_end = positions[0];
  *ret_edit_distance = (int64_t) score;
  ret_alignment.assign(alignment, (alignment + alignment_length));

  if (positions)
    free(positions);

  if (start_locations)
    free(start_locations);

  if (alignment)
      free(alignment);

  return ALIGNMENT_GOOD;
}

int MyersNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  int alphabet_length = 128;
  int score = 0;
  unsigned char* alignment = NULL;
  int alignment_length = 0;

  int *positions = NULL;
  int num_positions = 0;
  int *start_locations = NULL;
  int found_k = 0;

  int myers_return_code = edlibCalcEditDistance((const unsigned char *) read_data, read_length,
                        (const unsigned char *) reference_data, reference_length,
                        alphabet_length, -1, EDLIB_MODE_NW, false, true, &score, &positions, &start_locations, &num_positions,
                        &alignment, &alignment_length, &found_k);

  if (myers_return_code == EDLIB_STATUS_ERROR || num_positions == 0 || alignment_length == 0) {
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }

  int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

  *ret_alignment_position_start = positions[0] - (reconstructed_length - 1);
  *ret_alignment_position_end = positions[0];
  *ret_edit_distance = (int64_t) score;
  ret_alignment.assign(alignment, (alignment + alignment_length));

  if (positions)
    free(positions);

  if (start_locations)
    free(start_locations);

  if (alignment)
      free(alignment);

  return ALIGNMENT_GOOD;
}

#ifndef RELEASE_VERSION
int OpalNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

//  mismatch_penalty = abs(mismatch_penalty);
  gap_open_penalty = abs(gap_open_penalty);
  gap_extend_penalty = abs(gap_extend_penalty);

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  uint8_t *converted_data = new uint8_t[read_length];
  if (converted_data == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_data."));
    return ALIGNMENT_WRONG_DATA;
  }
  for (int64_t i=0; i<read_length; i++) {
    converted_data[i] = kBaseToBwaUnsigned[read_data[i]];
  }

  uint8_t *converted_ref = new uint8_t[reference_length];
  if (converted_ref == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_ref."));
    return ALIGNMENT_WRONG_DATA;
  }
  for (int64_t i=0; i<reference_length; i++) {
    converted_ref[i] = kBaseToBwaUnsigned[reference_data[i]];
  }

  int db_length = 1;
  uint8_t* db[1] = {converted_ref};
  int32_t db_seq_lengths[] = {((int32_t) reference_length)};
  OpalSearchResult* results[db_length];
  for (int i = 0; i < db_length; i++) {
      results[i] = new OpalSearchResult;
      opalInitSearchResult(results[i]);
  }

//  printf ("match_score = %ld\n", match_score);
//  printf ("mex_score = %ld\n", mex_score);
//  printf ("mismatch_penalty = %ld\n", mismatch_penalty);
//  printf ("gap_open_penalty = %ld\n", gap_open_penalty);
//  printf ("gap_extend_penalty = %ld\n", gap_extend_penalty);
//  printf ("Read:\n");
//  for (int64_t i=0; i<read_length; i++) {
//    printf ("%d", converted_data[i]);
//  }
//  printf ("\n\n");
//  printf ("Reference:\n");
//  for (int64_t i=0; i<reference_length; i++) {
//    printf ("%d", converted_ref[i]);
//  }
//  printf("\n");
//  fflush(stdout);

  int alphabet_length = 5;
  int scoreMatrix[25] = {
      (int) match_score,      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty,
      (int) mismatch_penalty, (int) match_score,      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty,
      (int) mismatch_penalty, (int) mismatch_penalty, (int) match_score,      (int) mismatch_penalty, (int) mismatch_penalty,
      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) match_score,      (int) mismatch_penalty,
      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty
  };

  int resultCode = resultCode = opalSearchDatabase(converted_data, read_length, db, db_length, db_seq_lengths,
                                                    gap_open_penalty, gap_extend_penalty, mex_score, scoreMatrix, alphabet_length, results,
                                                    OPAL_SEARCH_ALIGNMENT, OPAL_MODE_NW, OPAL_OVERFLOW_SIMPLE);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, "Finished running Opal.\n", "OpalNWWrapper");

  if (resultCode == OPAL_ERR_OVERFLOW) {
    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Opal returned with overflow error!"));
    return ALIGNMENT_OPAL_OVERFLOW_ERROR;
  }
  if (resultCode == OPAL_ERR_NO_SIMD_SUPPORT) {
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Opal returned with error, no SIMD support!"));
      return ALIGNMENT_OPAL_NO_SIMD;
  }

  if (results[0]->alignment == NULL) {
//    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Opal: no alignment was generated!"));
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Opal: no alignment was generated!"), std::string(__FUNCTION__));
    return ALIGNMENT_NOT_SANE;
  }

//  if (results[0]->endLocationQuery > 0) {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Opal: query end location is > 0! results[0]->endLocationQuery = %d, read_length = %ld\n", results[0]->endLocationQuery, read_length), "ProcessRead");
//  }
//
//  if (results[0]->endLocationQuery < read_length) {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Opal: query end location is < read_length! results[0]->endLocationQuery = %d, read_length = %ld\n", results[0]->endLocationQuery, read_length), "ProcessRead");
//  }

  *ret_alignment_position_start = results[0]->startLocationTarget;
  *ret_alignment_position_end = results[0]->endLocationTarget;
  *ret_edit_distance = (int64_t) results[0]->score;
  ret_alignment.assign(results[0]->alignment, (results[0]->alignment + results[0]->alignmentLength));

  /// Convert Opal alignment codes to EDLIB alignment codes.
//  const uint8_t alignment_lookup[] = {EDLIB_EQUAL, EDLIB_D, EDLIB_I, EDLIB_X};
//  for (int64_t i=0; i<ret_alignment.size(); i++) {
//    ret_alignment[i] = alignment_lookup[ret_alignment[i]];
//  }

//  printf ("Alignment:\n");
//  for (int64_t i=0; i<ret_alignment.size(); i++) {
//    printf ("%d", ret_alignment[i]);
//  }
//  printf ("\n\n");
//  fflush(stdout);

  if (converted_data)
    delete[] converted_data;
  if (converted_ref)
    delete[] converted_ref;

  return ALIGNMENT_GOOD;
}
#endif

#ifndef RELEASE_VERSION
int OpalSHWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  uint8_t *converted_data = new uint8_t[read_length];
  if (converted_data == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_data."));
    return 1;
  }
  for (int64_t i=0; i<read_length; i++) {
    converted_data[i] = kBaseToBwaUnsigned[read_data[i]];
  }

  uint8_t *converted_ref = new uint8_t[reference_length];
  if (converted_ref == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_ref."));
    return 1;
  }
  for (int64_t i=0; i<reference_length; i++) {
    converted_ref[i] = kBaseToBwaUnsigned[reference_data[i]];
  }

  int db_length = 1;
  uint8_t* db[1] = {converted_ref};
  int32_t db_seq_lengths[] = {((int32_t) reference_length)};
  OpalSearchResult* results[db_length];
  for (int i = 0; i < db_length; i++) {
      results[i] = new OpalSearchResult;
      opalInitSearchResult(results[i]);
  }

  int alphabet_length = 5;
  int scoreMatrix[25] = {
      (int) match_score,      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty,
      (int) mismatch_penalty, (int) match_score,      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty,
      (int) mismatch_penalty, (int) mismatch_penalty, (int) match_score,      (int) mismatch_penalty, (int) mismatch_penalty,
      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) match_score,      (int) mismatch_penalty,
      (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty, (int) mismatch_penalty
  };

  int resultCode = resultCode = opalSearchDatabase(converted_data, read_length, db, db_length, db_seq_lengths,
                                                    gap_open_penalty, gap_extend_penalty, mex_score, scoreMatrix, alphabet_length, results,
                                                    OPAL_SEARCH_ALIGNMENT, OPAL_MODE_NW, OPAL_OVERFLOW_BUCKETS);

  if (resultCode == OPAL_ERR_OVERFLOW) {
    LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Opal returned with overflow error!"));
    return 2;
  }
  if (resultCode == OPAL_ERR_NO_SIMD_SUPPORT) {
      LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Opal returned with error, no SIMD support!"));
      return 3;
  }

  if (results[0]->alignment == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Opal: no alignment was generated!"));
    return 4;
  }

  if (results[0]->startLocationQuery > 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, "Opal: query start location is > 0!", "ProcessRead");
  }

  if (results[0]->endLocationQuery > 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Opal: query end location is > 0! results[0]->endLocationQuery = %d", results[0]->endLocationQuery), "ProcessRead");
  }

  if (results[0]->endLocationQuery < read_length) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("Opal: query end location is < read_length! results[0]->endLocationQuery = %d", results[0]->endLocationQuery), "ProcessRead");
  }

  *ret_alignment_position_start = results[0]->startLocationTarget;
  *ret_alignment_position_end = results[0]->endLocationTarget;
  *ret_edit_distance = (int64_t) results[0]->score;
  ret_alignment.assign(results[0]->alignment, (results[0]->alignment + results[0]->alignmentLength));

  /// Convert Opal alignment codes to EDLIB alignment codes.
  const uint8_t alignment_lookup[] = {EDLIB_EQUAL, EDLIB_D, EDLIB_I, EDLIB_X};
  for (int64_t i=0; i<ret_alignment.size(); i++) {
    ret_alignment[i] = alignment_lookup[ret_alignment[i]];
  }

  if (converted_data)
    delete[] converted_data;
  if (converted_ref)
    delete[] converted_ref;

  return ALIGNMENT_GOOD;
}
#endif


int MyersSHWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mex_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  int alphabet_length = 128;
  int score = 0;
  unsigned char* alignment = NULL;
  int alignment_length = 0;
  int *positions = NULL;
  int num_positions = 0;
  int found_k = 0;
  int *start_locations = NULL;

  if (read_length <= reference_length) {
    int myers_return_code = edlibCalcEditDistance((const unsigned char *) read_data, read_length,
                                                  (const unsigned char *) reference_data, reference_length,
                                                  alphabet_length, -1, EDLIB_MODE_SHW, false, true, &score, &positions, &start_locations, &num_positions,
                                                  &alignment, &alignment_length, &found_k);

    if (myers_return_code == EDLIB_STATUS_ERROR) {
      return ALIGNMENT_MYERS_INTERNAL_ERROR;
    }
    if (num_positions == 0) {
      return ALIGNMENT_MYERS_INTERNAL_ERROR;
    }
    if (alignment_length == 0) {
      return ALIGNMENT_MYERS_INTERNAL_ERROR;
    }

    int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

    *ret_alignment_position_start = positions[0] - (reconstructed_length - 1);
    *ret_alignment_position_end = positions[0];
    *ret_edit_distance = (int64_t) score;
    ret_alignment.assign(alignment, (alignment + alignment_length));
  } else {
    int myers_return_code = edlibCalcEditDistance((const unsigned char *) reference_data, reference_length,
                                                  (const unsigned char *) read_data, read_length,
                                                  alphabet_length, -1, EDLIB_MODE_SHW, false, true, &score, &positions, &start_locations, &num_positions,
                                                  &alignment, &alignment_length, &found_k);

    if (myers_return_code == EDLIB_STATUS_ERROR || num_positions == 0 || alignment_length == 0) {
      if (positions)
        free(positions);
      if (alignment)
          free(alignment);
      return ALIGNMENT_MYERS_INTERNAL_ERROR;
    }

    int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

    ret_alignment.clear();
    int64_t start_on_read = positions[0] - (reconstructed_length - 1);
    // Check if the beginning of the alignment on the read is not at 0. This means we need to insert insertions at the beginning.
    if (start_on_read > 0) {
      std::vector<unsigned char> leading_insertions(start_on_read, EDLIB_I);
      ret_alignment.insert(ret_alignment.end(), leading_insertions.begin(), leading_insertions.end());
    }
    // Convert insertions to deletions and vice versa, becase we are changing contexts back to the original (reference-read).
    for (int64_t i=0; i<alignment_length; i++) {
      if (alignment[i] == EDLIB_I) { alignment[i] = EDLIB_D; }
      else if  (alignment[i] == EDLIB_D) { alignment[i] = EDLIB_I; }
    }
    // Check if there are leading/trailing deletions. These should be removed.
    int64_t leading_del_offset = 0;
    for (leading_del_offset=0; leading_del_offset<alignment_length; leading_del_offset++) { if (alignment[leading_del_offset] != 'D') break; }
    int64_t trailing_del_offset = (alignment_length-1);
    for (trailing_del_offset=(alignment_length-1); trailing_del_offset>=0; trailing_del_offset--) { if (alignment[trailing_del_offset] != 'D') break; }

    ret_alignment.insert(ret_alignment.end(), alignment+leading_del_offset, (alignment + trailing_del_offset + 1));
    // Check if there are any bases left unaligned on the read. If so, insert insertions at the end of the alignment.
    if (positions[0] < read_length) {
      std::vector<unsigned char> trailing_insertions((read_length - positions[0] - 1), EDLIB_I);
      ret_alignment.insert(ret_alignment.end(), trailing_insertions.begin(), trailing_insertions.end());
    }

    *ret_alignment_position_start = 0; // positions[0] - (reconstructed_length - 1);
    *ret_alignment_position_end = reference_length - 1; // positions[0];
    *ret_edit_distance = (int64_t) score;

  }

  if (positions)
    free(positions);

  if (alignment)
      free(alignment);

  return ALIGNMENT_GOOD;
}

int MyersEditDistanceWrapper(const int8_t *read_data, int64_t read_length,
                             const int8_t *reference_data, int64_t reference_length,
                             int64_t *ret_alignment_position_end,
                             int64_t *ret_edit_distance, int myers_mode_code) {

    if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
      return ALIGNMENT_WRONG_DATA;

    int alphabet_length = 128;
    int score = 0;
    unsigned char* alignment = NULL;
    int alignment_length = 0;

    int *positions = NULL;
    int num_positions = 0;
    int *start_locations = NULL;
    int found_k = 0;

    int myers_return_code = edlibCalcEditDistance((const unsigned char *) read_data, read_length,
                          (const unsigned char *) reference_data, reference_length,
                          alphabet_length, -1, myers_mode_code, false, false, &score, &positions, &start_locations, &num_positions,
                          &alignment, &alignment_length, &found_k);

    if (myers_return_code == EDLIB_STATUS_ERROR || num_positions == 0) {
      return ALIGNMENT_MYERS_INTERNAL_ERROR;
    }

    if (ret_alignment_position_end != NULL)
      *ret_alignment_position_end = positions[0];
    if (ret_edit_distance != NULL)
      *ret_edit_distance = (int64_t) score;

    if (positions)
      free(positions);

    if (start_locations)
      free(start_locations);

    if (alignment)
        free(alignment);

    return ALIGNMENT_GOOD;
}

// Checks if there is a strange (large) number of insertions and deletions, or consecutive insertion/deletion operations.
// SeqAn likes to make such alignments.
// Returns 0 if everything went ok.
int CheckAlignmentSaneSimple(std::vector<unsigned char> &alignment) {
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
////            printf ("CheckAlignmentSane returned false! return 1.\n");
////            fflush(stdout);
//            return 1;
//          }
          // If there are insertions following deletions (or other way around), something is wrong again.
          if ((last_move == EDLIB_I && alignment_char == EDLIB_D) || (last_move == EDLIB_D && alignment_char == EDLIB_I)) {
//            printf ("CheckAlignmentSane returned false! return 2.\n");
//            fflush(stdout);
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

  return 0;
}
