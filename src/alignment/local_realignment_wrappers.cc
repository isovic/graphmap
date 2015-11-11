/*
 * local_realignment_generic.cc
 *
 *  Created on: Jan 3, 2015
 *      Author: isovic
 */

#include "alignment/local_realignment_wrappers.h"



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
  int myers_return_code = myersCalcEditDistance(query, query_length, target, target_length,
                                                128, -1, MYERS_MODE_HW, &current_score, &current_positions, &current_num_positions,
                                                false, &current_alignment, &current_alignment_length);
  if (current_num_positions == 0 || myers_return_code != MYERS_STATUS_OK) {
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
  if (current_alignment) { free(current_alignment); } current_alignment = NULL;




  const unsigned char* reverse_query  = CreateReverseCopy(query, query_length);
  const unsigned char* reverse_target = CreateReverseCopy(target, (alignment_end + 1));

  int current_band_width = 0;

  myers_return_code = myersCalcEditDistance(reverse_query, query_length, reverse_target, (alignment_end + 1),
                                                128, -1, MYERS_MODE_SHW, &current_score, &current_positions, &current_num_positions,
                                                false, &current_alignment, &current_alignment_length, &current_band_width);
  if (current_num_positions == 0 || myers_return_code != MYERS_STATUS_OK) {
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
  typedef seqan::Dna5String TSequence;
  typedef seqan::StringSet<TSequence> TStringSet;
  typedef seqan::Gaps<TSequence, seqan::ArrayGaps> TGaps;
  typedef seqan::Iterator<TGaps>::Type TGapsIterator;
  TGaps gapsPattern(row(align, 0));
  TGaps gapsText(row(align, 1));
  TGapsIterator itGapsPattern = seqan::begin(row(align, 1));
  TGapsIterator itGapsEnd = seqan::end(row(align, 1));

  // Remove trailing gaps in pattern.
  int count = 0;
  while(isGap(--itGapsEnd))
      ++count;
  setClippedEndPosition(gapsPattern, length(gapsPattern) - count);
  // Remove leading gaps in pattern.
  if(isGap(itGapsPattern)) {
      setClippedBeginPosition(gapsPattern, countGaps(itGapsPattern));
      setClippedBeginPosition(gapsText, countGaps(itGapsPattern));
  }

  std::vector<unsigned char> alignment;
//  alignment.reserve(align.data_rows[0].

  // Reinitilaize the iterators.
  TGapsIterator itGapsText = seqan::begin(gapsText);
  itGapsPattern = begin(gapsPattern);
  itGapsEnd = seqan::end(gapsPattern);

  *ret_start_offset = gapsPattern._clippingBeginPos;
  *ret_end_offset = count;

  if (alignment_type == ALIGNMENT_TYPE_NW || alignment_type == ALIGNMENT_TYPE_SHW) {
    alignment.insert(alignment.begin(), (*ret_start_offset), (char) EDLIB_D);
    *ret_start_offset = 0;
  }

//  // Use a stringstream to construct the cigar string.
//  std::string cigar = "";

  int numChar = 0;
  while (itGapsPattern != itGapsEnd) {
      // Count insertions.
      if (isGap(itGapsText)) {
          int num_gaps = countGaps(itGapsText);
          alignment.insert(alignment.end(), num_gaps, (char) EDLIB_D);
          itGapsText += num_gaps;
          itGapsPattern += num_gaps;
          continue;
      }
      // Count deletions.
      if (isGap(itGapsPattern)) {
          int num_gaps = countGaps(itGapsPattern);
          alignment.insert(alignment.end(), num_gaps, (char) EDLIB_I);
          itGapsText += num_gaps;
          itGapsPattern += num_gaps;
          continue;
      }

      // Count matches.
      while (*itGapsText == *itGapsPattern && itGapsPattern != itGapsEnd) {
          ++numChar;
          alignment.insert(alignment.end(), 1, (char) EDLIB_EQUAL);
          ++itGapsText;
          ++itGapsPattern;
      }
      if (numChar != 0) {
//          cigar << numChar << "M";
          numChar = 0;
          continue;
      }

      // Count mismatches.
      while (*itGapsText != *itGapsPattern && itGapsPattern != itGapsEnd)
      {
          ++numChar;
          alignment.insert(alignment.end(), 1, (char) EDLIB_X);
          ++itGapsText;
          ++itGapsPattern;
      }
      if (numChar != 0) {
//          cigar << numChar << "X";
          numChar = 0;
          continue;
      }
  }

  if (alignment_type == ALIGNMENT_TYPE_NW) {
    alignment.insert(alignment.end(), (*ret_end_offset), (char) EDLIB_D);
    *ret_end_offset = 0;
  }

//  for (int64_t i=0; i<alignment.size(); i++) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }
//  for (int64_t i=(((int64_t) alignment.size()) - 1); i>=0; i--) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }

//  alignment.insert(alignment.end(), (*ret_end_offset), (char) EDLIB_D);
//  *ret_end_offset = 0;

  ret_alignment = alignment;

//  if (alignment.size() > 0)
//    cigar = AlignmentToCigar((unsigned char *) &alignment[0], alignment.size());
//  else
//    cigar = "*";
//
//  ret_cigar = cigar;

  return ALIGNMENT_GOOD;
}

int SeqAnSemiglobalWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
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
                           int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
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

//  if (true) {
//    std::string alignment_as_string = "";
//    alignment_as_string = PrintAlignmentToString((const unsigned char *) (read_data), read_length,
//                                               (const unsigned char *) (reference_data), (reference_length),
//                                               (unsigned char *) &(ret_alignment[0]), ret_alignment.size(),
//                                               (0), MYERS_MODE_NW);
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true,
//                                             FormatString("Alignment:\n%s\n\nalignment_position_start = %ld\n\n", alignment_as_string.c_str(), start_offset), "SeqAnNWWrapper");
//  }

//    printf ("start_offset = %d\n", start_offset);
//    fflush(stdout);

  *ret_alignment_position_start = start_offset;
  *ret_alignment_position_end = start_offset + (reconstructed_length - 1);
  *ret_edit_distance = (int64_t) seqan_edit_distance;

  return ALIGNMENT_GOOD;
}

int SeqAnSHWWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
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
  *ret_alignment_position_end = start_offset + (reconstructed_length - 1);
  *ret_edit_distance = (int64_t) seqan_edit_distance;

  return ALIGNMENT_GOOD;
}

int SeqAnSemiglobalWrapperWithMyersLocalization(const int8_t *read_data, int64_t read_length,
                                                const int8_t *reference_data, int64_t reference_length,
                                                int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                                                int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                                                int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Read length: %ld\n", read->get_sequence_length()), "SeqAnLocalRealignment");

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

  seqan::Infix<char *>::Type inf_target = seqan::infix((char *) (reference_data + localized_start), 0, (localized_end - localized_start + 1));
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

//  printf ("localized_start = %ld\n", localized_start);
//  printf ("localized_end = %ld\n", localized_end);
//  printf ("start_offset = %ld\n", start_offset);
//  printf ("end_offset = %ld\n", end_offset);
//  printf ("seqan_edit_distance = %ld\n", seqan_edit_distance);
//  fflush(stdout);

  *ret_alignment_position_start = start_offset + localized_start;
  *ret_alignment_position_end = start_offset + localized_start + (reconstructed_length - 1);
  *ret_edit_distance = (int64_t) localized_edit_distance;

  return ALIGNMENT_GOOD;
}

int MyersSemiglobalWrapper(const int8_t *read_data, int64_t read_length,
                           const int8_t *reference_data, int64_t reference_length,
                           int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                           int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                           int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Read length: %ld\n", read->get_sequence_length()), "SeqAnLocalRealignment");

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  int alphabet_length = 128;
  int score = 0;
  unsigned char* alignment = NULL;
  int alignment_length = 0;

  int *positions = NULL;
  int num_positions = 0;

  int myers_return_code = myersCalcEditDistance((const unsigned char *) read_data, read_length,
                        (const unsigned char *) reference_data, reference_length,
                        alphabet_length, band_width, MYERS_MODE_HW, &score, &positions, &num_positions,
                        true, &alignment, &alignment_length);

  if (myers_return_code == MYERS_STATUS_ERROR || num_positions == 0 || alignment_length == 0) {
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }

  int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

//  for (int64_t i=0; i<alignment_length; i++) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }
//  for (int64_t i=(((int64_t) alignment_length) - 1); i>=0; i--) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }

  *ret_alignment_position_start = positions[0] - (reconstructed_length - 1);
  *ret_alignment_position_end = positions[0];
  *ret_edit_distance = (int64_t) score;
  ret_alignment.assign(alignment, (alignment + alignment_length));
  //  *ret_cigar = AlignmentToCigar(alignment, alignment_length);

  if (positions)
    free(positions);

  if (alignment)
      free(alignment);

  return ALIGNMENT_GOOD;
}

int MyersNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Read length: %ld\n", read->get_sequence_length()), "SeqAnLocalRealignment");

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  int alphabet_length = 128;
  int score = 0;
  unsigned char* alignment = NULL;
  int alignment_length = 0;

  int *positions = NULL;
  int num_positions = 0;

  int myers_return_code = myersCalcEditDistance((const unsigned char *) read_data, read_length,
                        (const unsigned char *) reference_data, reference_length,
                        alphabet_length, band_width, MYERS_MODE_NW, &score, &positions, &num_positions,
                        true, &alignment, &alignment_length);

  if (myers_return_code == MYERS_STATUS_ERROR || num_positions == 0 || alignment_length == 0) {
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }

  int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

//  for (int64_t i=0; i<alignment_length; i++) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }
//  for (int64_t i=(((int64_t) alignment_length) - 1); i>=0; i--) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }

  *ret_alignment_position_start = positions[0] - (reconstructed_length - 1);
  *ret_alignment_position_end = positions[0];
  *ret_edit_distance = (int64_t) score;
  ret_alignment.assign(alignment, (alignment + alignment_length));
  //  *ret_cigar = AlignmentToCigar(alignment, alignment_length);

  if (positions)
    free(positions);

  if (alignment)
      free(alignment);

  return ALIGNMENT_GOOD;
}

int OpalNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t match_extend_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

//  mismatch_penalty = abs(mismatch_penalty);
  gap_open_penalty = abs(gap_open_penalty);
  gap_extend_penalty = abs(gap_extend_penalty);

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  uint8_t *converted_data = new uint8_t[read_length];
  if (converted_data == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_data."));
    return 1;
  }
  for (int64_t i=0; i<read_length; i++) {
    converted_data[i] = kBaseToBwaUnsigned[read_data[i]];
  }

  uint8_t *converted_ref = new uint8_t[reference_length];
  if (converted_ref == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_ref."));
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

//  printf ("match_score = %ld\n", match_score);
//  printf ("match_extend_score = %ld\n", match_extend_score);
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
                                                    gap_open_penalty, gap_extend_penalty, match_extend_score, scoreMatrix, alphabet_length, results,
                                                    OPAL_SEARCH_ALIGNMENT, OPAL_MODE_NW, OPAL_OVERFLOW_SIMPLE);
//  int resultCode = opalSearchDatabase(converted_data, read_length, db, db_length, db_seq_lengths,
//                                      gap_open_penalty, gap_extend_penalty, match_extend_score, scoreMatrix, alphabet_length, results,
//                                      OPAL_MODE_SW, OPAL_OVERFLOW_BUCKETS);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, "Finished running Opal.\n", "OpalNWWrapper");

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

int OpalSHWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t match_extend_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  uint8_t *converted_data = new uint8_t[read_length];
  if (converted_data == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_data."));
    return 1;
  }
  for (int64_t i=0; i<read_length; i++) {
    converted_data[i] = kBaseToBwaUnsigned[read_data[i]];
  }

  uint8_t *converted_ref = new uint8_t[reference_length];
  if (converted_ref == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: converted_ref."));
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
                                                    gap_open_penalty, gap_extend_penalty, match_extend_score, scoreMatrix, alphabet_length, results,
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

int MyersSHWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t band_width, int64_t match_score, int64_t mismatch_penalty, int64_t gap_open_penalty, int64_t gap_extend_penalty,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment) {

//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Read length: %ld\n", read->get_sequence_length()), "SeqAnLocalRealignment");

  if (read_data == NULL || reference_data == NULL || read_length <= 0 || reference_length <= 0)
    return ALIGNMENT_WRONG_DATA;

  int alphabet_length = 128;
  int score = 0;
  unsigned char* alignment = NULL;
  int alignment_length = 0;

  int *positions = NULL;
  int num_positions = 0;
  int found_k = 0;

//  int myers_return_code1 = myersCalcEditDistance((const unsigned char *) read_data, read_length,
//                        (const unsigned char *) reference_data, reference_length,
//                        alphabet_length, band_width, MYERS_MODE_SHW, &score, &positions, &num_positions,
//                        false, &alignment, &alignment_length, &found_k);
//  int myers_return_code = myersCalcEditDistance((const unsigned char *) read_data, read_length,
//                        (const unsigned char *) reference_data, reference_length,
//                        alphabet_length, found_k*2, MYERS_MODE_SHW, &score, &positions, &num_positions,
//                        true, &alignment, &alignment_length);
  int myers_return_code = myersCalcEditDistance((const unsigned char *) read_data, read_length,
                        (const unsigned char *) reference_data, reference_length,
                        alphabet_length, band_width, MYERS_MODE_SHW, &score, &positions, &num_positions,
                        true, &alignment, &alignment_length, &found_k);

//  printf ("read_length = %d\n", read_length);
//  printf ("reference_length = %d\n", reference_length);
//  printf ("band_width = %d\n", band_width);
//  printf ("found_k = %d\n", found_k);
//  printf ("score = %d\n", score);
//  for (int i=0; i<num_positions; i++) {
//    printf ("SHW position [%d] = %d\n", i, positions[i]);
//  }
//  fflush(stdout);

  if (myers_return_code == MYERS_STATUS_ERROR) {
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }
  if (num_positions == 0) {
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }
  if (alignment_length == 0) {
    return ALIGNMENT_MYERS_INTERNAL_ERROR;
  }

  int64_t reconstructed_length = CalculateReconstructedLength(alignment, alignment_length);

//  for (int64_t i=0; i<alignment_length; i++) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }
//  for (int64_t i=(((int64_t) alignment_length) - 1); i>=0; i--) {
//    if (alignment[i] == EDLIB_I)
//      alignment[i] = EDLIB_S;
//    else break;
//  }

  *ret_alignment_position_start = positions[0] - (reconstructed_length - 1);
  *ret_alignment_position_end = positions[0];
  *ret_edit_distance = (int64_t) score;
  ret_alignment.assign(alignment, (alignment + alignment_length));
  //  *ret_cigar = AlignmentToCigar(alignment, alignment_length);

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

    int myers_return_code = myersCalcEditDistance((const unsigned char *) read_data, read_length,
                          (const unsigned char *) reference_data, reference_length,
                          alphabet_length, -1, myers_mode_code, &score, &positions, &num_positions,
                          false, &alignment, &alignment_length);

    if (myers_return_code == MYERS_STATUS_ERROR || num_positions == 0) {
      return ALIGNMENT_MYERS_INTERNAL_ERROR;
    }

    if (ret_alignment_position_end != NULL)
      *ret_alignment_position_end = positions[0];
    if (ret_edit_distance != NULL)
      *ret_edit_distance = (int64_t) score;

    if (positions)
      free(positions);

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
