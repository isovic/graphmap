/*
 * alignment.cc
 *
 *  Created on: Jan 17, 2016
 *      Author: isovic
 */

#include "alignment/alignment.h"
#include "libs/opal.h"

int AlignRegion(const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, std::shared_ptr<is::Transcriptome> transcriptome, const ProgramParameters *parameters, const EValueParams *evalue_params, bool extend_to_end, PathGraphEntry *region_results) {
//  bool align_end_to_end = true;
//  bool spliced_alignment = true;
//  bool spliced_alignment = false;
//  bool align_end_to_end = parameters->alignment_algorithm != "spliced";
//  bool spliced_alignment = parameters->alignment_algorithm == "spliced";
//  bool align_end_to_end = !parameters->use_spliced || !parameters->use_split || (!parameters->disable_end_to_end);
  bool align_end_to_end = !parameters->use_spliced && !parameters->use_split && (!parameters->disable_end_to_end);
  bool spliced_alignment = parameters->use_spliced || parameters->use_split;

    if (parameters->alignment_algorithm == "sggotoh") {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using semiglobal alignment approach.\n", "Alignment");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Gotoh for alignment!\n", "Alignment");
      return SemiglobalAlignment(SeqAnSemiglobalWrapperWithMyersLocalization, read, index, parameters, evalue_params, region_results);

    } else if (parameters->alignment_algorithm == "sg") {

      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using semiglobal alignment approach.\n", "Alignment");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");

      return SemiglobalAlignment(MyersSemiglobalWrapper, read, index, parameters, evalue_params, region_results);

    } else if (parameters->alignment_algorithm == "anchor" || parameters->alignment_algorithm == "spliced") {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using anchored alignment approach.\n", "Alignment");
      bool is_linear = region_results->get_region_data().is_split == false || parameters->is_reference_circular == false;
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");
      return AnchoredAlignmentNew(MyersNWWrapper, MyersSHWWrapper, read, index, transcriptome, parameters, evalue_params, region_results, align_end_to_end, spliced_alignment);

    } else if (parameters->alignment_algorithm == "anchorgotoh" || parameters->alignment_algorithm == "splicedgotoh") {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using anchored alignment approach.\n", "Alignment");
      bool is_linear = region_results->get_region_data().is_split == false || parameters->is_reference_circular == false;
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Gotoh's algorithm for alignment!\n", "Alignment");
      return AnchoredAlignmentNew(SeqAnNWWrapper, SeqAnSHWWrapper, read, index, transcriptome, parameters, evalue_params, region_results, align_end_to_end, spliced_alignment);

#ifndef RELEASE_VERSION

    } else if (parameters->alignment_algorithm == "anchormex" || parameters->alignment_algorithm == "splicedmex") {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using anchored alignment approach.\n", "Alignment");
      bool is_linear = region_results->get_region_data().is_split == false || parameters->is_reference_circular == false;
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Match Extend algorithm for alignment!\n", "Alignment");
//      return AnchoredAlignmentNew(OpalNWWrapper, OpalSHWWrapper, read, index, parameters, evalue_params, region_results);
      return AnchoredAlignmentNew(OpalNWWrapper, NULL, read, index, transcriptome, parameters, evalue_params, region_results, false, spliced_alignment);
#endif

    } else {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using semiglobal alignment approach (default).\n", "Alignment");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters->debug_read, "Using Myers' bit-vector algorithm for alignment!\n", "Alignment");

      return SemiglobalAlignment(MyersSemiglobalWrapper, read, index, parameters, evalue_params, region_results);
    }

  return -10;
}

void VerboseAlignment(const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, const ProgramParameters *parameters, const AlignmentResults *aln) {
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    std::string alignment_as_string = "";

    alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data(), read->get_sequence_length(),
       (const unsigned char *) (&index->get_data()[0] + aln->raw_pos_start), (aln->raw_pos_end - aln->raw_pos_start + 1),
       (unsigned char *) &(aln->raw_alignment[0]), aln->raw_alignment.size(),
       aln->raw_pos_end - aln->raw_pos_start, aln->aln_mode_code);
//       (aln->raw_pos_end - aln->aln_window_start), aln->aln_mode_code);

    LOG_DEBUG_SPEC("raw_pos_start = %ld\nraw_pos_end = %ld\nedit_distance = %ld\n",
       aln->raw_pos_start, aln->raw_pos_end, aln->edit_distance);
    LOG_DEBUG_SPEC("Alignment:\n%s\n",
       alignment_as_string.c_str());
    LOG_DEBUG_SPEC("MD string:\n%s\n\n", aln->md.c_str());
    LOG_DEBUG_SPEC("CIGAR string:\n%s\n\n", aln->cigar.c_str());
    LOG_DEBUG_SPEC("\nref_start = %ld\nref_end = %ld\nraw_pos_start = %ld\nraw_pos_end = %ld\nquery_start = %ld\nquery_end = %ld\nreg_pos_start = %ld\nreg_pos_end = %ld\norientation = %s\nread_len = %ld\nis_aligned = %s\n",
                   aln->ref_start, aln->ref_end, aln->raw_pos_start, aln->raw_pos_end, aln->query_start, aln->query_end, aln->reg_pos_start, aln->reg_pos_end, (aln->orientation == kForward) ? "fwd" : "rev", read->get_sequence_length(), (aln->is_aligned == true) ? "true" : "false");
  }
}

/// Determines the start and end locations for semiglobal alignment locally to within the region. Region coordinates will be known from outside.
int GetL1PosInRegion(const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, const ProgramParameters *parameters, const PathGraphEntry *region_results,
                     int64_t *l1_start, int64_t *l1_end) {
  const Region &region = region_results->get_region_data();

  SeqOrientation orientation = (region_results->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t abs_ref_id = region_results->get_region_data().reference_id;
  int64_t ref_start = index->get_reference_starting_pos()[abs_ref_id];
  int64_t ref_length = index->get_reference_lengths()[abs_ref_id];

  int64_t l1_ref_start = 0, l1_ref_end = 0, reference_data_length = 0;

  if (region.is_split == false || parameters->is_reference_circular == false) {
    l1_ref_start = region_results->get_l1_data().l1_lmin;
    l1_ref_end = ((int64_t) (region_results->get_l1_data().l1_k * read->get_sequence_length())) + region_results->get_l1_data().l1_lmax;

    /// Check the bounds.
    if (l1_ref_start < ref_start)
      l1_ref_start = ref_start;
    if (l1_ref_end >= (ref_start + ref_length))
      l1_ref_end = ref_start + ref_length - 1;

    reference_data_length = l1_ref_end - l1_ref_start + 1;

    *l1_start = l1_ref_start - region.start;
    *l1_end = l1_ref_end - region.start;

  } else {
    l1_ref_start = region_results->get_l1_data().l1_lmin;
    l1_ref_end = ((int64_t) (region_results->get_l1_data().l1_k * read->get_sequence_length())) + region_results->get_l1_data().l1_lmax;
    int64_t region_length_joined = (region.end - region.start + 1) + (region.split_end - region.split_start + 1);

    /// Check the bounds.
    if (l1_ref_start < 0)
      l1_ref_start = 0;
    if (l1_ref_end >= region_length_joined)
      l1_ref_end = region_length_joined - 1;

    reference_data_length = l1_ref_end - l1_ref_start + 1;

    *l1_start = l1_ref_start;
    *l1_end = l1_ref_end;
  }

  /// Debug output.
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LOG_DEBUG_SPEC("\nread_length = %ld\nreference_length = %ld\nl1_start = %ld\nl1_end = %ld\nreference_data_length = %ld\nreference_start = %ld\nreference_length = %ld\nabsolute_reference_id = %ld\norientation = %s\n\n",
              read->get_sequence_length(), ref_length, l1_ref_start, l1_ref_end, reference_data_length, ref_start, ref_length, abs_ref_id, ((orientation == kForward) ? "forward" : "reverse"));
    LOG_DEBUG_SPEC("\n%s\n", region_results->VerboseInfoToString().c_str());
  }

  return 0;
}

int FindCircularEnd(const std::vector<uint8_t> &alignment, int64_t pos_of_ref_end,
                    int64_t *ret_end_on_aln, int64_t *ret_end_on_read, int64_t *ret_end_on_ref,
                    int64_t *ret_start_on_aln, int64_t *ret_start_on_read, int64_t *ret_start_on_ref) {
  if (pos_of_ref_end < 0) return 1;

  int64_t pos_on_ref = 0, pos_on_read = 0;
  for (int64_t i=0; i<alignment.size(); i++) {
    if (alignment[i] == EDLIB_S) {
      pos_on_read += 1;
    } else if (alignment[i] == EDLIB_M || alignment[i] == EDLIB_EQUAL || alignment[i] == EDLIB_X) {
      if (pos_on_ref == pos_of_ref_end) {
        *ret_end_on_ref = pos_on_ref;
        *ret_end_on_read = pos_on_read;
        *ret_end_on_aln = i;
      } else if (pos_on_ref > pos_of_ref_end) {
        *ret_start_on_ref = pos_on_ref;
        *ret_start_on_read = pos_on_read;
        *ret_start_on_aln = i;
        break;
      }
      pos_on_read += 1;
      pos_on_ref += 1;
    } else if (alignment[i] == EDLIB_I) {
      pos_on_read += 1;
    } else if (alignment[i] == EDLIB_D) {
      pos_on_ref += 1;
    }
  }

  return 0;
}

// pos_of_ref_end is the distance from the beginning of the region to the last base of the reference, before the beginning of the reference is concatenated.
// ref_start is used to 'rewind' the right part of the circular alignment.
int SplitCircularAlignment(const AlignmentResults *aln, int64_t pos_of_ref_end, int64_t ref_start, int64_t ref_len, AlignmentResults *aln_l, AlignmentResults *aln_r) {

  if (aln->reg_pos_start > pos_of_ref_end) {
    *aln_r = *aln;

    // The entire alignment is on the right part of the circular region.
    aln_l->is_aligned = false;

    // Fix the starting (and ending) coordinates of the right alignment part, because they would currently
    // be placed outside of the reference coordinate space.
    aln_r->ref_start = aln_r->ref_start - ref_len;             // Update the reference alignment coordinates. This one will be converted to local ref coord.
    aln_r->ref_end = aln_r->ref_end - ref_len;                                // Update the reference alignment coordinates. This one will be converted to local ref coord.
    aln_r->raw_pos_start = aln_r->raw_pos_start - ref_len;     // Update the reference alignment coordinates. This one will be converted to local ref coord.
    aln_r->raw_pos_end = aln_r->raw_pos_end - ref_len;                       // Update the reference alignment coordinates. This one will be converted to local ref coord.

    return 0;

  } else if (aln->reg_pos_end <= pos_of_ref_end) {
    *aln_l = *aln;

    // The entire alignment is on the left part of the circular region.
    aln_r->is_aligned = false;

    return 0;
  }

  *aln_l = *aln;
  *aln_r = *aln;

  int64_t end_on_aln = 0, end_on_read = 0, end_on_ref = 0;
  int64_t start_on_aln = 0, start_on_read = 0, start_on_ref = 0;

  FindCircularEnd(aln->raw_alignment, pos_of_ref_end - aln->reg_pos_start,
                  &end_on_aln, &end_on_read, &end_on_ref,
                  &start_on_aln, &start_on_read, &start_on_ref);

  // Trim the left part of the alignment.
  aln_l->ref_end = aln_r->ref_start + end_on_ref;             // Update the reference alignment coordinates. This one will be converted to local ref coord.
  aln_l->raw_pos_end = aln_l->raw_pos_start + end_on_ref;     // Update the raw reference alignment coordinates. This one will be left unchanged.
  aln_l->reg_pos_end = aln_l->reg_pos_start + end_on_ref;     // Update the region alignment coordinates.
  aln_l->query_end = aln_l->query_start + end_on_read;        // Update the query coordinate.
  int64_t clip_l = aln_l->query_len - (end_on_read + 1);
  uint8_t *alignment_l = (uint8_t *) malloc(sizeof(uint8_t) * ((end_on_aln + 1) + clip_l + 1));
  memmove(&alignment_l[0], &aln->raw_alignment[0], end_on_aln + 1);
  memset(&alignment_l[end_on_aln+1], EDLIB_S, clip_l);
  // Copy the alignment to the raw_alignment member.
  aln_l->raw_alignment.clear();
  aln_l->raw_alignment.assign(alignment_l, alignment_l + ((end_on_aln + 1) + clip_l));
  // Update the final alignment of the left part.
  aln_l->alignment = aln_l->raw_alignment;
  if (aln_l->orientation == kReverse) ReverseArray(aln_l->alignment);

  // Trim the right part of the alignment.
  aln_r->ref_start = aln_r->ref_start + start_on_ref - ref_len;             // Update the reference alignment coordinates. This one will be converted to local ref coord.
  aln_r->ref_end = aln_r->ref_end - ref_len;                                // Update the reference alignment coordinates. This one will be converted to local ref coord.
  aln_r->raw_pos_start = aln_r->raw_pos_start + start_on_ref - ref_len;     // Update the reference alignment coordinates. This one will be converted to local ref coord.
  aln_r->raw_pos_end = aln_r->raw_pos_end  - ref_len;                       // Update the reference alignment coordinates. This one will be converted to local ref coord.
  aln_r->reg_pos_start = aln_r->reg_pos_start + start_on_ref;               // Update the region alignment coordinates.
  aln_r->query_start = aln_r->query_start + start_on_read;                  // Update the query coordinate.
  int64_t clip_r = start_on_read;
  int64_t copy_r = (aln->raw_alignment.size() - start_on_aln);
  uint8_t *alignment_r = (uint8_t *) malloc(sizeof(uint8_t) * (copy_r + clip_r + 1));
  memset(&alignment_r[0], EDLIB_S, clip_r);
  memmove(&alignment_r[clip_r], &aln->raw_alignment[start_on_aln], copy_r);
  alignment_r[copy_r + clip_r] = '\0';
  // Copy the alignment to the raw_alignment member.
  aln_r->raw_alignment.clear();
  aln_r->raw_alignment.assign(alignment_r, alignment_r + copy_r + clip_r);
  // Update the final alignment of the left part.
  aln_r->alignment = aln_r->raw_alignment;
  if (aln_r->orientation == kReverse) ReverseArray(aln_r->alignment);

  int64_t l_ref = 0, l_read = 0;
  for (int64_t i=0; i<aln_r->raw_alignment.size(); i++) {
    if (aln_r->raw_alignment[i] == EDLIB_M || aln_r->raw_alignment[i] == EDLIB_X || aln_r->raw_alignment[i] == EDLIB_EQUAL || aln_r->raw_alignment[i] == EDLIB_D)
      l_ref += 1;
    if (aln_r->raw_alignment[i] == EDLIB_M || aln_r->raw_alignment[i] == EDLIB_X || aln_r->raw_alignment[i] == EDLIB_EQUAL || aln_r->raw_alignment[i] == EDLIB_I || aln_r->raw_alignment[i] == EDLIB_S)
      l_read += 1;
  }

  return 0;
}

// Checks if there is a strange (large) number of insertions and deletions, or consecutive insertion/deletion operations.
// SeqAn likes to make such alignments.
// Returns 0 if everything went ok.
int CheckAlignmentSane(std::vector<unsigned char> &alignment, const SingleSequence* read, std::shared_ptr<is::MinimizerIndex> index, int64_t reference_hit_id, int64_t reference_hit_pos) {
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
          // If there are insertions following deletions (or other way around), something is wrong again.
          if ((last_move == EDLIB_I && alignment_char == EDLIB_D) || (last_move == EDLIB_D && alignment_char == EDLIB_I)) {
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, true, FormatString("CheckAlignmentSane returned false! return 2. Consecutive I and D operations! last_move = %c, alignment_char = %c. num_same_moves = %ld, qname: '%s', read_length: %ld, ref_length: %ld\n", (last_move == EDLIB_I) ? 'I' : 'D', (alignment_char == EDLIB_I) ? 'I' : 'D', num_same_moves, read->get_header(), read_length, ref_length), "CheckAlignmentSane");
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
    LOG_DEBUG("CheckAlignmentSane returned false! return 3. Mismatch in length of the read determined from the alignmend and the actual length. Calculated read length = %ld (from alignment), read->get_sequence_length() = %ld.\n", read_length, read->get_sequence_length());
    return 3;
  }
  if ((index != NULL && reference_hit_id >= 0 && reference_hit_pos >= 0) && ref_length > index->get_reference_lengths()[reference_hit_id]) {
    LOG_DEBUG("CheckAlignmentSane returned false! return 4. Calculated reference length (from alignment) is longer than the actual reference. (ref_length = %ld, index->get_reference_lengths()[reference_hit_id] = %ld", ref_length, index->get_reference_lengths()[reference_hit_id]);
    return 4;
  }
  if ((index != nullptr && reference_hit_id >= 0 && reference_hit_pos >= 0) &&
      (reference_hit_pos + ref_length) > (index->get_reference_starting_pos()[reference_hit_id] + index->get_reference_lengths()[reference_hit_id])) {
    LOG_DEBUG("CheckAlignmentSane returned false! return 5. Alignment steps out of bounds of the reference.\n");
    return 5;
  }

  return 0;
}
