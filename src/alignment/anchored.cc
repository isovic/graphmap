/*
 * anchored.cc
 *
 *  Created on: Feb 16, 2016
 *      Author: isovic
 */

#include <alignment/alignment_wrappers.h>
#include "alignment/alignment.h"

#include "alignment/transcriptome_mod.h"

enum ClusterAlnType {
  kAlnBeginning = 0,
  kAlnCluster = 1,
  kAlnInBetween = 2,
  kAlnEnding = 3
};

int AlignFront(AlignmentFunctionType AlignmentFunctionSHW,
                const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, const ProgramParameters *parameters,
                const PathGraphEntry *region_results, int64_t ref_index_start, int64_t region_ref_start,
                const int8_t *ref_data, int64_t ref_len, int64_t clip_count_front,
                int64_t alignment_position_start, bool align_end_to_end, AlignmentResults &aln) {
  aln.aln_mode_code = EDLIB_MODE_SHW;
  aln.query_start = 0;
  aln.query_end = clip_count_front - 1;
  aln.raw_pos_start = alignment_position_start;
  aln.raw_pos_end = alignment_position_start;
  aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
  aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
  aln.ref_start = aln.raw_pos_start - ref_index_start;
  aln.ref_end = aln.raw_pos_end - ref_index_start;
  aln.is_aligned = true;
  aln.raw_alignment.clear();

//  int64_t absolute_reference_id = region_results->get_region_data().reference_id;
//  int64_t reference_id = region_results->get_region_data().reference_id;
//  int64_t reference_start = index->get_reference_starting_pos()[absolute_reference_id];
//  int64_t reference_length = index->get_reference_lengths()[absolute_reference_id];
  int64_t reference_start = 0;

  LOG_DEBUG_SPEC("Aligning the beginning of the read (overhang).\n");

  /// Check if we need to extend the alignment to the left boundary. Also, even if the user specified it, if we are to close to the boundary, just clip it.
  if (align_end_to_end == false) { // || ((alignment_position_start - clip_count_front*2) > 0 && (alignment_position_start - clip_count_front*2) < region_ref_start)) {
    std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
    aln.raw_alignment.insert(aln.raw_alignment.end(), insertions_front.begin(), insertions_front.end());
    aln.edit_distance = insertions_front.size();

//    LOG_DEBUG_SPEC("Skipping alignment of the beginning. align_end_to_end = %s, (alignment_position_start - clip_count_front*2) = %ld should be less than region_ref_start = %ld\n", (align_end_to_end == false) ? "false" : "true", ((alignment_position_start - clip_count_front*2)), region_ref_start);
    LOG_DEBUG_SPEC("Skipping alignment because end_to_end == false.\n");
  } else {
    LOG_DEBUG_SPEC("Performing the alignment.\n");

    /// Reversing the sequences to make the semiglobal alignment of the trailing and leading parts.
//    int8_t *reversed_query_front = reverse_data(read->get_data(), clip_count_front);
    LOG_DEBUG_SPEC("Reversing the query.\n");
    std::vector<int8_t> reversed_query_front;
    reverse_data2(read->get_data(), clip_count_front, reversed_query_front);

//    int8_t *reversed_ref_front = NULL;
    LOG_DEBUG_SPEC("Reversing the ref.\n");
    std::vector<int8_t> reversed_ref_front;
    int64_t reversed_ref_len = 0;
    if (clip_count_front*2  > (alignment_position_start - reference_start)) {
      LOG_DEBUG_SPEC("Option 1.\n");
      LOG_DEBUG_SPEC("alignment_position_start = %ld, reference_start = %ld\n", alignment_position_start, reference_start);
      reverse_data2(ref_data + reference_start, (alignment_position_start - reference_start), reversed_ref_front);
      reversed_ref_len = alignment_position_start - reference_start;
    } else {
      LOG_DEBUG_SPEC("Option 2.\n");
      reverse_data2(ref_data + (alignment_position_start - 1) - (clip_count_front*2 - 1), clip_count_front*2, reversed_ref_front);
      reversed_ref_len = clip_count_front*2;
    }

    int64_t bandwidth = -1;
//    bandwidth = 0.30f*read->get_sequence_length();

    LOG_DEBUG_SPEC("Running AlignmentFunctionSHW.\n");
    int64_t leftover_left_start = 0, leftover_left_end = 0, leftover_left_edit_distance = 0;
    std::vector<unsigned char> leftover_left_alignment;
    int ret_code_right = AlignmentFunctionSHW(&reversed_query_front[0], (clip_count_front),
                                              (int8_t *) (&reversed_ref_front[0]), reversed_ref_len,
                                              bandwidth, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
                                              &leftover_left_start, &leftover_left_end,
                                              &leftover_left_edit_distance, leftover_left_alignment);

    LOG_DEBUG_SPEC("Checking the retcode.\n");
    if (ret_code_right == 0) {
      if (parameters->verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters->debug_read) {
        std::string alignment_as_string = "";

        LOG_DEBUG_SPEC("End of the beginning part of the read: %ld\n", (alignment_position_start - 1) - (clip_count_front*2 - 1));
        alignment_as_string = PrintAlignmentToString((const unsigned char *) &reversed_query_front[0], clip_count_front,
                                                     (const unsigned char *) (&reversed_ref_front[0]), reversed_ref_len,
                                                     (unsigned char *) &(leftover_left_alignment[0]), leftover_left_alignment.size(),
                                                     (0), EDLIB_MODE_SHW);
        LOG_DEBUG_SPEC("Aligning the beginning of the read:\n%s\n", alignment_as_string.c_str());
        for (int64_t i=0; i<leftover_left_alignment.size(); i++) {
          LOG_DEBUG_SPEC_NO_HEADER("%ld", leftover_left_alignment[i]);
        }
        LOG_DEBUG_SPEC_NO_HEADER("\n\n\n");
      }
    }

    /// Check if the return code is ok. Otherwise, just clip the front.
    /// Added on 15.11.2015.: check if the edit distance of the front part is too high. EDlib will automatically return an error code, but SeqAn won't.
    /// An example is when the entire front part does not match (e.g. alignment of a read to a part of the reference consisted of N bases).
    if (ret_code_right != 0 || leftover_left_edit_distance > leftover_left_alignment.size()/1.9f) {
      // TODO: This is a nasty hack. EDlib used to crash when query and target are extremely small, e.g. query = "C" and target = "TC".
      // In this manner we just ignore the leading part, and clip it.
      std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
      aln.raw_alignment.insert(aln.raw_alignment.end(), insertions_front.begin(), insertions_front.end());
      aln.edit_distance = insertions_front.size();

      if (ret_code_right == 0) {
        LOG_DEBUG_SPEC("Not using the alignment because edit distance too big (edit_distance = %ld, max = %ld), clipping the end instead!\n", leftover_left_edit_distance, clip_count_front/2);
      } else {
        LOG_DEBUG_SPEC("Not using the alignment because AlignmentFunctionSHW returned with an error! ret_code_right = %ld\n", ret_code_right);
      }

    } else {
      if (leftover_left_alignment.size() == 0) {
        std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
        aln.raw_alignment.insert(aln.raw_alignment.end(), insertions_front.begin(), insertions_front.end());

      } else {
        unsigned char *reversed_alignment = reverse_data(&(leftover_left_alignment[0]), leftover_left_alignment.size());
        aln.raw_alignment.insert(aln.raw_alignment.end(), reversed_alignment, reversed_alignment + leftover_left_alignment.size());
        aln.edit_distance = leftover_left_edit_distance;
        alignment_position_start -= leftover_left_end + 1;
        aln.query_start = 0;
        aln.query_end = clip_count_front - 1;
        aln.raw_pos_start = alignment_position_start;     // alignment_position_start has been modified slightly above to point to the correct base.
        aln.raw_pos_end = alignment_position_start + leftover_left_end;     // TODO: Check this coordinate!! Might be missing a '+ 1';
        aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
        aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
        aln.ref_start = aln.raw_pos_start - ref_index_start;
        aln.ref_end = aln.raw_pos_end - ref_index_start;
        if (reversed_alignment)
          free(reversed_alignment);
      }

//      if (reversed_query_front)
//        free(reversed_query_front);
//      if (reversed_ref_front)
//        free(reversed_ref_front);
    }
  }

  return 0;
}

int AlignAnchor(AlignmentFunctionType AlignmentFunctionNW,
                const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, const ProgramParameters *parameters,
                const PathGraphEntry *region_results, int64_t ref_index_start, int64_t region_ref_start,
                const int8_t *ref_data, int64_t ref_len, int64_t cluster_id, AlignmentResults &aln) {

  int64_t cluster_query_start = region_results->get_mapping_data().clusters[cluster_id].query.start;
  int64_t cluster_query_end = region_results->get_mapping_data().clusters[cluster_id].query.end;
  int64_t cluster_ref_start = region_results->get_mapping_data().clusters[cluster_id].ref.start;
  int64_t cluster_ref_end = region_results->get_mapping_data().clusters[cluster_id].ref.end;
  int64_t query_alignment_length = cluster_query_end - cluster_query_start + 1;
  int64_t ref_alignment_length = cluster_ref_end - cluster_ref_start + 1;

  aln.aln_mode_code = EDLIB_MODE_NW;
  aln.query_start = cluster_query_start;
  aln.query_end = cluster_query_end;
  aln.raw_pos_start = cluster_ref_start;
  aln.raw_pos_end = cluster_ref_end;
  aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
  aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
  aln.ref_start = aln.raw_pos_start - ref_index_start;
  aln.ref_end = aln.raw_pos_end - ref_index_start;
  aln.raw_alignment.clear();

//  if (query_alignment_length == 0 || ref_alignment_length == 0) { return 1; }

  int64_t bandwidth = -1;
//  bandwidth = 0.30f*read->get_sequence_length();

  int64_t anchor_alignment_position_start = 0, anchor_alignment_position_end = 0, anchor_edit_distance = 0;
  std::vector<unsigned char> anchor_alignment;
  int ret_code1 = AlignmentFunctionNW(read->get_data() + cluster_query_start, (query_alignment_length),
                                      (int8_t *) (ref_data + cluster_ref_start), (ref_alignment_length),
                                   bandwidth, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
                                   &anchor_alignment_position_start, &anchor_alignment_position_end,
                                   &anchor_edit_distance, anchor_alignment);
  if (ret_code1 != 0 || anchor_alignment.size() == 0) {
    LOG_DEBUG_SPEC("Alignment returned with error! ret_code1 = %d\n", ret_code1);
    return ret_code1;
  }

  if (parameters->verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters->debug_read) {
    std::string alignment_as_string = "";
    alignment_as_string = PrintAlignmentToString((const unsigned char *) (read->get_data() + cluster_query_start), query_alignment_length,
                                                 (const unsigned char *) (ref_data + cluster_ref_start), (ref_alignment_length),
                                                 (unsigned char *) &(anchor_alignment[0]), anchor_alignment.size(),
                                                 (0), EDLIB_MODE_NW);
    LOG_DEBUG_SPEC("Aligned anchor %d:\n%s\n", cluster_id, alignment_as_string.c_str());
    LOG_DEBUG_SPEC("aln.query_start = %ld\n", aln.query_start);
    LOG_DEBUG_SPEC("aln.query_end = %ld\n", aln.query_end);
    LOG_DEBUG_SPEC("aln.ref_start = %ld\n", aln.ref_start);
    LOG_DEBUG_SPEC("aln.ref_end = %ld\n", aln.ref_end);
    LOG_DEBUG_SPEC("aln.raw_start = %ld\n", aln.raw_pos_start);
    LOG_DEBUG_SPEC("aln.raw_end = %ld\n", aln.raw_pos_end);
  }

  aln.raw_alignment.insert(aln.raw_alignment.end(), anchor_alignment.begin(), anchor_alignment.end());
  aln.edit_distance = anchor_edit_distance;
  aln.is_aligned = true;

  return 0;
}

int AlignInBetweenAnchors(AlignmentFunctionType AlignmentFunctionNW,
                          const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, const ProgramParameters *parameters,
                          const PathGraphEntry *region_results, int64_t ref_index_start, int64_t region_ref_start,
                          const int8_t *ref_data, int64_t ref_len, int64_t cluster_id, AlignmentResults &aln) {

  int64_t cluster_query_start = region_results->get_mapping_data().clusters[cluster_id].query.start;
  int64_t cluster_query_end = region_results->get_mapping_data().clusters[cluster_id].query.end;
  int64_t cluster_ref_start = region_results->get_mapping_data().clusters[cluster_id].ref.start;
  int64_t cluster_ref_end = region_results->get_mapping_data().clusters[cluster_id].ref.end;
  int64_t query_alignment_length = cluster_query_end - cluster_query_start + 1;
  int64_t ref_alignment_length = cluster_ref_end - cluster_ref_start + 1;

  if ((cluster_id + 1) < region_results->get_mapping_data().clusters.size()) {
    int64_t next_query_start = region_results->get_mapping_data().clusters[cluster_id+1].query.start;
    int64_t next_query_end = region_results->get_mapping_data().clusters[cluster_id+1].query.end;
    int64_t next_ref_start = region_results->get_mapping_data().clusters[cluster_id+1].ref.start;
    int64_t next_ref_end = region_results->get_mapping_data().clusters[cluster_id+1].ref.end;
    int64_t inbetween_query_length = (next_query_start - (cluster_query_end)) - 1;  /////I
    int64_t inbetween_ref_length = (next_ref_start - (cluster_ref_end)) - 1;  /////I

    aln.aln_mode_code = EDLIB_MODE_NW;
    aln.query_start = cluster_query_end + 1;
    aln.query_end = next_query_start - 1;
    aln.raw_pos_start = cluster_ref_end + 1;
    aln.raw_pos_end = next_ref_start - 1;
    aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
    aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
    aln.ref_start = aln.raw_pos_start - ref_index_start;
    aln.ref_end = aln.raw_pos_end - ref_index_start;
    aln.raw_alignment.clear();

    LOG_DEBUG_SPEC("Aligning in between anchors.\n");

    /// Check if there is actually any distance between the queries, or between the references.
    /// If there is no difference, that means there is a clean insertion/deletion.
    if (inbetween_query_length <= 0 && inbetween_ref_length > 0) {
      std::vector<unsigned char> deletions_inbetween(inbetween_ref_length, EDLIB_D);
      aln.raw_alignment.insert(aln.raw_alignment.end(), deletions_inbetween.begin(), deletions_inbetween.end());
      LOG_DEBUG_SPEC("Skipping alignment because: inbetween_query_length <= 0 && inbetween_ref_length > 0. inbetween_query_length = %ld, inbetween_ref_length = %ld\n", inbetween_query_length, inbetween_ref_length);
      aln.is_aligned = true;

    } else if (inbetween_ref_length <= 0 && inbetween_query_length > 0) {
      std::vector<unsigned char> insertions_inbetween(inbetween_query_length, EDLIB_I);
      aln.raw_alignment.insert(aln.raw_alignment.end(), insertions_inbetween.begin(), insertions_inbetween.end());
      LOG_DEBUG_SPEC("Skipping alignment because: inbetween_ref_length <= 0 && inbetween_query_length > 0. inbetween_query_length = %ld, inbetween_ref_length = %ld\n", inbetween_query_length, inbetween_ref_length);
      aln.is_aligned = true;

    } else if (inbetween_query_length < 0 && inbetween_ref_length < 0) {
      LOG_DEBUG_SPEC("Problem aligning in between anchors!\n");
      return ALIGNMENT_DISTANCE_BETWEEN_ANCHORS_PROBLEM;
      LOG_DEBUG_SPEC("Skipping alignment because: inbetween_query_length < 0 && inbetween_ref_length < 0. inbetween_query_length = %ld, inbetween_ref_length = %ld\n", inbetween_query_length, inbetween_ref_length);
      aln.is_aligned = true;

    } else if (inbetween_query_length > 0 && inbetween_ref_length > 0) {
      LOG_DEBUG_SPEC("Performing the alignment.\n");
      LOG_DEBUG_SPEC("Coordinates: [%ld, %ld]-[%ld, %ld]\n",
                     (cluster_query_end + 1), (cluster_ref_end + 1),
                     (cluster_query_end + 1 + inbetween_query_length - 1), (cluster_ref_end + 1 + inbetween_ref_length - 1));

//      if (inbetween_query_length == 0 || inbetween_ref_length == 0) { return 1; }
      int64_t bandwidth = -1;
//      bandwidth = 0.30f*read->get_sequence_length();

      int64_t between_alignment_position_start = 0, between_alignment_position_end = 0, between_anchor_edit_distance = 0;
      std::vector<unsigned char> between_anchor_alignment;
      int ret_code2 = AlignmentFunctionNW(read->get_data() + (cluster_query_end) + 1, inbetween_query_length,
                                          (int8_t *) (ref_data + cluster_ref_end) + 1, inbetween_ref_length,
                                       bandwidth, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
                                       &between_alignment_position_start, &between_alignment_position_end,
                                       &between_anchor_edit_distance, between_anchor_alignment);

      if (ret_code2 != 0 || between_anchor_alignment.size() == 0) {
        LOG_DEBUG_SPEC("Alignment returned with error! ret_code2 = %d\n", ret_code2);
        LOG_DEBUG_SPEC("inbetween_query_length = %ld\ninbetween_ref_length = %ld\nnext_ref_start = %ld\nref_end = %ld\n", inbetween_query_length, inbetween_ref_length, next_ref_start, cluster_ref_end);
        return ret_code2;
      }

      aln.raw_alignment.insert(aln.raw_alignment.end(), between_anchor_alignment.begin(), between_anchor_alignment.end());
      aln.edit_distance = between_anchor_edit_distance;
      aln.is_aligned = true;

      if (parameters->verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters->debug_read) {
        std::string alignment_as_string = "";
        alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + (cluster_query_end) + 1, inbetween_query_length,
                                                     (const unsigned char *) (ref_data + cluster_ref_end) + 1, inbetween_ref_length,
                                                     (unsigned char *) &(between_anchor_alignment[0]), between_anchor_alignment.size(),
                                                     (0), EDLIB_MODE_NW);
        LOG_DEBUG_SPEC("Aligning in between anchors %d and %d:\n%s\n", cluster_id, (cluster_id+1), alignment_as_string.c_str());
        for (int64_t i=0; i<between_anchor_alignment.size(); i++) {
          LOG_DEBUG_SPEC_NO_HEADER("%ld", between_anchor_alignment[i]);
        }
        LOG_DEBUG_SPEC_NO_HEADER("\n\n\n");
      }

    }
  }

  return 0;
}

int AlignBack(AlignmentFunctionType AlignmentFunctionSHW,
              const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, const ProgramParameters *parameters,
              const PathGraphEntry *region_results, int64_t ref_index_start, int64_t region_ref_start,
              const int8_t *ref_data, int64_t ref_len, int64_t clip_count_back, int64_t ref_start,
              int64_t alignment_position_end, bool align_end_to_end, AlignmentResults &aln) {
  aln.aln_mode_code = EDLIB_MODE_SHW;
  aln.query_start = read->get_sequence_length() - clip_count_back;
  aln.query_end = read->get_sequence_length() - 1;
  aln.raw_pos_start = alignment_position_end + 1;
  aln.raw_pos_end = alignment_position_end + 1;
  aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
  aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
  aln.ref_start = aln.raw_pos_start - ref_index_start;
  aln.ref_end = aln.raw_pos_end - ref_index_start;
  aln.raw_alignment.clear();

//  int64_t absolute_reference_id = region_results->get_region_data().reference_id;
//  int64_t reference_id = region_results->get_region_data().reference_id;
//  int64_t reference_start = index->get_reference_starting_pos()[absolute_reference_id];
//  int64_t reference_length = index->get_reference_lengths()[absolute_reference_id];
  int64_t reference_start = 0;
  int64_t reference_length = ref_len;

  LOG_DEBUG_SPEC("Aligning the end of the read (overhang).\n");

  /// Handle the clipping at the end, or extend alignment to the end of the sequence.
  if (align_end_to_end == false) { // || (alignment_position_end + 1 + clip_count_back * 2) >= (ref_start + ref_len)) {
      std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
      aln.raw_alignment.insert(aln.raw_alignment.end(), insertions_back.begin(), insertions_back.end());
      aln.edit_distance = insertions_back.size();
      aln.is_aligned = true;

      LOG_DEBUG_SPEC("Skipping alignment of the end because either end_to_end == false.\n");
//      LOG_DEBUG_SPEC("Skipping alignment of the end. align_end_to_end = %s, (alignment_position_end + 1 + clip_count_back * 2) = %ld should be less than (ref_start + ref_len) = %ld\n", (align_end_to_end == false) ? "false" : "true", (alignment_position_end + 1 + clip_count_back * 2), (ref_start + ref_len));
  } else {
    LOG_DEBUG_SPEC("Performing the alignment.\n");

    int64_t query_end = region_results->get_mapping_data().clusters.back().query.end; // + clip_count_back;  /////I
    int64_t leftover_right_start = 0, leftover_right_end = 0, leftover_right_edit_distance = 0;
    int64_t ref_len_for_aln = std::min(clip_count_back*2, (reference_start + reference_length - alignment_position_end - 1));

//    if (ref_len_for_aln == 0) {
//    }

    if (clip_count_back == 0) {
    }

    LOG_DEBUG_SPEC("clip_count_back = %ld, ref_len_for_aln = %ld\n", clip_count_back, ref_len_for_aln);

    int64_t bandwidth = -1;
//    bandwidth = 0.30f*read->get_sequence_length();

    std::vector<unsigned char> leftover_right_alignment;
    int ret_code_right = AlignmentFunctionSHW(read->get_data() + query_end + 1, (clip_count_back),
                                     (int8_t *) (ref_data + alignment_position_end + 1), ref_len_for_aln,
                                     bandwidth, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
                                     &leftover_right_start, &leftover_right_end,
                                     &leftover_right_edit_distance, leftover_right_alignment);

    if (ret_code_right == 0) {
      if (parameters->verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters->debug_read) {
        std::string alignment_as_string = "";
        alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + query_end + 1, clip_count_back,
                                                     (const unsigned char *) (ref_data + alignment_position_end + 1),
                                                     std::min(clip_count_back*2, (reference_start + reference_length - alignment_position_end - 1)),
                                                     (unsigned char *) &(leftover_right_alignment[0]), leftover_right_alignment.size(),
                                                     (0), EDLIB_MODE_SHW);
        LOG_DEBUG_SPEC("Aligning the end of the read:\n%s\n", alignment_as_string.c_str());
      }
    }

    /// Check if the return code is ok. Otherwise, just clip the back.
    /// Added on 15.11.2015.: check if the edit distance of the back part is too high. EDlib will automatically return an error code, but SeqAn won't.
    /// An example is when the entire back part does not match (e.g. alignment of a read to a part of the reference consisted of N bases).
    if (ret_code_right != 0 || leftover_right_edit_distance > clip_count_back/1.9f) {
      // TODO: This is a nasty hack. EDlib used to crash when query and target are extremely small, e.g. query = "C" and target = "TC".
      // In this manner we just ignore the trailing part, and clip it.
      std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
      aln.raw_alignment.insert(aln.raw_alignment.end(), insertions_back.begin(), insertions_back.end());
      aln.edit_distance = insertions_back.size();
      // Start and end positions need to be moved to one base prior the ones defined at the beginning of the function.
      // This is because there are no actual bases to move the start position to, all are insertions.
      aln.raw_pos_start = alignment_position_end;
      aln.raw_pos_end = alignment_position_end;
      aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
      aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
      aln.ref_start = aln.raw_pos_start - ref_index_start;
      aln.ref_end = aln.raw_pos_end - ref_index_start;
      aln.is_aligned = true;

      if (ret_code_right == 0) {
        LOG_DEBUG_SPEC("Not using the alignment because edit distance too big (edit_distance = %ld, max = %ld, clip_count_back = %ld), clipping the end instead!\n", leftover_right_edit_distance, clip_count_back/2, clip_count_back);
      } else {
        LOG_DEBUG_SPEC("Not using the alignment because AlignmentFunctionSHW returned with an error! ret_code_right = %ld\n", ret_code_right);
      }

    } else {
      if (leftover_right_alignment.size() == 0) {
        std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
        aln.raw_alignment.insert(aln.raw_alignment.end(), insertions_back.begin(), insertions_back.end());
        aln.edit_distance = insertions_back.size();
        // Start and end positions need to be moved to one base prior the ones defined at the beginning of the function.
        // This is because there are no actual bases to move the start position to, all are insertions.
        aln.raw_pos_start = alignment_position_end;
        aln.raw_pos_end = alignment_position_end;
        aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
        aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
        aln.ref_start = aln.raw_pos_start - ref_index_start;
        aln.ref_end = aln.raw_pos_end - ref_index_start;
        aln.is_aligned = true;

      } else {

        aln.raw_alignment.insert(aln.raw_alignment.end(), leftover_right_alignment.begin(), leftover_right_alignment.end());
        aln.edit_distance = leftover_right_edit_distance;
        alignment_position_end += leftover_right_end + 1;

        aln.query_start = read->get_sequence_length() - clip_count_back;
        aln.query_end = read->get_sequence_length() - 1;
        aln.raw_pos_start = alignment_position_end - leftover_right_end;
        aln.raw_pos_end = alignment_position_end;
        aln.reg_pos_start = aln.raw_pos_start - region_ref_start;
        aln.reg_pos_end = aln.raw_pos_end - region_ref_start;
        aln.ref_start = aln.raw_pos_start - ref_index_start;
        aln.ref_end = aln.raw_pos_end - ref_index_start;
      }

      aln.is_aligned = true;
    }
  }

  return 0;
}

int JoinAlignmentResults(std::vector<AlignmentResults> &alns, const SingleSequence *read, const int8_t *ref_data) {
  if (alns.size() == 0) { return 1; }

  int32_t first = 0, last = alns.size() - 1;
  for (first=0; first<alns.size(); first++) {
    if (alns[first].is_aligned == true) { break; }
  }
  for (last=(alns.size()-1); last>=0; last--) {
    if (alns[last].is_aligned == true) { break; }
  }

  AlignmentResults final_aln;
  final_aln = alns[first];
  final_aln.query_end = alns[last].query_end;
  final_aln.raw_pos_end = alns[last].raw_pos_end;
  final_aln.ref_end = alns[last].ref_end;
  final_aln.reg_pos_end = alns[last].reg_pos_end;
  final_aln.raw_alignment.clear();
  final_aln.alignment.clear();

  for (int32_t i=0; i<alns.size(); i++) {
    if (alns[i].is_aligned == false) { continue; }

    /// Check for a special case when previous global alignment ended with deletions or insertions, and the new one starts with deletions or insertions.
    /// Switching from deletions to insertions is basically a mismatch streak.
    if (final_aln.raw_alignment.size() > 0 && alns[i].raw_alignment.size() > 0 &&
        ((final_aln.raw_alignment.back() == EDLIB_D && alns[i].raw_alignment[0] == EDLIB_I) || (final_aln.raw_alignment.back() == EDLIB_I && alns[i].raw_alignment[0] == EDLIB_D))) {
      int64_t num_trailing_indels = 0;
      int64_t num_leading_indels = 0;
      int64_t current_op1 = final_aln.raw_alignment.size() - 1;
      while (current_op1 >= 0) {
        if ((current_op1 + 1) < final_aln.raw_alignment.size() && final_aln.raw_alignment[current_op1] != final_aln.raw_alignment[current_op1+1])
          break;
        num_trailing_indels += 1;
        current_op1 -= 1;
      }
      int64_t current_op2 = 0;
      while (current_op2 < alns[i].raw_alignment.size()) {
        if (current_op2 > 0 && alns[i].raw_alignment[current_op2] != alns[i].raw_alignment[current_op2-1])
          break;
        num_leading_indels += 1;
        current_op2 += 1;
      }

      int64_t min_count = std::min(num_trailing_indels, num_leading_indels);
      for (current_op1 = 0; current_op1 < min_count; current_op1++) {
//        if ((ref_data + alignment_position_end + 1 + leftover_right_start - current_op1 - 1) == (read->get_data() + (query_end + 1) - current_op1))
        if ((ref_data + alns[i].raw_pos_end - current_op1 - 1) == (read->get_data() + (alns[i].query_end + 1) - current_op1))
          final_aln.raw_alignment[final_aln.raw_alignment.size() - current_op1 - 1] = EDLIB_EQUAL;
        else
          final_aln.raw_alignment[final_aln.raw_alignment.size() - current_op1 - 1] = EDLIB_X;
      }

      final_aln.raw_alignment.insert(final_aln.raw_alignment.end(), alns[i].raw_alignment.begin() + min_count, alns[i].raw_alignment.end());
    } else {
      final_aln.raw_alignment.insert(final_aln.raw_alignment.end(), alns[i].raw_alignment.begin(), alns[i].raw_alignment.end());
    }
  }

  alns.clear();
  alns.push_back(final_aln);

  return 0;
}

int FillFrontAndBackInsertions(std::vector<AlignmentResults> &alns) {
  if (alns.size() == 0) { return 1; }

  for (int32_t i=0; i<alns.size(); i++) {
    if (alns[i].query_start > 0) {
      std::vector<int8_t> insertions_front(alns[i].query_start, EDLIB_I);
      alns[i].raw_alignment.insert(alns[i].raw_alignment.begin(), insertions_front.begin(), insertions_front.end());
    }
    if ((alns[i].query_end) < alns[i].query_len) {
      std::vector<int8_t> insertions_back(alns[i].query_len - alns[i].query_end - 1, EDLIB_I);
      alns[i].raw_alignment.insert(alns[i].raw_alignment.end(), insertions_back.begin(), insertions_back.end());
    }
  }

  return 0;
}

//int AnchoredAlignment2(bool is_linear, bool end_to_end, AlignmentFunctionType AlignmentFunctionNW, AlignmentFunctionType AlignmentFunctionSHW,
//                           const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
//                           int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part, int64_t *ret_nonclipped_left_part,
//                           int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part, int64_t *ret_nonclipped_right_part,
//                           SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
//                           int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op) {
int AnchoredAlignmentNew(AlignmentFunctionType AlignmentFunctionNW, AlignmentFunctionType AlignmentFunctionSHW,
                         const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, std::shared_ptr<is::Transcriptome> transcriptome,
                         const ProgramParameters *parameters,
                         const EValueParams *evalue_params, PathGraphEntry *region_results, bool align_end_to_end, bool spliced_alignment) {

  if (region_results->get_mapping_data().clusters.size() <= 0) {
    LOG_DEBUG_SPEC("No valid anchors exist!");
    return ALIGNMENT_WRONG_CLUSTER_SIZE;
  }

  LOG_DEBUG_SPEC("Entering anchored alignment.\n");

  const Region &region = region_results->get_region_data();
  SeqOrientation orientation = (region.reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);

  int64_t abs_ref_id = region.reference_id;
  int64_t ref_id = region.reference_id % index->get_num_sequences_forward();
  int64_t ref_data_start = index->get_reference_starting_pos()[abs_ref_id];
  int64_t ref_data_len = index->get_reference_lengths()[abs_ref_id];
  int64_t reference_length = index->get_reference_lengths()[abs_ref_id];
  int64_t region_length_joined = 0, start_offset = 0, position_of_ref_end = 0;

  int8_t *ref_data  = (int8_t *) &index->get_data()[0];       // The data of the region.
  int64_t region_ref_start = region.start;            // Position of the start of the region on the original Index data. E.g. reg_data[0] is the same base as index->get_data()[index_pos].
  int64_t pos_of_ref_end = 0; // If the region was circular, it crosses the boundary between the end and the start of the data. ref_data[index_pos_of_ref_end] is the last base of the reference before the split part is concatenated.
  bool is_cleanup_required = NULL;  // If true, the region_data will have to be freed manually.

  if (region.is_split == true) {
    LOG_DEBUG_SPEC("Concatenating regions for circular alignment.\n");
//    ConcatenateSplitRegion(index, (Region *) &(region), &reg_data, &region_length_joined, &start_offset, &position_of_ref_end);
//    ref_len = region_length_joined;

    // Here, a pointer to the beginning of the region is obtained. If the region was linear, reg_data just points to a location in the Index.
    // If the region was circular, reg_data points to a newly allocated memory for the concatenated region, and needs to be cleared.
    int retval_region = GetRegionData(index, &region, &ref_data, &ref_data_len, &region_ref_start, &pos_of_ref_end, &is_cleanup_required);
    ref_data_start = 0;

//    LOG_DEBUG_HIGH("\nConcatenating regions for circular alignment.\n");
  }



  std::vector<AlignmentResults> alns;     // Holds alignments for each cluster/anchor and inbetween clusters/anchors + the beginning and ending of the sequence.
  std::vector<ClusterAlnType> alns_anchor_type;      // Specifies whether an alignment is a cluster/anchor or an in-between part, or an overhang.
  int64_t num_alns = 0;
  alns.resize(region_results->get_mapping_data().clusters.size()*2 + 2);      // + 2 for the beginning/end overhang alignments.
  alns_anchor_type.resize(region_results->get_mapping_data().clusters.size()*2 + 2);  // + 2 for the beginning/end overhang alignments.

  AlignmentResults default_aln;
  default_aln.orientation = orientation;
  default_aln.is_reverse = (orientation == kReverse) ? true : false;
  default_aln.ref_id = ref_id;
  default_aln.ref_header = index->get_headers()[ref_id];
  default_aln.ref_len = reference_length;
  default_aln.query_id = read->get_sequence_absolute_id();
  default_aln.query_header = read->get_header();
  default_aln.query_len = read->get_sequence_length();
  default_aln.is_aligned = false;
  default_aln.aln_mode_code = EDLIB_MODE_NW;
  default_aln.edit_distance = -1;
  default_aln.raw_alignment.clear();
  default_aln.alignment.clear();
  default_aln.ref_start = -1;     // Starting position of the alignment within a reference (offset from the beginning of the reference).
  default_aln.ref_end = -1;       //
  default_aln.query_start = -1;   // Beginning of the alignment on the query. Clipped bases are not included.
  default_aln.query_end = -1;     //
  default_aln.raw_pos_start = -1; // The raw_pos_start then holds the absolute coordinate of the alignment in such joined sequence data.
  default_aln.raw_pos_end = -1;   //
  default_aln.reg_pos_start = -1; // Starting position of the alignment within a region (offset from the beginning of the region).
  default_aln.reg_pos_end = -1;   // End position of the alignment within a region (offset from the beginning of the region).
  // These values need to be set individually:
  //   Handled by CountAlignmentOperations below: &aln.num_eq_ops, &aln.num_x_ops, &aln.num_i_ops, &aln.num_d_ops, &aln.alignment_score, &aln.nonclipped_length
  //   Handled by each Align* function: aln.raw_alignment, aln.alignment, where aln.alignment is the same as raw_alignment if orientation == kForward, and reverse-complemented if orientation == kReverse.
  //                                    aln.ref_start, aln.ref_end, aln.query_start, aln.query_end, aln.raw_pos_start, aln.raw_pos_end, aln.reg_pos_start, aln.reg_pos_end

  for (int32_t i=0; i<alns.size(); i++) {
    alns[i] = default_aln;
  }

  int64_t clip_count_front = region_results->get_mapping_data().clusters.front().query.start;
  int64_t clip_count_back = read->get_sequence_length() - (region_results->get_mapping_data().clusters.back().query.end) - 1;  /////I

  int64_t alignment_position_start = region_results->get_mapping_data().clusters.front().ref.start; // - clip_count_front;
  int64_t alignment_position_end = region_results->get_mapping_data().clusters.back().ref.end; // + clip_count_back;  /////I

  // Verbose debug.
  if (parameters->verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters->debug_read) {
    for (int64_t i=0; i<region_results->get_mapping_data().clusters.size(); i++) {
      LOG_DEBUG_SPEC("[anchor %d] [%ld, %ld]-[%ld, %ld]\n", i, region_results->get_mapping_data().clusters[i].query.start, region_results->get_mapping_data().clusters[i].ref.start, region_results->get_mapping_data().clusters[i].query.end, region_results->get_mapping_data().clusters[i].ref.end);
    }
  }



  /// Aligning the begining of the read (in front of the first anchor).
  if (clip_count_front > 0) {
    int ret_code_front = AlignFront(AlignmentFunctionSHW, read, index, parameters, region_results, ref_data_start,
                                    region_ref_start, ref_data, ref_data_len, clip_count_front,
                                    alignment_position_start, align_end_to_end, alns[num_alns]);
    alns_anchor_type[num_alns] = kAlnBeginning;
    num_alns += 1;
    if (ret_code_front) { return ret_code_front; }
  }

  for (int64_t i=0; i<region_results->get_mapping_data().clusters.size(); i++) {
    if (parameters->verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters->debug_read) {
      LOG_DEBUG_SPEC("Aligning an anchor.\n");
      LOG_DEBUG_SPEC("[anchor %d] [%ld, %ld]-[%ld, %ld]\n", i, region_results->get_mapping_data().clusters[i].query.start, region_results->get_mapping_data().clusters[i].ref.start, region_results->get_mapping_data().clusters[i].query.end, region_results->get_mapping_data().clusters[i].ref.end);
    }

    ///////////////////////////
    /// Align the anchor.
    int ret_code_anchor = AlignAnchor(AlignmentFunctionNW, read, index, parameters, region_results, ref_data_start, region_ref_start, ref_data, ref_data_len, i, alns[num_alns]);
    alns_anchor_type[num_alns] = kAlnCluster;
    num_alns += 1;
    if (ret_code_anchor) { return ret_code_anchor; }

    ///////////////////////////
    ///////////////////////////
    /// Align in between the anchors.
    if ((i + 1) < region_results->get_mapping_data().clusters.size()) {
      int ret_code_inbetween = AlignInBetweenAnchors(AlignmentFunctionNW, read, index, parameters, region_results, ref_data_start, region_ref_start, ref_data, ref_data_len, i, alns[num_alns]);
      alns_anchor_type[num_alns] = kAlnInBetween;
      num_alns += 1;
      if (ret_code_inbetween) { return ret_code_inbetween; }
    }
  }

  /// Aligning the end of the read.
  if (clip_count_back > 0) {
    LOG_DEBUG_SPEC("Trying to align the end of the read. clip_count_back = %ld\n", clip_count_back);
    int ret_code_back = AlignBack(AlignmentFunctionSHW, read, index, parameters, region_results, ref_data_start, region_ref_start, ref_data, ref_data_len, clip_count_back, ref_data_start, alignment_position_end, align_end_to_end, alns[num_alns]);
//    if (!ret_code_back) {
      alns_anchor_type[num_alns] = kAlnEnding;
      num_alns += 1;
//    }
    if (ret_code_back) { return ret_code_back; }
  }

  // Check whether the alignment is supposed to be spliced or all together. If it's spliced, omit the non-cluster alignments.
  if (spliced_alignment == false) {
//    LOG_DEBUG_SPEC("Before joining:\n");
//    for (int32_t i5=alns.size()-3; i5<alns.size(); i5++) {
//      LOG_DEBUG_SPEC("i5 = %ld / %ld\n", (i5+1), alns.size());
//      LOG_DEBUG_SPEC("aln.query_start = %ld\n", alns[i5].query_start);
//      LOG_DEBUG_SPEC("aln.query_end = %ld\n", alns[i5].query_end);
//      LOG_DEBUG_SPEC("aln.ref_start = %ld\n", alns[i5].ref_start);
//      LOG_DEBUG_SPEC("aln.ref_end = %ld\n", alns[i5].ref_end);
//      LOG_DEBUG_SPEC("aln.raw_start = %ld\n", alns[i5].raw_pos_start);
//      LOG_DEBUG_SPEC("aln.raw_end = %ld\n", alns[i5].raw_pos_end);
//      LOG_DEBUG_SPEC("alns[alns.size()-3].is_aligned = %s\n", (alns[i5].is_aligned == true) ? "true" : "false");
//      LOG_DEBUG_SPEC("anchor_type = %s\n", (alns_anchor_type[i5] == kAlnBeginning) ? "beginning" : (alns_anchor_type[i5] == kAlnCluster) ? "cluster" : (alns_anchor_type[i5] == kAlnInBetween) ? "inbetween" : "ending");
//      LOG_DEBUG_SPEC_NEWLINE;
//    }
//    LOG_DEBUG_SPEC("alns.size() = %ld, num_alns = %ld\n", alns.size(), num_alns);

    JoinAlignmentResults(alns, read, ref_data);

//    LOG_DEBUG_SPEC("aln.query_start = %ld\n", alns[0].query_start);
//    LOG_DEBUG_SPEC("aln.query_end = %ld\n", alns[0].query_end);
//    LOG_DEBUG_SPEC("aln.ref_start = %ld\n", alns[0].ref_start);
//    LOG_DEBUG_SPEC("aln.ref_end = %ld\n", alns[0].ref_end);
//    LOG_DEBUG_SPEC("aln.raw_start = %ld\n", alns[0].raw_pos_start);
//    LOG_DEBUG_SPEC("aln.raw_end = %ld\n", alns[0].raw_pos_end);
//    LOG_DEBUG_SPEC_NEWLINE;

  } else {
    FillFrontAndBackInsertions(alns);

    for (int32_t i=0; i<alns.size(); i++) {
      if (alns_anchor_type[i] != kAlnCluster) {
        alns[i].is_aligned = false;
      }
    }
  }



//  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//    LOG_DEBUG_SPEC("alignment_position_start = %ld\nalignment_position_end = %ld\n", alignment_position_start, alignment_position_end);
//    std::string alignment_as_string = "";
//    alignment_as_string = PrintAlignmentToString((const unsigned char *) (read->get_data()), read->get_sequence_length(),
//                                               (const unsigned char *) (ref_data + alignment_position_start), (alignment_position_end - alignment_position_start + 1),
//                                               (unsigned char *) &(alignment[0]), alignment.size(),
//                                               (0), MYERS_MODE_NW);
//    LOG_DEBUG_SPEC("Alignment:\n%s\n\nalignment_position_start = %ld\n\n", alignment_as_string.c_str(), alignment_position_start);
//    LOG_DEBUG_SPEC("Alignment array:\n");
//    for (int i1=0; i1<alignment.size(); i1++) {
//      LOG_DEBUG_SPEC("%d", alignment[i1]);
//    }
//    LOG_DEBUG_SPEC_NEWLINE;
//  }

  for (int32_t i=0; i<alns.size(); i++) {
    AlignmentResults &aln = alns[i];
    aln.alignment = aln.raw_alignment;

  //  VerboseAlignment(read, index, parameters, &aln);
    CountAlignmentOperations((std::vector<unsigned char> &) aln.raw_alignment, read->get_data(), ref_data, abs_ref_id, aln.raw_pos_start, orientation,
                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend, true,
                             &aln.num_eq_ops, &aln.num_x_ops, &aln.num_i_ops, &aln.num_d_ops, &aln.alignment_score, &aln.edit_distance, &aln.nonclipped_length);

    LOG_DEBUG_SPEC("Calculating the E-value.\n");
    CalculateEValueDNA(aln.alignment_score, aln.nonclipped_length, index->get_data_length_forward(), evalue_params, &aln.evalue);
    LOG_DEBUG_SPEC("AS = %ld\n", aln.alignment_score);
    LOG_DEBUG_SPEC("aln.num_eq_ops = %ld\n", aln.num_eq_ops);
    LOG_DEBUG_SPEC("aln.num_x_ops = %ld\n", aln.num_x_ops);
    LOG_DEBUG_SPEC("aln.num_i_ops = %ld\n", aln.num_i_ops);
    LOG_DEBUG_SPEC("aln.num_d_ops = %ld\n", aln.num_d_ops);

    /// If alignment is linear, just add it. Otherwise, it needs to be split first.
    if (region.is_split == false || parameters->is_reference_circular == false) {
      LOG_DEBUG_SPEC("Alignment is linear.\n");
      if (orientation == kReverse) ReverseArray(aln.alignment);
      region_results->AddAlignmentData(aln);

    } else {
      LOG_DEBUG_SPEC("Alignment is circular and will be split. pos_of_ref_end = %ld, ref_len = %ld\n", pos_of_ref_end, ref_data_len);

      // These need to be fixed. The rest of the above code thinks that the region is linear and not circular.
      // In case the region is linear, anchors/clusters will have absolute reference coordinates.
      // The absolute coordinates allow the reference data to be accessed seemlesly.
      // For the circular alignment, the data needed to be copied into a new array, which is indexed from zero.
      // That's why the ref coordinates of anchors/clusters do not correspond to their reference positions. The region
      // information still contains still contains the region's start position on the reference, and this one will be
      // subtracted from the ref coordinates of anchors/clusters when calculating aln.reg_pos_start and aln.reg_pos_end.
      // For this reason we need to increase it back to obtain correct coodinates.
      aln.reg_pos_start += region_ref_start;
      aln.reg_pos_end += region_ref_start;
      aln.raw_pos_start += region_ref_start;
      aln.raw_pos_end += region_ref_start;
      aln.ref_start += region_ref_start;
      aln.ref_end += region_ref_start;

      AlignmentResults aln_l;
      AlignmentResults aln_r;
      SplitCircularAlignment(&aln, pos_of_ref_end, 0, reference_length, &aln_l, &aln_r);
//      LOG_DEBUG_SPEC("aln:\n");
//      LOG_DEBUG_SPEC("  aln.query_start = %ld\n", aln.query_start);
//      LOG_DEBUG_SPEC("  aln.query_end = %ld\n", aln.query_end);
//      LOG_DEBUG_SPEC("  aln.ref_data_start = %ld\n", aln.ref_start);
//      LOG_DEBUG_SPEC("  aln.ref_end = %ld\n", aln.ref_end);
//      LOG_DEBUG_SPEC("  aln.raw_start = %ld\n", aln.raw_pos_start);
//      LOG_DEBUG_SPEC("  aln.raw_end = %ld\n", aln.raw_pos_end);
//      LOG_DEBUG_SPEC("  aln.reg_pos_start = %ld\n", aln.reg_pos_start);
//      LOG_DEBUG_SPEC("  aln.reg_pos_end = %ld\n", aln.reg_pos_end);
//      LOG_DEBUG_SPEC("  aln.is_aligned = %s\n", (aln.is_aligned == true) ? "true" : "false");
//      LOG_DEBUG_SPEC("reference_length = %ld\n", reference_length);
//      LOG_DEBUG_SPEC_NEWLINE;
//      LOG_DEBUG_SPEC("aln_l:\n");
//      LOG_DEBUG_SPEC("  aln_l.query_start = %ld\n", aln_l.query_start);
//      LOG_DEBUG_SPEC("  aln_l.query_end = %ld\n", aln_l.query_end);
//      LOG_DEBUG_SPEC("  aln_l.ref_start = %ld\n", aln_l.ref_start);
//      LOG_DEBUG_SPEC("  aln_l.ref_end = %ld\n", aln_l.ref_end);
//      LOG_DEBUG_SPEC("  aln_l.raw_start = %ld\n", aln_l.raw_pos_start);
//      LOG_DEBUG_SPEC("  aln_l.raw_end = %ld\n", aln_l.raw_pos_end);
//      LOG_DEBUG_SPEC("  aln_l.reg_pos_start = %ld\n", aln_l.reg_pos_start);
//      LOG_DEBUG_SPEC("  aln_l.reg_pos_end = %ld\n", aln_l.reg_pos_end);
//      LOG_DEBUG_SPEC("  aln_l.is_aligned = %s\n", (aln_l.is_aligned == true) ? "true" : "false");
//      LOG_DEBUG_SPEC_NEWLINE;
//      LOG_DEBUG_SPEC("aln_r:\n");
//      LOG_DEBUG_SPEC("  aln_r.query_start = %ld\n", aln_r.query_start);
//      LOG_DEBUG_SPEC("  aln_r.query_end = %ld\n", aln_r.query_end);
//      LOG_DEBUG_SPEC("  aln_r.ref_start = %ld\n", aln_r.ref_start);
//      LOG_DEBUG_SPEC("  aln_r.ref_end = %ld\n", aln_r.ref_end);
//      LOG_DEBUG_SPEC("  aln_r.raw_start = %ld\n", aln_r.raw_pos_start);
//      LOG_DEBUG_SPEC("  aln_r.raw_end = %ld\n", aln_r.raw_pos_end);
//      LOG_DEBUG_SPEC("  aln_r.reg_pos_start = %ld\n", aln_r.reg_pos_start);
//      LOG_DEBUG_SPEC("  aln_r.reg_pos_end = %ld\n", aln_r.reg_pos_end);
//      LOG_DEBUG_SPEC("  aln_r.is_aligned = %s\n", (aln_r.is_aligned == true) ? "true" : "false");

      region_results->AddAlignmentData(aln_l);
      region_results->AddAlignmentData(aln_r);
//      if (aln_l.is_aligned == true && aln_r.is_aligned == true) {
//        LOG_DEBUG_HIGH("\nCircular!\n");
//      }
    }
  }

  /// Fill out statistics for each alignment (E-value calculation, couting of CIGAR operations, etc.) and check if the alignments are sane.
  for (int32_t i=0; i<region_results->get_alignments().size(); i++) {
    AlignmentResults *curr_aln = &region_results->get_alignments()[i];

    if (curr_aln->is_aligned == false) { continue; }

    /// Verbose debug.
    LOG_DEBUG_SPEC("Alignment part %d / %d:\n", (i + 1), region_results->get_alignments().size());

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
    curr_aln->ref_end = final_aln_pos_end; // % ref_data_len;

    FixAlignmentLeadingTrailingID(curr_aln->alignment, &curr_aln->ref_start, &curr_aln->ref_end);

    LOG_DEBUG_SPEC("Checking if the alignment is sane.\n");
    if (CheckAlignmentSane((std::vector<unsigned char> &) curr_aln->raw_alignment, read, index, curr_aln->ref_id + (curr_aln->is_reverse) ? index->get_num_sequences_forward() : 0, curr_aln->raw_pos_start) != 0) {
      curr_aln->is_aligned = false;
      LOG_DEBUG_SPEC("Alignment is insane!\n");
    } else {
      LOG_DEBUG_SPEC("Alignment is ok!\n");
    }

    if (parameters->is_transcriptome) {
      LOG_DEBUG_SPEC("Converting alignment from transcriptome space to genome space.\n");
      ConvertFromTranscriptomeToGenomeAln(parameters, index, transcriptome, curr_aln);
    }

    LOG_DEBUG_SPEC("Converting alignment to CIGAR string.\n");
    curr_aln->cigar = AlignmentToCigar((unsigned char *) &(curr_aln->alignment[0]), curr_aln->alignment.size(), parameters->use_extended_cigar);

    LOG_DEBUG_SPEC("Converting alignment to MD string.\n");
    curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) curr_aln->alignment, &index->get_data()[0], final_aln_pos_start);

//    printf ("final_aln_pos_start = %ld\n", final_aln_pos_start);
//    printf ("Query:\n%s\nTarget:\n%s\n", GetSubstring((char *) read->get_data(), read->get_data_length()).c_str(), GetSubstring((char *) index->get_data() + final_aln_pos_start, read->get_data_length()).c_str());
//    fflush(stdout);

    LOG_DEBUG_SPEC("Counting alignment operations.\n");
    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), &index->get_data()[0], abs_ref_id, curr_aln->raw_pos_start, orientation,
                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend, true,
                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->edit_distance, &curr_aln->nonclipped_length);
//    LOG_DEBUG_SPEC("Calculating the E-value.\n");
//    CalculateEValueDNA(curr_aln->alignment_score, curr_aln->nonclipped_length, index->get_data_length_forward(), evalue_params, &curr_aln->evalue);

    LOG_DEBUG_SPEC("Calculating alignment statistics.\n");
    double error_rate = ((double) curr_aln->num_x_ops + curr_aln->num_i_ops + curr_aln->num_d_ops) / ((double) curr_aln->nonclipped_length);
    double indel_error_rate = ((double) curr_aln->num_i_ops + curr_aln->num_d_ops) / ((double) curr_aln->nonclipped_length);
    if (error_rate > parameters->max_error_rate) {
      curr_aln->is_aligned = false;
      LOG_DEBUG_SPEC("Error rate is too high! error_rate = %lf, parameters->max_error_rate = %lf\n", error_rate, parameters->max_error_rate);
    }
    if (indel_error_rate > parameters->max_indel_error_rate) {
      curr_aln->is_aligned = false;
      LOG_DEBUG_SPEC("Indel error rate is too high! error_rate = %lf, parameters->max_error_rate = %lf\n", error_rate, parameters->max_error_rate);
    }

    VerboseAlignment(read, index, parameters, curr_aln);
  }

  if (is_cleanup_required && ref_data) { delete[] ref_data; }

  /// Count the number of 'unaligned' alignments.
  int32_t num_unaligned = 0;
  for (int32_t i=0; i<region_results->get_alignments().size(); i++) { if (region_results->get_alignments()[i].is_aligned == false) num_unaligned += 1; }
  if (num_unaligned == region_results->get_alignments().size()) { return ALIGNMENT_NOT_SANE; }

  return 0;
}

//int AnchoredAlignmentNew(AlignmentFunctionType AlignmentFunctionNW, AlignmentFunctionType AlignmentFunctionSHW, const SingleSequence *read, const Index *index, const ProgramParameters *parameters, const EValueParams *evalue_params, PathGraphEntry *region_results) {
/*
  /// General useful things.
  const Region &region = region;
  SeqOrientation orientation = (region.reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);
  int64_t abs_ref_id = region.reference_id;
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

  AlignmentResults aln;

//  parameters->extend_aln_to_end



//  GetL1PosInRegion(read, index, parameters, region_results, &l1_start, &l1_end);
//
//  if (retval_region) {
//    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "GetRegionData returned with error!"));
//    if (is_cleanup_required) { free(reg_data); }
//    return 1;
//  }
//
//  if (is_cleanup_required) { LOG_DEBUG_SPEC("Alignment is circular and manual cleanup will be required.\n"); }
//
//  /// Perform alignment and store results.
//  AlignmentResults aln;
//  std::vector<unsigned char> alignment1;
//  int64_t aln_pos_start = 0, aln_pos_end = 0;
//  int ret_code = AlignmentFunction(read->get_data(), read->get_sequence_length(),
//                                   reg_data + l1_start, (int64_t) (l1_end - l1_start),
//                                   -1, parameters->match_score, parameters->mex_score, -parameters->mismatch_penalty, -parameters->gap_open_penalty, -parameters->gap_extend_penalty,
//                                   &aln_pos_start, &aln_pos_end,
//                                   &aln.edit_distance,
//                                   aln.raw_alignment);
//
//  /// Sanity check.
//  if (ret_code != 0 || aln.raw_alignment.size() == 0) {
//    if (is_cleanup_required) { free(reg_data); }
//    LOG_DEBUG_SPEC("Something went wrong with AlignmentFunction.\n");
//    return ret_code;
//  }
//
//  int64_t aln_start_in_reg = l1_start + aln_pos_start;
//  int64_t aln_end_in_reg = l1_start + aln_pos_end;
//
//  int64_t aln_abs_start = aln_start_in_reg + index_pos;
//  int64_t aln_abs_end = aln_end_in_reg + index_pos;





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

//  VerboseAlignment(read, index, parameters, &aln);
  CountAlignmentOperations((std::vector<unsigned char> &) aln.raw_alignment, read->get_data(), reg_data, ref_id, aln.reg_pos_start, orientation,
                           parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend,
                           &aln.num_eq_ops, &aln.num_x_ops, &aln.num_i_ops, &aln.num_d_ops, &aln.alignment_score, &aln.nonclipped_length);
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
    region_results->AddAlignmentData(aln_l);
    region_results->AddAlignmentData(aln_r);
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
    curr_aln->ref_start = final_aln_pos_start; // % ref_len;
    curr_aln->ref_end = final_aln_pos_end; // % ref_len;

    LOG_DEBUG_SPEC("Converting alignment to CIGAR string.\n");
    curr_aln->cigar = AlignmentToCigar((unsigned char *) &(curr_aln->alignment[0]), curr_aln->alignment.size(), parameters->use_extended_cigar);

    LOG_DEBUG_SPEC("Converting alignment to MD string.\n");
    curr_aln->md = AlignmentToMD((std::vector<unsigned char> &) curr_aln->alignment, read->get_data(), index->get_data(), ref_id, curr_aln->ref_start + ref_start + index->get_reference_starting_pos()[ref_id]);

    LOG_DEBUG_SPEC("Counting alignment operations.\n");
    CountAlignmentOperations((std::vector<unsigned char> &) curr_aln->raw_alignment, read->get_data(), reg_data, ref_id, curr_aln->reg_pos_start, orientation,
                             parameters->evalue_match, parameters->evalue_mismatch, parameters->evalue_gap_open, parameters->evalue_gap_extend,
                             &curr_aln->num_eq_ops, &curr_aln->num_x_ops, &curr_aln->num_i_ops, &curr_aln->num_d_ops, &curr_aln->alignment_score, &curr_aln->nonclipped_length);
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

*/
//  return 0;
//}



