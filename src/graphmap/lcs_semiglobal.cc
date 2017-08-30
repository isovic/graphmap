/*
 * postprocess_lcs.cc
 *
 *  Created on: Oct 5, 2014
 *      Author: ivan
 */

#include <algorithm>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "graphmap/graphmap.h"
#include "algorithm/fenwick.h"

// Index of an anchor is packed in the lower bytes of a 128-bit integer.
#define MASK_128_IDX        ((unsigned __int128) (0x000000000000000000000000FFFFFFFF))
// Then starting positions on the query.
#define MASK_128_STARTQ     ((unsigned __int128) (0x0000000000000000FFFFFFFF00000000))
// Then starting positions on the reference.
#define MASK_128_STARTR     ((unsigned __int128) (0x00000000FFFFFFFF0000000000000000))

/*
  Assumes vertices are sorted.
  If use_l1_filtering is true, then all vertices/anchors that have coordinates further than allowed_dist from the L1 line are filtered out.
  Otherwise, all vertices will be used.
  The L1 line is specified with k = 1 and l parameters (y = k*x + l).

  @allowed_anchor_overlap is the allowed number of overlapping bases of two neighboring anchors.
                          If > 0, each anchor's end coordinate will be subtracted by this number to achieve the result.
                          If < 0, each anchor's end coordinate will be equal to (start + 1). This is important for performing LCSk directly
                          on seed hits, as many of them can overlap, and LCSk would normally filter those out and reduce coverage.
*/
void GraphMap::CalcLCSFromLocalScoresCacheFriendly_(const Vertices *vertices, bool use_l1_filtering, int64_t l, int64_t allowed_dist,
                                                    int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices, int64_t allowed_anchor_overlap) {
  uint32_t num_vertices = vertices->num_vertices;

  if (num_vertices <= 0)
    return;

  unsigned __int128 *events = NULL;
  uint64_t *matches_starts = NULL;
  uint64_t *matches_dists_ref = NULL;
  uint64_t *matches_indices = NULL;

  uint32_t n = 0;
  int64_t num_matches = 0;

  int64_t num_events = 0;
  int64_t lcskpp_length = 0;

  uint32_t allowed_anchor_overlap_uint = (uint32_t) allowed_anchor_overlap;

  if (use_l1_filtering == false) {
    int64_t min_ref = vertices->reference_starts[0], min_query = vertices->query_starts[0];

    for (uint32_t i=0; i<num_vertices; i++) {

      min_ref = std::min(min_ref, vertices->reference_starts[i]);
      min_query = std::min(min_query, vertices->query_starts[i]);
    }

    events = (unsigned __int128 *) malloc(sizeof(unsigned __int128) * num_vertices * 2);
    matches_starts = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);
    matches_dists_ref = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);
    matches_indices = (uint64_t *) malloc(sizeof(uint64_t) * num_vertices);

    for (uint32_t i=0; i<num_vertices; i++) {
      uint32_t ref_start = vertices->reference_starts[i] - min_ref;
      uint32_t ref_end = vertices->reference_ends[i] - min_ref;
      uint32_t query_start = vertices->query_starts[i] - min_query;
      uint32_t query_end = vertices->query_ends[i] - min_query;
      uint32_t dist_ref = ref_end - ref_start;
      uint32_t dist_query = query_end - query_start;

      if (allowed_anchor_overlap > 0) {
        // Allowed anchor overlap by subtraction of the fuzz amount. Performs a sanity check for the size of elements because of unsigned int.
        query_end = (allowed_anchor_overlap_uint <= query_end) ? (std::max(query_start + 1, query_end - allowed_anchor_overlap_uint)) : query_end;
        ref_end = (allowed_anchor_overlap_uint <= ref_end) ? (std::max(ref_start + 1, ref_end - allowed_anchor_overlap_uint)) : ref_end;
      } else {  // This is useful for a large number of short seeds which may overlap. LCSk would otherwise remove many of overlapping seeds.
        query_end = query_start + 1;
        ref_end = ref_start + 1;
      }

      unsigned __int128 event1 = (((unsigned __int128) ref_start) << (8 * 8)) | (((unsigned __int128) query_start) << (4 * 8)) | (((unsigned __int128) (i + num_vertices)));
      events[num_events] = event1;
      num_events += 1;

      unsigned __int128 event2 = (((unsigned __int128) ref_end) << (8 * 8)) | (((unsigned __int128) query_end) << (4 * 8)) | ((((unsigned __int128) i)));
      events[num_events] = event2;
      num_events += 1;

      matches_starts[num_matches] = (uint64_t) (event1 >> (4 * 8));
      matches_dists_ref[num_matches] = std::max(dist_ref, dist_query);
      matches_indices[num_matches] = i;

      num_matches += 1;

      n = std::max(n, ref_end);
      n = std::max(n, query_end);
    }

  } else {
    int64_t num_filtered_vertices = 0;
    int64_t min_ref = 0, min_query = 0;
    bool is_min_initialized = false;

    for (uint32_t i=0; i<num_vertices; i++) {
      int64_t current_l_1 = vertices->reference_starts[i] - vertices->query_starts[i];
      int64_t current_l_2 = vertices->reference_ends[i] - vertices->query_ends[i];
      float distance1 = abs((float) ((current_l_1 - l)) * sqrt(2.0f)) / 2.0f;
      float distance2 = abs((float) ((current_l_2 - l)) * sqrt(2.0f)) / 2.0f;

      if (distance1 <= allowed_dist && distance2 <= allowed_dist) {
        num_filtered_vertices += 1;

        if (is_min_initialized == true) {
          min_ref = std::min(min_ref, vertices->reference_starts[i]);
          min_query = std::min(min_query, vertices->query_starts[i]);
        } else {
          min_ref = vertices->reference_starts[i];
          min_query = vertices->query_starts[i];
          is_min_initialized = true;
        }
      }
    }

    events = (unsigned __int128 *) malloc(sizeof(unsigned __int128) * num_filtered_vertices*2);
    matches_starts = (uint64_t *) malloc(sizeof(uint64_t) * num_filtered_vertices);
    matches_dists_ref = (uint64_t *) malloc(sizeof(uint64_t) * num_filtered_vertices);
    matches_indices = (uint64_t *) malloc(sizeof(uint64_t) * num_filtered_vertices);

    for (uint32_t i=0; i<num_vertices; i++) {
      int64_t current_l_1 = vertices->reference_starts[i] - vertices->query_starts[i];
      int64_t current_l_2 = vertices->reference_ends[i] - vertices->query_ends[i];
      float distance1 = abs((float) ((current_l_1 - l)) * sqrt(2.0f)) / 2.0f;
      float distance2 = abs((float) ((current_l_2 - l)) * sqrt(2.0f)) / 2.0f;

      if (distance1 <= allowed_dist && distance2 <= allowed_dist) {
        uint32_t ref_start = vertices->reference_starts[i] - min_ref;
        uint32_t ref_end = vertices->reference_ends[i] - min_ref;
        uint32_t query_start = vertices->query_starts[i] - min_query;
        uint32_t query_end = vertices->query_ends[i] - min_query;
        uint32_t dist_ref = ref_end - ref_start;
        uint32_t dist_query = query_end - query_start;

        unsigned __int128 event1 = (((unsigned __int128) ref_start) << (8 * 8)) | (((unsigned __int128) query_start) << (4 * 8)) | (((unsigned __int128) (num_matches + num_filtered_vertices)));
        events[num_events] = event1;
        num_events += 1;

        unsigned __int128 event2 = (((unsigned __int128) ref_end) << (8 * 8)) | (((unsigned __int128) query_end) << (4 * 8)) | ((((unsigned __int128) num_matches)));
        events[num_events] = event2;
        num_events += 1;

        matches_starts[num_matches] = (uint64_t) (event1 >> (4 * 8));
        matches_dists_ref[num_matches] = std::max(dist_ref, dist_query);
        matches_indices[num_matches] = i;

        num_matches += 1;

        n = std::max(n, ref_end);
        n = std::max(n, query_end);
      }
    }
  }

  std::sort(events, (events + num_events));

  // Indexed by column, first:dp value, second:index in matches.
  FenwickMax<std::pair<int, int> > dp_col_max(n);
  std::vector<int> dp(num_matches);
  std::vector<int> recon(num_matches);
  std::vector<int> continues(num_matches, -1);

//    matches bi trebao biti sortiran na pocetku!
//    sva sreca, pa su mi do sada entry-ji bili sortirani (vise-manje) po referenci...
//    ali, iznimka:
//
//    [59] (timestamp = 62287; q[8799, 8804]; r[18649998, 18650010]; d[5, 12]; length = 3; cov_bases_query = 11; cov_bases_ref = 13; registry_num = 59)
//    [60] (timestamp = 62289; q[9015, 9020]; r[18650004, 18650012]; d[5, 8]; length = 3; cov_bases_query = 11; cov_bases_ref = 13; registry_num = 60)
//    [61] (timestamp = 62306; q[10230, 10238]; r[18650023, 18650029]; d[8, 6]; length = 2; cov_bases_query = 12; cov_bases_ref = 12; registry_num = 61)
//    [62] (timestamp = 62308; q[1684, 1696]; r[18650026, 18650031]; d[12, 5]; length = 6; cov_bases_query = 16; cov_bases_ref = 11; registry_num = 62)
//    [63] (timestamp = 62310; q[1684, 1694]; r[18650026, 18650033]; d[10, 7]; length = 6; cov_bases_query = 16; cov_bases_ref = 13; registry_num = 63)
//    [64] (timestamp = 62311; q[10230, 10239]; r[18650023, 18650034]; d[9, 11]; length = 3; cov_bases_query = 13; cov_bases_ref = 13; registry_num = 64)
//    [65] (timestamp = 62310; q[1684, 1690]; r[18650026, 18650033]; d[6, 7]; length = 2; cov_bases_query = 12; cov_bases_ref = 12; registry_num = 65)



  for (int64_t curr = 0; curr < num_matches; curr++) {
    uint64_t G = 0;
    uint64_t G1 = (((matches_starts[curr] >> (4 * 8)) & (0x00000000FFFFFFFF)) - 1);
    uint64_t G2 = (((matches_starts[curr]) & (0x00000000FFFFFFFF)) - 1);
    G = (G1 << (4 * 8)) | G2;

    auto prev = std::lower_bound(matches_starts, (matches_starts + num_matches), G);

    if (prev != (matches_starts + num_matches) && *prev == G) {
      continues[curr] = prev - matches_starts;
    }
  }

  int best_idx = 0;
  lcskpp_length = 0;

  for (int64_t current_event=0; current_event<num_events; current_event++) {
    int64_t raw_idx = (int64_t) ((events[current_event] & (0x00000000FFFFFFFF)));

    int64_t idx = (raw_idx >= ((int64_t) num_matches)) ? (raw_idx - ((int64_t) num_matches)) : (raw_idx);
    bool is_beginning = (raw_idx >= num_matches);
    uint64_t i = (uint64_t) ((events[current_event] >> (8 * 8)) & (0x00000000FFFFFFFF));
    uint64_t j = (uint64_t) ((events[current_event] >> (4 * 8)) & (0x00000000FFFFFFFF));
    int primary_diagonal = n - 1 + i - j;

    if (is_beginning) { // begin
      std::pair<int, int> prev_dp = dp_col_max.get(j);
      uint64_t k_length = matches_dists_ref[idx];
      dp[idx] = k_length;      // k
      recon[idx] = -1;

      if (prev_dp.first > 0) {
        dp[idx] = prev_dp.first + k_length;
        recon[idx] = prev_dp.second;
      }
    } else {
      if (continues[idx] != -1) {
        if (dp[continues[idx]] + 1 > dp[idx]) {
          dp[idx] = dp[continues[idx]] + 1;
          recon[idx] = continues[idx];
        }
      }

      dp_col_max.update(j, std::make_pair(dp[idx], idx));

      if (dp[idx] > lcskpp_length) {
        lcskpp_length = dp[idx];
        best_idx = idx;
      }
    }
  }

  ret_lcskpp_indices->clear();
  ret_lcskpp_indices->reserve(num_matches);

  if (best_idx != -1 && recon.size() > 0) {
    ret_lcskpp_indices->push_back(matches_indices[best_idx]);

    for (int i1 = best_idx; i1 != -1; i1 = recon[i1]) {
      if (recon[i1] != -1) {
        ret_lcskpp_indices->push_back(matches_indices[recon[i1]]);
      }
    }
  }

  *ret_lcskpp_length = lcskpp_length;

  if (events)
    free(events);
  events = NULL;
  if (matches_starts)
    free(matches_starts);
  matches_starts = NULL;
  if (matches_dists_ref)
    free(matches_dists_ref);
  matches_dists_ref = NULL;
  if (matches_indices)
    free(matches_indices);
  matches_indices = NULL;
}

int GraphMap::CalculateL1ParametersWithMaximumDeviation_(ScoreRegistry *local_score, std::vector<int> &lcskpp_indices, float maximum_allowed_deviation, int64_t *ret_k, int64_t *ret_l, float *ret_sigma_L2, float *ret_confidence_L1) {
  int64_t k = 1;

  if (lcskpp_indices.size() == 0)
    return 1;

  std::vector<int64_t> l_array(lcskpp_indices.size() * 2);

  for (uint64_t i = 0; i < lcskpp_indices.size(); i++) {
    l_array[i * 2 + 0] = local_score->get_registry_entries().reference_starts[lcskpp_indices[i]] - k * local_score->get_registry_entries().query_starts[lcskpp_indices[i]];
    l_array[i * 2 + 1] = local_score->get_registry_entries().reference_ends[lcskpp_indices[i]] - k * local_score->get_registry_entries().query_ends[lcskpp_indices[i]];
  }

  std::sort(l_array.begin(), l_array.end());

  // This is a good definition! I was in dillema should it be '- 1' on the left or '+ 1' on the right side of the average,
  // but the arrays are 0-based, so that would cause an exception.
  int64_t l_median = ((l_array.size() % 2) == 1) ? (l_array.at((l_array.size() + 1) / 2)) : ((l_array.at((l_array.size() / 2) - 1) + l_array.at((l_array.size() / 2))) / 2);
  *ret_k = k;
  *ret_l = l_median;

  float confidence_L1 = 0.0f, confidence_L2 = 0.0f;

  // Just count the number of points that are within the maximum allowed limits.
  // These limits are determined using the read length and the error rate, for the worst case
  // if all errors are insertions or deletions.
  int64_t num_points_under_max_dev_threshold = 0;
  for (uint64_t i = 0; i < l_array.size(); i++) {
    float distance = abs((float) (l_array[i] - l_median) * (sqrt(2.0f)) / 2.0f);

    if (distance <= maximum_allowed_deviation)
      num_points_under_max_dev_threshold += 1;
  }

  if (num_points_under_max_dev_threshold == 0)
    return 2;

  // Calculate the standard deviation and the confidence interval for only those
  // points that are within the maximum allowed limits.
  for (uint64_t i = 0; i < l_array.size(); i++) {
    float distance = abs((float) (l_array[i] - l_median) * (sqrt(2.0f)) / 2.0f);

    if (distance > maximum_allowed_deviation)
      continue;

    confidence_L1 += ((float) abs(distance)) / ((float) num_points_under_max_dev_threshold);
    confidence_L2 += (distance * distance);
  }

  l_array.clear();

  float sigma = sqrt((1.0f / (num_points_under_max_dev_threshold - 1.0f)) * confidence_L2);

  *ret_sigma_L2 = sigma;
  *ret_confidence_L1 = confidence_L1;

  return 0;
}

int GraphMap::SemiglobalPostProcessRegionWithLCS_(ScoreRegistry* local_score, MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || ((int64_t) read->get_sequence_id()) == parameters->debug_read), FormatString("Entering function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB] current_readid = %ld, current_local_score = %ld\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024), read->get_sequence_id(), local_score->get_scores_id()), "PostProcessRegionWithLCS_");

  int lcskpp_length = 0;
  std::vector<int> lcskpp_indices;

  #ifndef RELEASE_VERSION
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      VerboseLocalScoresToFile(FormatString("temp/local_scores/scores-%ld.csv", local_score->get_scores_id()), read, local_score, NULL, 0, 0, false);
    }
  #endif

//  CalcLCSFromLocalScores2(&(local_score->get_registry_entries()), false, 0, 0, &lcskpp_length, &lcskpp_indices);
  CalcLCSFromLocalScoresCacheFriendly_(&(local_score->get_registry_entries()), false, 0, 0, &lcskpp_length, &lcskpp_indices, 0);

  if (lcskpp_length == 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Current local scores: %ld, lcskpp_length == 0 || best_score == NULL\n", local_score->get_scores_id()), "PostProcessRegionWithLCS_");
    return 1;
  }

  // Find the L1 parameters (median line and the confidence intervals).
  double l_diff = read->get_sequence_length() * parameters->error_rate;
  double maximum_allowed_deviation = l_diff * sqrt(2.0f) / 2.0f;
  float sigma_L2 = 0.0f, confidence_L1 = 0.0f;
  int64_t k = 0, l = 0;
  // Actual L1 calculation.
  int ret_L1 = CalculateL1ParametersWithMaximumDeviation_(local_score, lcskpp_indices, maximum_allowed_deviation, &k, &l, &sigma_L2, &confidence_L1);
  // Sanity check.
  if (ret_L1) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("An error occured, L1 function returned with %ld!\n", ret_L1), "L1-PostProcessRegionWithLCS_");
    return 1;
  }
  float allowed_L1_deviation = 3.0f * confidence_L1;

  #ifndef RELEASE_VERSION
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("l_median = %ld\n", l), "PostProcessRegionWithLCS_-DoubleLCSk");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("allowed_L1_deviation = %f\n", allowed_L1_deviation), "PostProcessRegionWithLCS_-DoubleLCSk");
      VerboseLocalScoresToFile(FormatString("temp/local_scores/LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, 0, 0, false);
      VerboseLocalScoresToFile(FormatString("temp/local_scores/LCSL1-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, l, allowed_L1_deviation, true);
    }
  #endif

  lcskpp_indices.clear();

  // Call the LCSk again, only on the bricks within the L1 bounded window.
  CalcLCSFromLocalScoresCacheFriendly_(&(local_score->get_registry_entries()), true, l, allowed_L1_deviation, &lcskpp_length, &lcskpp_indices, 0);

  // Count the number of covered bases, and find the first and last element of the LCSk.
  int64_t indexfirst = -1;
  int64_t indexlast = -1;

  int64_t covered_bases = 0;
  int64_t covered_bases_query = 0, covered_bases_reference = 0;
  int64_t num_covering_kmers = 0;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Counting the covered bases and finding the first and the last brick index.\n"), "PostProcessRegionWithLCS_-DoubleLCSk");
  for (uint64_t i = 0; i < lcskpp_indices.size(); i++) {
    covered_bases_query += local_score->get_registry_entries().covered_bases_queries[lcskpp_indices[i]];
    covered_bases_reference += local_score->get_registry_entries().covered_bases_references[lcskpp_indices[i]];
    num_covering_kmers += local_score->get_registry_entries().num_kmers[lcskpp_indices[i]];
  }
  covered_bases = std::max(covered_bases_query, covered_bases_reference);

  if (lcskpp_indices.size() > 0) {
    indexfirst = lcskpp_indices.front();
    indexlast = lcskpp_indices.back();
  }

  // There are no valid graph paths! All scores were dismissed because of high deviation.
  // This is most likely a false positive.
  if (indexfirst == -1 || indexlast == -1) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("An error occured, indexfirst = %ld, indexlast = %ld\n", indexfirst, indexlast), "L1-PostProcessRegionWithLCS_");
    return 1;
  }

  ret_L1 = CalculateL1ParametersWithMaximumDeviation_(local_score, lcskpp_indices, maximum_allowed_deviation, &k, &l, &sigma_L2, &confidence_L1);


  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("After second LCS calculation (with L1 filtering):\n", local_score->get_scores_id()), "PostProcessRegionWithLCS_2_");
    for (int64_t i = 0; i < lcskpp_indices.size(); i++) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld] %s\n", i, local_score->get_registry_entries().VerboseToString(lcskpp_indices[i]).c_str()), "[]");
    }
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Checking index indices.\n"), "PostProcessRegionWithLCS_-DoubleLCSk");

    #ifndef RELEASE_VERSION
      VerboseLocalScoresToFile(FormatString("temp/local_scores/double_LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, l, 3.0f * confidence_L1, true);
    #endif
  }

  MappingResults mapping_info;
  mapping_info.lcs_length = lcskpp_length;
  mapping_info.cov_bases_query = covered_bases_query;
  mapping_info.cov_bases_ref = covered_bases_reference;
  mapping_info.cov_bases_max = covered_bases;
  mapping_info.query_coords.start = local_score->get_registry_entries().query_starts[indexlast];
  mapping_info.query_coords.end = local_score->get_registry_entries().query_ends[indexfirst];
  mapping_info.ref_coords.start = local_score->get_registry_entries().reference_starts[indexlast];
  mapping_info.ref_coords.end = local_score->get_registry_entries().reference_ends[indexfirst];
  mapping_info.num_covering_kmers = num_covering_kmers;
  mapping_info.deviation = confidence_L1;
  mapping_info.is_reverse = (local_score->get_region().reference_id >= indexes[0]->get_num_sequences_forward());
  mapping_info.local_score_id = local_score->get_scores_id();

  L1Results l1_info;
  l1_info.l1_l = l;
  l1_info.l1_k = 1.0f;
  l1_info.l1_lmin = ((double) l) - ((double) l_diff);
  l1_info.l1_lmax = ((double) l) + ((double) l_diff);
  l1_info.l1_confidence_abs = confidence_L1;
  l1_info.l1_std = sigma_L2;
  l1_info.l1_rough_start = l1_info.l1_k * 0 + l1_info.l1_lmin;
  l1_info.l1_rough_end = ((double) l1_info.l1_k * read->get_sequence_length()) + l1_info.l1_lmax;
  if (l1_info.l1_rough_start < indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id])
    l1_info.l1_rough_start = indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id];
  if (l1_info.l1_rough_end >= (indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id] + indexes[0]->get_reference_lengths()[local_score->get_region().reference_id]))
    l1_info.l1_rough_end = (indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id] + indexes[0]->get_reference_lengths()[local_score->get_region().reference_id]) - 1;

  CheckMinimumMappingConditions_(&mapping_info, &l1_info, indexes[0], read, parameters);

  PathGraphEntry *new_entry = new PathGraphEntry(indexes[0], read, parameters, (Region &) local_score->get_region(), &mapping_info, &l1_info);

  float ratio = new_entry->CalcDistanceRatio();
  float ratio_suppress = new_entry->CalcDistanceRatioSuppress();
//  if (ratio_suppress > (parameters->error_rate * 2)) {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Called unmapped, because ratio suppress is too high. new_entry->ratio_suppress = %.2f, max = %.2f.", ratio_suppress, (parameters->error_rate * 2)), "L1-PostProcessRegionWithLCS_");
//    delete new_entry;
//    return 1;
//  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "Adding new entry.\n", "L1-PostProcessRegionWithLCS_");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
  mapping_data->intermediate_mappings.push_back(new_entry);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "l = %ld, ", "L1-PostProcessRegionWithLCS_");


  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || read->get_sequence_id() == parameters->debug_read), FormatString("Exiting function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "PostProcessRegionWithLCS_");

  return 0;
}

int GraphMap::CheckMinimumMappingConditions_(MappingResults *mapping_data, L1Results *l1_data, std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters) {
  mapping_data->is_mapped = true;

//  if (mapping_data->cov_bases_query < 0.005f* read->get_sequence_length())
//    mapping_data->is_mapped = false;

  if (PathGraphEntry::CalcFPFilterStatic(*mapping_data, index, read, parameters) < 0.02f) {
    mapping_data->is_mapped = false;
    return 1;
  }

  return 0;
}

int GraphMap::VerboseLocalScoresToFile(std::string file_path, const SingleSequence *read, const ScoreRegistry *local_score, const std::vector<int> *indices, int64_t l_median, float maximum_allowed_deviation, bool check_median_filtering, std::vector<int32_t> *cluster_ids) {
  FILE *fp = NULL;

  fp = fopen(file_path.c_str(), "w");
  if (fp == NULL) {
    return 1;
  }

  fprintf (fp, "%s\t%ld\t%ld\t%d\t%f\t%f\n", read->get_header(), read->get_sequence_id(), read->get_sequence_length(), ((check_median_filtering == true) ? 1 : 0), ((float) l_median), maximum_allowed_deviation);

  if (indices != NULL) {
    for (int64_t i = 0; i < indices->size(); i++) {
      if (check_median_filtering == false) {
        if (cluster_ids == NULL) {
          fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_starts[indices->at(i)], local_score->get_registry_entries().reference_starts[indices->at(i)]);
          fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_ends[indices->at(i)], local_score->get_registry_entries().reference_ends[indices->at(i)]);
        } else {
          fprintf (fp, "%ld\t%ld\t%d\n", local_score->get_registry_entries().query_starts[indices->at(i)], local_score->get_registry_entries().reference_starts[indices->at(i)], cluster_ids->at(i));
          fprintf (fp, "%ld\t%ld\t%d\n", local_score->get_registry_entries().query_ends[indices->at(i)], local_score->get_registry_entries().reference_ends[indices->at(i)], cluster_ids->at(i));
        }

      } else {
        float distance1 = abs((float) ((local_score->get_registry_entries().reference_starts[indices->at(i)] - local_score->get_registry_entries().query_starts[indices->at(i)]) - l_median) * (sqrt(2.0f)) / 2.0f);
        if (distance1 <= maximum_allowed_deviation) {
          if (cluster_ids == NULL) {
            fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_starts[indices->at(i)], local_score->get_registry_entries().reference_starts[indices->at(i)]);
          } else {
            fprintf (fp, "%ld\t%ld\t%d\n", local_score->get_registry_entries().query_starts[indices->at(i)], local_score->get_registry_entries().reference_starts[indices->at(i)], cluster_ids->at(i));
          }
        }

        float distance2 = abs((float) ((local_score->get_registry_entries().reference_ends[indices->at(i)] - local_score->get_registry_entries().query_ends[indices->at(i)]) - l_median) * (sqrt(2.0f)) / 2.0f);
        if (distance2 <= maximum_allowed_deviation) {
          if (cluster_ids == NULL) {
            fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_ends[indices->at(i)], local_score->get_registry_entries().reference_ends[indices->at(i)]);
          } else {
            fprintf (fp, "%ld\t%ld\t%d\n", local_score->get_registry_entries().query_ends[indices->at(i)], local_score->get_registry_entries().reference_ends[indices->at(i)], cluster_ids->at(i));
          }
        }
      }
    }

  } else {
    for (int64_t i = 0; i < local_score->get_registry_entries().num_vertices; i++) {
      if (check_median_filtering == false) {
        fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_starts[i], local_score->get_registry_entries().reference_starts[i]);
        fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_ends[i], local_score->get_registry_entries().reference_ends[i]);

      } else {
        float distance1 = abs((float) ((local_score->get_registry_entries().reference_starts[i] - local_score->get_registry_entries().query_starts[i]) - l_median) * (sqrt(2.0f)) / 2.0f);
        if (distance1 <= maximum_allowed_deviation) {
          fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_starts[i], local_score->get_registry_entries().reference_starts[i]);
        }

        float distance2 = abs((float) ((local_score->get_registry_entries().reference_ends[i] - local_score->get_registry_entries().query_ends[i]) - l_median) * (sqrt(2.0f)) / 2.0f);
        if (distance2 <= maximum_allowed_deviation) {
          fprintf (fp, "%ld\t%ld\n", local_score->get_registry_entries().query_ends[i], local_score->get_registry_entries().reference_ends[i]);
        }
      }
    }
  }

  fclose(fp);

  return 0;
}
