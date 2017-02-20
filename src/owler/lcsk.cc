/*
 * lcsk.cc
 *
 *  Created on: Feb 16, 2017
 *      Author: isovic
 */

#include "owler.h"
#include "algorithm/fenwick.h"

void Owler::WriteHits_(std::string out_path, const std::vector<uint128_t> &hits, int64_t hits_start, int64_t hits_end,
               int64_t ref_id, std::string read_header, int64_t read_length,
               std::string reference_header, int64_t reference_length,
               const std::vector<int32_t> *indices_to_output, const std::vector<int32_t> *cluster_ids) {

  std::vector<uint128_t> filtered_seed_hits;

  if (indices_to_output != NULL) {
    filtered_seed_hits.resize(indices_to_output->size());
    for (int64_t j = 0; j < indices_to_output->size(); j++) {
      filtered_seed_hits[j] = hits[(*indices_to_output)[j] + hits_start];
    }
  } else {
    filtered_seed_hits.resize((hits_end - hits_start + 1));
    for (int64_t j = 0; j < (hits_end - hits_start + 1); j++) {
      filtered_seed_hits[j] = hits[j + hits_start];
    }
  }

  FILE *fp = fopen(out_path.c_str(), "w");
  if (fp != NULL) {
    fprintf (fp, "%s\t0\t%ld\t%s\t0\t%ld\t0.0\n", read_header.c_str(), read_length, reference_header.c_str(), reference_length);
    for (int64_t j=0; j<filtered_seed_hits.size(); j++) {
      uint32_t seed_len = 12; // owler_data->seed_types[filtered_seed_hits[j].seed_type].length();
      int32_t cluster_id = 0;
      if (cluster_ids) {
        cluster_id = cluster_ids->at(j);
      }
      fprintf (fp, "%d\t%d\t%d\n", Owler::HitPosRead_(filtered_seed_hits[j]), Owler::HitPosRef_(filtered_seed_hits[j]), cluster_id);
      fprintf (fp, "%d\t%d\t%d\n", Owler::HitPosRead_(filtered_seed_hits[j]) + seed_len, Owler::HitPosRef_(filtered_seed_hits[j]) + seed_len, cluster_id);
    }
    fclose(fp);
  }
}

int Owler::CalcCoveredBases_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters,
                            const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int64_t seed_len, PairwiseOverlap &overlap) {

  int64_t cov_bases_A = 0, cov_bases_B = 0;

  int64_t A_end_prev = -1;
  int64_t B_end_prev = -1;
  for (int64_t i = 0; i < overlap.lcsk_indices.size(); i++) {
    int64_t A_start = Owler::HitPosRead_(hits[begin_hit+overlap.lcsk_indices[i]]);
    int64_t A_end = A_start + seed_len - 1;
    int64_t B_start = Owler::HitPosRef_(hits[begin_hit+overlap.lcsk_indices[i]]);
    int64_t B_end = B_start + seed_len - 1;

    cov_bases_A += std::min((A_end - A_end_prev), seed_len);
    cov_bases_B += std::min((B_end - B_end_prev), seed_len);

    if (cov_bases_A < 0 || cov_bases_B < 0) {
      LOG_DEBUG("ERROR: Number of covered bases should not be less than 0! i = %ld, std::min((A_end - A_end_prev), seed_length) = %ld < 0, A_end = %ld, A_end_prev = %ld, B_end = %ld, B_end_prev = %ld, seed_length = %ld\n", i, std::min((A_end - A_end_prev), seed_len), A_end, A_end_prev, B_end, B_end_prev, seed_len);
    }

    A_end_prev = A_end;
    B_end_prev = B_end;
  }

  overlap.cov_bases_query = cov_bases_A;
  overlap.cov_bases_target = cov_bases_B;

  return 0;
}


// Assumes vertices are sorted.
// If use_l1_filtering is true, then all vertices/anchors that have coordinates further than allowed_dist from the L1 line are filtered out.
// Otherwise, all vertices will be used.
// The L1 line is specified with k = 1 and l parameters (y = k*x + l).
//void Owler::CalcLCSFromLocalScoresCacheFriendly2_(OwlerData* owler_data, int64_t ref_hits_start, int64_t ref_hits_end, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices) {
int Owler::WrapLCSk_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters,
                       const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int32_t seed_len, PairwiseOverlap &overlap) {
  std::vector<uint128_t> events;
  std::vector<uint64_t> matches_starts;
  std::vector<uint64_t> matches_indices;
  int64_t max_seq_len = 0;

  LOG_DEBUG_SPEC("qid = %ld, tid = %ld, tid_fwd = %ld\n", overlap.qid, overlap.tid, overlap.tid_fwd);
  LOG_DEBUG_SPEC("num_hits = %ld\n", end_hit - begin_hit + 1);

  PrepareEvents_(hits, begin_hit, end_hit, 1, events, matches_starts, matches_indices, max_seq_len);

  int64_t lcsk_len = 0;
  std::vector<int32_t> raw_lcsk_indices;
  LCSk_(events, max_seq_len, 1, matches_starts, matches_indices, raw_lcsk_indices, lcsk_len);

  LOG_DEBUG_SPEC("raw_lcsk_indices.size() = %ld\n", raw_lcsk_indices.size());

//  std::vector<int32_t> lcsk_indices;
//  std::vector<int32_t> cluster_ids;
//  int32_t num_sv = 0;
  FilterColinear_(index, read, parameters, hits, begin_hit, end_hit, seed_len, raw_lcsk_indices, overlap.lcsk_indices, &(overlap.cluster_ids), overlap.num_sv);

//  overlap.lcsk_len = std::max((int64_t) 0, lcsk_len - seed_len);     // Reduce the length by the size of a seed, because the end point
//                                                                     // Also do not include the last seed.

  CalcCoveredBases_(index, read, parameters, hits, begin_hit, end_hit, seed_len, overlap);
  overlap.lcsk_len = overlap.cov_bases_target;

  LOG_DEBUG_SPEC("lcsk_indices.size() = %ld\n", overlap.lcsk_indices.size());

  // Debug output.
  #ifndef RELEASE_VERSION
      if (parameters->verbose_level > 5 && read->get_sequence_absolute_id() == parameters->debug_read) {
        std::vector<int32_t> *cluster_ids_pntr = (overlap.cluster_ids.size() > 0) ? (&(overlap.cluster_ids)) : (NULL);

        std::string reference_header = index->get_headers()[(overlap.tid % index->get_num_sequences_forward())];

        int64_t tid = overlap.tid % index->get_num_sequences_forward();
        std::string rev_suffix = (overlap.tid >= index->get_num_sequences_forward()) ? (std::string("-rev")) : (std::string("-fwd"));

        // All hits.
        WriteHits_(FormatString("temp/overlaps/scores-%ld%s.csv", tid, rev_suffix.c_str()), hits, begin_hit, end_hit, overlap.tid,
                  std::string(read->get_header()), read->get_sequence_length(), reference_header, index->get_reference_lengths()[overlap.tid], NULL, NULL);
        // After LCSk.
        WriteHits_(FormatString("temp/overlaps/LCS-%ld%s.csv", tid, rev_suffix.c_str()), hits, begin_hit, end_hit, overlap.tid,
                  std::string(read->get_header()), read->get_sequence_length(), reference_header, index->get_reference_lengths()[overlap.tid], &raw_lcsk_indices, NULL);
        // After LCSk and clustering.
        WriteHits_(FormatString("temp/overlaps/LCSL1-%ld%s.csv", tid, rev_suffix.c_str()), hits, begin_hit, end_hit, overlap.tid,
                  std::string(read->get_header()), read->get_sequence_length(), reference_header, index->get_reference_lengths()[overlap.tid], &(overlap.lcsk_indices), cluster_ids_pntr);
        // Just a space filler. Same as before.
        WriteHits_(FormatString("temp/overlaps/double_LCS-%ld%s.csv", tid, rev_suffix.c_str()), hits, begin_hit, end_hit, overlap.tid,
                  std::string(read->get_header()), read->get_sequence_length(), reference_header, index->get_reference_lengths()[overlap.tid], &(overlap.lcsk_indices), cluster_ids_pntr);
      }
  #endif

  if (overlap.lcsk_indices.size() == 0) {
    return 1;
  }

  return 0;
}

int Owler::PrepareEvents_(const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int64_t seed_len,
                          std::vector<uint128_t> &events, std::vector<uint64_t> &matches_starts, std::vector<uint64_t> &matches_indices, int64_t &max_seq_len) {
  uint32_t num_vertices = end_hit - begin_hit + 1;

  if (num_vertices <= 0) {
    return 1;
  }

  int64_t num_matches = 0;
  int64_t num_events = 0;
  int64_t lcskpp_length = 0;

  int32_t min_ref = HitPosRef_(hits[begin_hit]);
  int32_t min_query = HitPosRead_(hits[begin_hit]);

  for (uint32_t i=begin_hit; i<(end_hit+1); i++) {
    min_ref = std::min(min_ref, HitPosRef_(hits[i]));
    min_query = std::min(min_query, HitPosRead_(hits[i]));
  }

  uint32_t temp_max_seq_len = 0;

  events.clear();
  events.resize(num_vertices * 2);

  matches_starts.clear();
  matches_starts.resize(num_vertices);

  matches_indices.clear();
  matches_indices.resize(num_vertices);

  for (uint32_t i=0; i<num_vertices; i++) {
    int64_t hit_id = i + begin_hit;

    uint32_t ref_start = HitPosRef_(hits[hit_id]) - min_ref;
    uint32_t ref_end = ref_start + seed_len;
    uint32_t query_start = HitPosRead_(hits[hit_id]) - min_query;
    uint32_t query_end = query_start + seed_len;

    unsigned __int128 event1 = (((unsigned __int128) ref_start) << (8 * 8)) | (((unsigned __int128) query_start) << (4 * 8)) | (((unsigned __int128) (i + num_vertices)));
    events[num_events] = event1;
    num_events += 1;

    unsigned __int128 event2 = (((unsigned __int128) ref_end) << (8 * 8)) | (((unsigned __int128) query_end) << (4 * 8)) | ((((unsigned __int128) i)));
    events[num_events] = event2;
    num_events += 1;

    matches_starts[num_matches] = (uint64_t) (event1 >> (4 * 8));
    matches_indices[num_matches] = i;

    num_matches += 1;

    temp_max_seq_len = std::max(temp_max_seq_len, ref_end);
    temp_max_seq_len = std::max(temp_max_seq_len, query_end);
  }

  max_seq_len = (int64_t) temp_max_seq_len;
  std::sort(events.begin(), events.end());

  return 0;
}

// Events and matches should be sorted.
void Owler::LCSk_(std::vector<uint128_t> &events, int64_t n, int64_t k, std::vector<uint64_t> &matches_starts, std::vector<uint64_t> &matches_indices, std::vector<int32_t> &lcsk_indices, int64_t &lcsk_len) {
  int64_t num_events = events.size();
  int64_t num_matches = matches_starts.size();

  // Indexed by column, first:dp value, second:index in matches.
  FenwickMax<std::pair<int, int> > dp_col_max(n);
  std::vector<int> dp(num_matches);
  std::vector<int> recon(num_matches);
  std::vector<int> continues(num_matches, -1);

  for (int64_t curr = 0; curr < num_matches; curr++) {
    uint64_t G = 0;
    uint64_t G1 = (((matches_starts[curr] >> (4 * 8)) & (0x00000000FFFFFFFF)) - 1);
    uint64_t G2 = (((matches_starts[curr]) & (0x00000000FFFFFFFF)) - 1);
    G = (G1 << (4 * 8)) | G2;

    auto prev = std::lower_bound(matches_starts.begin(), (matches_starts.end()), G);

    if (prev != (matches_starts.end()) && *prev == G) {
      continues[curr] = prev - matches_starts.begin();
    }
  }

  int best_idx = 0;
  lcsk_len = 0;

  for (int64_t current_event=0; current_event<num_events; current_event++) {
    int64_t raw_idx = (int64_t) ((events[current_event] & (0x00000000FFFFFFFF)));

    int64_t idx = (raw_idx >= ((int64_t) num_matches)) ? (raw_idx - ((int64_t) num_matches)) : (raw_idx);
    bool is_beginning = (raw_idx >= num_matches);
    uint64_t i = (uint64_t) ((events[current_event] >> (8 * 8)) & (0x00000000FFFFFFFF));
    uint64_t j = (uint64_t) ((events[current_event] >> (4 * 8)) & (0x00000000FFFFFFFF));
//    int primary_diagonal = n - 1 + i - j;

    if (is_beginning) { // begin
      std::pair<int, int> prev_dp = dp_col_max.get(j);
      dp[idx] = k;
      recon[idx] = -1;

      if (prev_dp.first > 0) {
        dp[idx] = prev_dp.first + k;
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

      if (dp[idx] > lcsk_len) {
        lcsk_len = dp[idx];
        best_idx = idx;
      }
    }
  }

  lcsk_indices.clear();
  lcsk_indices.reserve(num_matches);

  if (best_idx != -1 && recon.size() > 0) {
    lcsk_indices.push_back(matches_indices[best_idx]);

    for (int i1 = best_idx; i1 != -1; i1 = recon[i1]) {
      if (recon[i1] != -1) {
        lcsk_indices.push_back(matches_indices[recon[i1]]);
      }
    }
  }
}

