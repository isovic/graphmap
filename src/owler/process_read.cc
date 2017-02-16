/*
 * process_read.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: isovic
 */

#include <limits>
#include <algorithm>

#include "owler/owler.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include "algorithm/fenwick.h"



int Owler::ProcessRead_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data) {
  LOG_DEBUG_SPEC_NEWLINE;
  LOG_DEBUG_SPEC("Entered function.\n");

  // If the read length is too short, call it unmapped.
  if (read->get_sequence_length() < parameters->min_read_len) {
    std::stringstream ss;
    ss << "Unmapped_5__readlength_too_short" << "__readlength=" << read->get_sequence_length() << "__limit=" << 80;
    owler_data.unmapped_reason += ss.str();
    return 0;
  }

  int32_t diag_epsilon = 500;     // TODO: Parametrize this.

  TicToc tt_collect;
  tt_collect.start();
  CollectHits_(index_, read, parameters, owler_data);
  tt_collect.stop();
  LOG_DEBUG_SPEC("Time collecting hits: %f\n", tt_collect.get_msecs());

  TicToc tt_cluster;
  tt_cluster.start();
  ClusterHits_(index, read, parameters, diag_epsilon, owler_data);
  tt_cluster.stop();
  LOG_DEBUG_SPEC("Time clustering hits: %f\n", tt_cluster.get_msecs());

//  GenerateOutput_(index_, read, parameters, owler_data);

//  owler_data->seed_hits2.reserve(500000);
//
//  /// Check if it's a case of self-overlap. In this case, overlap can be performed faster, because the index will already have pre-processed seeds of all reads.
//  if (parameters->reads_path == parameters->reference_path)
//    CollectSeedHitsExperimentalSubseededIndex(owler_data, indexes, read, parameters);
//  else
//    CollectSeedHitsExperimentalCalcSubseedsFast(owler_data, indexes, read, parameters);
//
//  ApplyLCS2(owler_data, indexes, read, parameters);

  // Just verbose.
  LOG_DEBUG_SPEC_NEWLINE;
  LOG_DEBUG_SPEC("Exiting function.\n");

  return 0;
}

int Owler::CollectHits_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data) {
  TicToc tt_collect;

  int64_t read_len = read->get_sequence_length();
  int64_t num_fwd_seqs = index->get_num_sequences_forward();
  bool is_overlapper = (parameters->reference_path == parameters->reads_path);

  owler_data.hits.clear();
  owler_data.hits.reserve(index->avg_seed_occurrence() * read_len);

  std::vector<uint128_t> seeds;       // All seeds for a given index and a given sequence.
  std::vector<int8_t> seed_key_lens;  // Lengths of seed keys. Since multiple lookup keys can be used, each can be of different length.
  std::vector<int8_t> seed_lens;    // Since there are several
  int64_t qid = read->get_sequence_id();

//  index->CollectLookupSeeds(read->get_data(), read->get_quality(), read->get_sequence_length(), 0.0f, false, false, 1, seeds);
  index->CollectLookupSeeds(read->get_data(), read->get_quality(), read->get_sequence_length(), 0.0f, false, parameters->use_minimizers, parameters->minimizer_window, seeds);

  owler_data.hits.clear();

  for (int64_t j=0; j<seeds.size(); j++) {
    AppendSeedHits_(seeds[j], index, parameters->threshold_hits, index->count_cutoff(), is_overlapper, qid, owler_data.hits);
  }

  std::sort(owler_data.hits.begin(), owler_data.hits.end());

  return 0;
}

int Owler::ClusterHits_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, int32_t diag_epsilon, OwlerData &owler_data) {

  int64_t max_shape_width = index->get_shape_max_width();

  int64_t begin_hit = 0;
  int64_t total_num_hits = owler_data.hits.size();

  for (int64_t end_hit=0; end_hit<total_num_hits; end_hit++) {
    if ((end_hit + 1) == total_num_hits ||
        HitSeqId_(owler_data.hits[end_hit+1]) != HitSeqId_(owler_data.hits[end_hit]) ||
        (HitDiag_(owler_data.hits[end_hit+1]) - HitDiag_(owler_data.hits[end_hit])) >= diag_epsilon) {

      Range rhits(begin_hit, end_hit);     // Store the current range of hits.
      begin_hit = end_hit + 1;                // Update the beginning of the range right away in case there will be any 'continue' statements.

      int64_t num_hits = rhits.end - rhits.start + 1;
      int64_t cov_bases_estimate = num_hits * max_shape_width;

      if (cov_bases_estimate < 100) { continue; }

      // Add a cluster.
      int64_t ref_id = HitSeqId_(owler_data.hits[rhits.end]);
      int64_t ref_start = index->get_reference_starting_pos()[ref_id];
      int64_t ref_len = index->get_reference_lengths()[ref_id];
      int64_t read_len = read->get_sequence_length();

//      LOG_ALL("begin_hit = %ld, end_hit = %ld\n", begin_hit, end_hit);

      PairwiseOverlap overlap;
      int rv_lcsk = LCSkFilter_(owler_data.hits, rhits.start, rhits.end, read->get_sequence_id(), ref_id, max_shape_width, overlap);

      if (rv_lcsk) {      // Nothing survived the filter.
        continue;
      }

      if (!CheckOverlap_(index, read, parameters, overlap)) {
        std::string overlap_line = GenerateMHAPLine(index, read, parameters, overlap);
        owler_data.out_lines.push_back(overlap_line);
      }

//      int64_t start_on_read = 0, end_on_read = 0;
//      int64_t start_on_ref = 0, end_on_ref = 0;

//      FindRegionBounds_(owler_data.hits, begin_hit, end_hit, start_on_read, end_on_read, start_on_ref, end_on_ref);



//      PairwiseOverlap o;
//      o.query.start =
//
//      Region region;
//      region.start = std::max(ref_start + start - read_len, ref_start);
//      region.end = ref_start + std::min(end + max_shape_width + read_len, ref_len);; // TODO: instead of max_seed_len, the exact seed length (for a sepecific lookup key) should be used here.
//      region.reference_id = ref_id;
//      region.region_index = 0;
//      region.region_votes = (end_hit - begin_hit + 1);
//      region.is_split = false;
//      region.split_start = 0;
//      region.split_end = 0;
//
//      begin_hit = end_hit + 1;
//
//      if (region.reference_id == 0xFFFFFFFFFFFFFFFF ||
//          region.region_votes < 5) {
//        continue;
//      }
//      regions.emplace_back(region);
    }
  }

  std::sort(owler_data.overlaps.begin(), owler_data.overlaps.end(), [](const PairwiseOverlap& o1, const PairwiseOverlap& o2) { return o1.num_seeds > o2.num_seeds; });

  return 0;
}

void Owler::AppendSeedHits_(const uint128_t& seed, std::shared_ptr<is::MinimizerIndex> index, bool threshold_hits, double count_cutoff, bool is_overlapper, int64_t qid, std::vector<uint128_t> &all_hits) {
  int64_t key = is::MinimizerIndex::seed_key(seed);
  int32_t pos_read_int32 = is::MinimizerIndex::seed_position(seed);
  uint128_t pos_read = pos_read_int32;

  const uint128_t *found_seeds = NULL;
  int64_t num_found_seeds = 0;

  int lookup_ret = index->KeyLookup(key, &found_seeds, &num_found_seeds);

  // Something went wrong.
  if (lookup_ret) {
    return;
  }

  // Filter seeds with excessive hits if necessary.
  if (threshold_hits &&
      num_found_seeds >= count_cutoff) {
    return;
  }

  all_hits.insert(all_hits.end(), num_found_seeds, 0);

  // Reformat the new hits. The key part is no longer needed,
  // but source position on read is.
  auto *phits = &all_hits[0] + (all_hits.size() - num_found_seeds);
  int64_t k_added = 0;

  for (int64_t k=0; k<num_found_seeds; k++) {
    if (found_seeds[k] == kInvalidSeed) {
      continue;
    }

    int32_t seq_id_int32 = is::MinimizerIndex::seed_seq_id(found_seeds[k]);
    uint128_t seq_id = seq_id_int32;
    if (is_overlapper && seq_id_int32 <= qid) {
      continue;
    }

    int32_t pos_ref_int32 = is::MinimizerIndex::seed_position(found_seeds[k]);
    uint128_t diag = ((uint128_t) (pos_ref_int32 - pos_read_int32)) & kSeedMask32_1;      // Diagonals can be negative values.
    uint128_t pos_ref = pos_ref_int32;

    phits[k_added++] = MakeHit_(seq_id, diag, pos_ref, pos_read);
  }

  all_hits.resize(all_hits.size() - num_found_seeds + k_added);
}

inline uint128_t Owler::MakeHit_(const uint128_t& seq_id, const uint128_t& diag, const uint128_t& pos_ref, const uint128_t& pos_read) {
  return ((seq_id << 96) | (diag << 64) | (pos_ref << 32) | (pos_read));
}

void Owler::GenerateOutput_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data) {
}

inline int32_t Owler::HitPosRead_(const uint128_t& hit) {
  return (int32_t) (hit & kSeedMask32_1);
}

inline int32_t Owler::HitPosRef_(const uint128_t& hit) {
  return (int32_t) ((hit & kSeedMask32_2) >> 32);
}

inline int32_t Owler::HitDiag_(const uint128_t& hit) {
  return (int32_t) ((hit & kSeedMask32_3) >> 64);
}

inline int32_t Owler::HitSeqId_(const uint128_t& hit) {
  return (int32_t) ((hit & kSeedMask32_4) >> 96);
}

void Owler::FindRegionBounds_(const std::vector<uint128_t> &all_hits, int64_t begin_hit, int64_t end_hit, int64_t &start_read, int64_t &end_read, int64_t &start_ref, int64_t &end_ref) {
  if ((end_hit - begin_hit + 1) <= 0) {
    return;
  }

  start_ref = HitPosRef_(all_hits[begin_hit]);
  end_ref = start_ref;
  for (int64_t i=begin_hit; i<end_hit; i++) {
    int64_t curr_pos = HitPosRef_(all_hits[i]);
    start_ref = std::min(start_ref, curr_pos);
    end_ref = std::max(end_ref, curr_pos);
  }
}









// Assumes vertices are sorted.
// If use_l1_filtering is true, then all vertices/anchors that have coordinates further than allowed_dist from the L1 line are filtered out.
// Otherwise, all vertices will be used.
// The L1 line is specified with k = 1 and l parameters (y = k*x + l).
//void Owler::CalcLCSFromLocalScoresCacheFriendly2_(OwlerData* owler_data, int64_t ref_hits_start, int64_t ref_hits_end, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices) {
int Owler::LCSkFilter_(const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int64_t qid, int64_t tid, int32_t seed_len, PairwiseOverlap &overlap) {
  std::vector<uint128_t> events;
  std::vector<uint64_t> matches_starts;
  std::vector<uint64_t> matches_indices;
  int64_t max_seq_len = 0;

//  printf ("(1) LCSk:\n");
//  for (int32_t i=0; i<hits.size(); i++) {
//    int32_t p = i;
//    printf ("HitPosRead_[%d] = %d, HitPosRef_[%d] = %d, p = %d, hits[p] = %s\n", i, HitPosRead_(hits[p]), i, HitPosRef_(hits[p]), p, is::PrintSeed(hits[p]).c_str());
//  }
//  printf ("\n");
//  fflush(stdout);

  PrepareEvents_(hits, begin_hit, end_hit, seed_len, events, matches_starts, matches_indices, max_seq_len);

  int64_t lcsk_len = 0;
  std::vector<int32_t> lcsk_indices;
  LCSk_(events, max_seq_len, seed_len, matches_starts, matches_indices, lcsk_indices, lcsk_len);

  lcsk_len = std::max((int64_t) 0, lcsk_len - seed_len);     // Reduce the length by the size of a seed, because the end point
                                                   // Also do not include the last seed.

  if (lcsk_indices.size() == 0) {
    return 1;
  }

//  printf ("(2) LCSk:\n");
//  for (int32_t i=begin_hit; i<(end_hit-begin_hit); i++) {
//    int32_t p = i;
//    printf ("HitPosRead_[%d] = %d, HitPosRef_[%d] = %d, p = %d, hits[p] = %s\n", i, HitPosRead_(hits[p]), i, HitPosRef_(hits[p]), p, is::PrintSeed(hits[p]).c_str());
//  }
//  printf ("\n");
//  fflush(stdout);
//
//  printf ("(3) LCSk indices (size = %ld):\n", lcsk_indices.size());
//  for (int32_t i=0; i<lcsk_indices.size(); i++) {
//    int32_t p = begin_hit + lcsk_indices[i];
////    p = i;
//    printf ("HitPosRead_[%d] = %d, HitPosRef_[%d] = %d, p = %d, hits[p] = %s\n", i, HitPosRead_(hits[p]), i, HitPosRef_(hits[p]), p, is::PrintSeed(hits[p]).c_str());
//  }
//  printf ("\n");
//  printf ("qid = %ld, tid = %ld\n", qid, tid);
//  fflush(stdout);
//
//  if (tid > 0)
//    exit(1);

  int64_t front = begin_hit + lcsk_indices.front();
  int64_t back = begin_hit + lcsk_indices.back();

  overlap.query.start = HitPosRead_(hits[back]);
  overlap.query.end = HitPosRead_(hits[front]) + 1; // +1 means that the end is not inclusive. // + seed_len;  TODO: The end does not cover the last seed entirely!
  overlap.target.start = HitPosRef_(hits[back]);
  overlap.target.end = HitPosRef_(hits[front]) + 1; // + seed_len;
  overlap.num_seeds = lcsk_indices.size();
  overlap.cov_bases = lcsk_len;
  overlap.qid = qid;
  overlap.tid = tid;

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
    int primary_diagonal = n - 1 + i - j;

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

std::string Owler::GenerateMHAPLine(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap) {
  std::stringstream ss;

  float jaccard_score = std::min(1.0f, ((float) overlap.cov_bases) / ((float) (overlap.target.end - overlap.target.start)));
  int64_t shared_minmers = overlap.num_seeds;

  bool read_is_reverse = overlap.tid > index->get_num_sequences_forward();
  bool ref_is_reverse = false;

  int64_t read_id = overlap.qid;
  int64_t read_start = overlap.query.start;
  int64_t read_end = overlap.query.end;
  int64_t read_length = read->get_sequence_length();
  if (read_is_reverse) {
    read_start = read_length - overlap.query.end - 1;
    read_end = read_length - overlap.query.start;
  }

  int64_t ref_id = overlap.tid % index->get_num_sequences_forward();
  int64_t ref_start = overlap.target.start;
  int64_t ref_end = overlap.target.end;
  int64_t ref_length = index->get_reference_lengths()[overlap.tid];

  ss << (read_id + 1) << " ";      /// read1_id
  ss << (ref_id + 1) << " ";      /// read2_id
  ss << jaccard_score << " ";      /// Jaccard score
  ss << shared_minmers << " ";        /// Shared minmers
  ss << (read_is_reverse ? 1 : 0) << " ";  /// A is reverse
  ss << read_start << " ";
  ss << read_end << " ";
  ss << read_length << " ";
  ss << (ref_is_reverse ? 1 : 0) << " ";
  ss << ref_start << " ";
  ss << ref_end << " ";
  ss << ref_length;

//  ss << "\tcov_bases = " << overlap.cov_bases;

  return ss.str();
}

int Owler::CheckOverlap_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap) {
  if (overlap.num_seeds == 0) {
    return 1;
  }

  if (overlap.cov_bases < 50) {
    return 1;
  }

  return 0;
}
