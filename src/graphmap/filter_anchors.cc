/*
 * filter.cc
 *
 *  Created on: Mar 15, 2016
 *      Author: isovic
 */

#include "filter_anchors.h"
#include "log_system/log_system.h"

//#define DEBUG_TEST1

#define Dist(x1, y1, x2, y2) sqrt((double) (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))

int64_t CalcScore(int32_t qpos, int32_t rpos, int32_t next_qpos, int32_t next_rpos, double indel_bandwidth_margin, int32_t fwd_length, double *score_gap, double *score_dist) {
//  int32_t min_dist = (next_qpos <= qpos) ? (qpos - next_qpos) :
//                     (next_qpos >= (qpos + seed_length_q)) ? (next_qpos - qpos - seed_length_q) :
//                     -1;

  int32_t l1 = rpos - qpos, l2 = next_rpos - next_qpos;
  double a = abs(l2 - l1);      // Distance between diagonals. This is the actual gap score should we extend the diagonal of the (qpos, rpos) to the same X coordinate of the next point (next_qpos).
  double band = a * (sqrt(2.0) / 2.0);
//  double dist = sqrt((double) (rpos - next_rpos) * (rpos - next_rpos) + (qpos - next_qpos) * (qpos - next_qpos));
  double dist = Dist(qpos, rpos, next_qpos, next_rpos);
  double mismatch = Dist(qpos, rpos, next_qpos, (next_rpos - (l2 - l1)));     // Calculate the length of the diagonal until we hit the next_rpos coordinate.

  *score_gap = -a;
  *score_dist = -mismatch / 2.0;     // Divide by two because not necessarily all bases on the diagonal are mismatches, some might be matches. Just a heuristic.

  int32_t bandwidth = indel_bandwidth_margin * dist;
//  int32_t bandwidth = indel_bandwidth_margin * fwd_length;
//  printf ("  -> fwd_length = %d, dist = %f, mismatch = %f, band = %f\n", fwd_length, dist, mismatch, band);
//  fflush(stdout);

  int32_t x_dist = abs(qpos - next_qpos);

  if ((fwd_length <= 0 || x_dist <= fwd_length) && band <= bandwidth) {
    // All is good.
    return 0;

  } else if ((fwd_length > 0 && x_dist > fwd_length) && band <= bandwidth) {
    // Bandwidth is fine but distance is larger than expected. Perhaps consider this one and then break the chain?
    return 1;

  } else if ((fwd_length <= 0 || x_dist <= fwd_length) && band > bandwidth) {
    // Distance is good, but it's out of bandwidth. Break the chain.
    return 2;

  } else if ((fwd_length > 0 && x_dist > fwd_length) && band > bandwidth) {
//    printf ("band = %f, bandwidth = %d\n", band, bandwidth);
//    fflush(stdout);
    // All numbers are too high, break the chain.
    return 3;
  }

  return 4;
}

// Used for packed 128-bit hits, mainly for Owler.
void GetPositionsFrom128bit(const std::vector<uint128_t> &hits, const std::vector<int> &lcskpp_indices, int64_t lcskpp_id, int32_t seed_len, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end) {
  *qpos_start = get128_qpos(hits[lcskpp_indices[lcskpp_id]]);
  *rpos_start = get128_rpos(hits[lcskpp_indices[lcskpp_id]]);
  *qpos_end = *qpos_start + seed_len;
  *rpos_end = *rpos_start + seed_len;
}

// Used for the vertices in GraphMap. The vertices are the ugly representation of several arrays which makes it more cache friendly.
void GetPositionsFromRegistry(const Vertices& registry_entries, const std::vector<int> &lcskpp_indices, int64_t lcskpp_id, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end) {
  int32_t vertex_id = lcskpp_indices[lcskpp_id];
  *qpos_start = registry_entries.query_starts[vertex_id];
  *rpos_start = registry_entries.reference_starts[vertex_id];
  *qpos_end = registry_entries.query_ends[vertex_id];
  *rpos_end = registry_entries.reference_ends[vertex_id];
}

// Used for the vertices in GraphMap. The vertices are the ugly representation of several arrays which makes it more cache friendly.
void GetPositionsFromRegistry2(const Vertices& registry_entries, int64_t vertex_id, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end) {
  *qpos_start = registry_entries.query_starts[vertex_id];
  *rpos_start = registry_entries.reference_starts[vertex_id];
  *qpos_end = registry_entries.query_ends[vertex_id];
  *rpos_end = registry_entries.reference_ends[vertex_id];
}

//int FilterAnchors1(const SingleSequence* seq, ScoreRegistry* local_score, const ProgramParameters *parameters,
//                  const std::vector<int> &lcskpp_indices, double indel_bandwidth_margin, int32_t max_dist, int32_t lookahead_dist_factor, int64_t min_covered_bases, int32_t cluster_size_cutoff,
//                  std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {
//
////  int32_t bandwidth = 10;
//  //  int32_t max_dist = 0.10f * seq->get_sequence_length();
//
////  double bandwidth_margin = 0.10f;
////  int32_t max_dist = 200;
////  int32_t lookahead_dist_factor = 2;
//  int64_t chain_len_lookahead = max_dist;
//
//  const Vertices& registry_entries = local_score->get_registry_entries();
//
//  std::vector<std::vector<int32_t> > chains;
//  std::vector<int32_t> chain_cov_bases;
//
//  {
//    std::vector<int32_t> new_chain;
//    new_chain.reserve(lcskpp_indices.size());
//    chains.push_back(new_chain);
//    chain_cov_bases.push_back(0);
//  }
//  for (int64_t i=0; i<lcskpp_indices.size(); i++) {
//    int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
//    GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
//    int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));
////    int32_t qpos = get128_qpos(hits[lcskpp_indices[i]]);
////    int32_t rpos = get128_rpos(hits[lcskpp_indices[i]]);
//
//#ifdef DEBUG_TEST1
//    printf ("[i = %ld] Taking the next seed. qpos = %d, rpos = %d\n", i, qpos_start, rpos_start);
//    fflush(stdout);
////    if (chains.back().size() > 0) {
////      printf ("Chain info:\n");
////      printf ("  get128_qpos(hits[chains.back().front() = %d, get128_qpos(hits[chains.back().back()]) = %d\n", get128_qpos(hits[chains.back().front()]), get128_qpos(hits[chains.back().back()]));
////      printf ("  get128_rpos(hits[chains.back().front() = %d, get128_rpos(hits[chains.back().back()]) = %d\n", get128_rpos(hits[chains.back().front()]), get128_rpos(hits[chains.back().back()]));
////    }
//#endif
//
//    double max_score = 1;
//    int64_t max_score_id = -1;
//    int64_t next_id = (i + 1);
//    while (next_id < lcskpp_indices.size()) {
////      int32_t next_qpos = get128_qpos(hits[lcskpp_indices[next_id]]);
////      int32_t next_rpos = get128_rpos(hits[lcskpp_indices[next_id]]);
//      int32_t next_qpos_start = 0, next_rpos_start = 0, next_qpos_end = 0, next_rpos_end = 0;
//      GetPositionsFromRegistry(registry_entries, lcskpp_indices, next_id, &next_qpos_start, &next_rpos_start, &next_qpos_end, &next_rpos_end);
//      int32_t next_cov_bases = std::max((next_qpos_end - next_qpos_start), (next_rpos_end - next_rpos_start));
//
//      #ifdef DEBUG_TEST1
//            printf ("  [next_id = %ld] Considering next seed for chaining. qpos = %d, rpos = %d\n", next_id, next_qpos_start, next_rpos_start);
//            fflush(stdout);
//      #endif
//      double score_gap = 0.0, score_mismatch = 0.0;
//      int64_t ret_val = CalcScore(qpos_start, rpos_start, next_qpos_end, next_rpos_end, indel_bandwidth_margin, chain_len_lookahead, &score_gap, &score_mismatch);
//      double score = score_gap; // + score_mismatch;
//      #ifdef DEBUG_TEST1
//            if (ret_val > 0) { printf ("  Break.\n"); break; }
//      #endif
//      if (ret_val > 0) { break; }
//      if (max_score > 0 || score > max_score) {
//        #ifdef DEBUG_TEST1
//                printf ("  New max score found. Previously: score = %d, score_id = %d. Now: score = %d, score_id = %d\n", max_score, max_score_id, score, next_id);
//                fflush(stdout);
//        #endif
//        max_score = score;
//        max_score_id = next_id;
//      }
//      next_id += 1;
//    }
//    if (max_score_id == -1) {
//      chains[chains.size()-1].push_back(lcskpp_indices[i]);
//      chain_cov_bases.back() += cov_bases;
//
//      #ifdef DEBUG_TEST1
//            printf ("  -> Num chains: %ld, last chain size: %ld\n", chains.size(), chains.back().size());
//            fflush(stdout);
//      #endif
//      std::vector<int32_t> new_chain;
//      new_chain.reserve(lcskpp_indices.size() - i);
//      chains.push_back(new_chain);
//      chain_cov_bases.push_back(0);
////      i = max_score_id - 1;
//      #ifdef DEBUG_TEST1
//            printf ("  (I) Going for next i. i = %ld\n", i);
//            printf ("\n");
//            fflush(stdout);
//      #endif
//      continue;
//    }
//    int32_t max_qpos_start = 0, max_rpos_start = 0, max_qpos_end = 0, max_rpos_end = 0;
//    GetPositionsFromRegistry(registry_entries, lcskpp_indices, max_score_id, &max_qpos_start, &max_rpos_start, &max_qpos_end, &max_rpos_end);
//    int32_t max_cov_bases = std::max((max_qpos_end - max_qpos_start), (max_rpos_end - max_rpos_start));
//    chains[chains.size()-1].push_back(lcskpp_indices[max_score_id]);
//    chain_cov_bases.back() += max_cov_bases;
//    i = max_score_id - 1;
//    #ifdef DEBUG_TEST1
//        printf ("  -> Num chains: %ld, last chain size: %ld\n", chains.size(), chains.back().size());
//        printf ("  (II) Going for next i. i = %ld\n", i);
//        printf ("\n");
//        fflush(stdout);
//    #endif
//  }
//
//#ifdef DEBUG_TEST1
//  for (int64_t i=0; i<chains.size(); i++) {
//    printf ("Chain %ld, size: %ld\n", i, chains[i].size());
//    fflush(stdout);
//    for (int64_t j=0; j<chains[i].size(); j++) {
//      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
//      GetPositionsFromRegistry2(registry_entries, chains[i][j], &qpos_start, &rpos_start, &qpos_end, &rpos_end);
//      printf ("  [%ld] %d %d\n", j, qpos_start, rpos_start);
//      fflush(stdout);
//    }
//    printf ("\n");
//    fflush(stdout);
//  }
//#endif
//
//  ret_filtered_lcskpp_indices.size();
//  if (ret_cluster_ids) {
//    ret_cluster_ids->clear();
//  }
//
//  int64_t num_filtered_indices = 0;
//  for (int64_t i=0; i<chains.size(); i++) { num_filtered_indices += chains[i].size(); }
//  ret_filtered_lcskpp_indices.reserve(num_filtered_indices);
//  if (ret_cluster_ids) { ret_cluster_ids->reserve(num_filtered_indices); }
//
////  // Attempt to join neighbouring clusters if they are on the same (similar) diagonal.
////  printf ("Attempting join neighbouring chains.\n");
////  fflush(stdout);
////  int64_t prev_chain = 0;
////  for (int64_t i=1; i<chains.size(); i++) {
////    if (chains[i].size() == 0) {
////      continue;
////    }
////    printf ("i = %ld / %ld\n", (i + 1), chains.size());
////    printf ("  prev_chain = %ld, prev_chain.length() = %ld\n", prev_chain, chains[prev_chain].size());
////    fflush(stdout);
////    int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
////    GetPositionsFromRegistry2(registry_entries, chains[i].back(), &qpos_start, &rpos_start, &qpos_end, &rpos_end);
////
////    int32_t prev_qpos_start = 0, prev_rpos_start = 0, prev_qpos_end = 0, prev_rpos_end = 0;
////    GetPositionsFromRegistry2(registry_entries, chains[prev_chain].back(), &prev_qpos_start, &prev_rpos_start, &prev_qpos_end, &prev_rpos_end);
////
////    double score_gap = 0.0, score_mismatch = 0.0;
////    int64_t ret_val = CalcScore(prev_qpos_end, prev_rpos_end, qpos_start, rpos_start, 0.15f, -1, &score_gap, &score_mismatch);
////
////    if (ret_val == 0) {
////      chains[prev_chain].insert(chains[prev_chain].end(), chains[i].begin(), chains[i].end());
////      chains[i].clear();
////    } else {
////      prev_chain = i;
////    }
//////    printf ("[%ld] (front) %d %d, (back) %d %d\n", i, qpos_start, rpos_start, prev_qpos_start, prev_rpos_start);
//////    fflush(stdout);
////  }
////  exit(1);
//
//  // Remove chains with too low of a count.
//  for (int64_t i=0; i<chains.size(); i++) {
//    if (chains[i].size() <= cluster_size_cutoff && chain_cov_bases[i] < min_covered_bases) { continue; }
//
////    if (chain_cov_bases[i] < min_covered_bases) { continue; }
////    for (int64_t j=(chains[i].size()-1); j>=0; j--) {
//    for (int64_t j=0; j<chains[i].size(); j++) {
//      ret_filtered_lcskpp_indices.push_back(chains[i][j]);
//      if (ret_cluster_ids) {
//        ret_cluster_ids->push_back(chains.size() - i);
//      }
//    }
//  }
//
//#ifdef DEBUG_TEST1
//  printf ("\n");
//  printf ("Chain sizes:\n");
//  for (int64_t i=0; i<chains.size(); i++) {
//    printf ("%ld ", chains[i].size());
//  }
//  printf ("\n");
////  for (int64_t i=0; i<chain_sizes.size(); i++) {
////    printf ("%ld ", chain_sizes[i]);
////  }
//  printf ("\n");
//#endif
//
//#ifdef DEBUG_TEST1
//  printf ("ret_filtered_lcskpp_indices.size() = %ld\n", ret_filtered_lcskpp_indices.size());
//#endif
//
//  return 0;
//
//}






double CalcAnchorDiff(double qpos_start, double qpos_end, double rpos_start, double rpos_end) {
  double r_diff = rpos_end - rpos_start;
  double q_diff = qpos_end - qpos_start;
  double perc_diff = std::min(q_diff, r_diff) / std::max(q_diff, r_diff);
  return perc_diff;
}

void CalcAnchorLengthStatistics(const ScoreRegistry* local_score, const std::vector<int> &lcskpp_indices, double *ret_avg, double *ret_std_dev) {
  if (lcskpp_indices.size() == 0) { return; }

  const Vertices& registry_entries = local_score->get_registry_entries();

  // Calculate the average.
  double avg = 0.0f;
  for (int64_t i=0; i<lcskpp_indices.size(); i++) {
    int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
    GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
    avg += CalcAnchorDiff(qpos_start, qpos_end, rpos_start, rpos_end);
  }
  avg /= ((double) lcskpp_indices.size());

  // Calculate the standard deviation.
  double std_dev = 0.0f;
  for (int64_t i=0; i<lcskpp_indices.size(); i++) {
    int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
    GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
    double perc_diff = CalcAnchorDiff(qpos_start, qpos_end, rpos_start, rpos_end);
    std_dev += (perc_diff - avg) * (perc_diff - avg);
  }
  std_dev /= ((double) lcskpp_indices.size());
  std_dev = sqrt(std_dev);

  *ret_avg = avg;
  *ret_std_dev = std_dev;
}

int FilterAnchorsByDiff(const SingleSequence* read, ScoreRegistry* local_score, const ProgramParameters *parameters,
                  const std::vector<int> &lcskpp_indices, std::vector<int> &ret_filtered_lcskpp_indices) {
  double avg = 0.0f, std_dev = 0.0f;
  CalcAnchorLengthStatistics(local_score, lcskpp_indices, &avg, &std_dev);
  double min_diff = avg - 3*std_dev;

  const Vertices& registry_entries = local_score->get_registry_entries();

  ret_filtered_lcskpp_indices.clear();
  ret_filtered_lcskpp_indices.reserve(lcskpp_indices.size());

  for (int64_t i=0; i<lcskpp_indices.size(); i++) {
    int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
    GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
    double perc_diff = CalcAnchorDiff(qpos_start, qpos_end, rpos_start, rpos_end);
    if (perc_diff >= min_diff) {
      ret_filtered_lcskpp_indices.push_back(lcskpp_indices[i]);
    }
  }

  LOG_DEBUG_SPEC("Filtered %ld anchors.\n", ret_filtered_lcskpp_indices.size() - lcskpp_indices.size());

  return 0;
}

int FilterAnchorsByChaining(const SingleSequence* read, ScoreRegistry* local_score, const ProgramParameters *parameters,
                  const std::vector<int> &lcskpp_indices, double indel_bandwidth_margin, int32_t max_dist, int32_t lookahead_dist_factor, int64_t min_covered_bases, int32_t cluster_size_cutoff,
                  std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {

  int64_t chain_len_lookahead = max_dist;
  const Vertices& registry_entries = local_score->get_registry_entries();

  std::vector<std::vector<int32_t> > chains;
  std::vector<int32_t> chain_cov_bases;

  {
    std::vector<int32_t> new_chain;
    new_chain.reserve(lcskpp_indices.size());
    chains.push_back(new_chain);
    chain_cov_bases.push_back(0);
  }

  double avg_diff = 0.0f, std_diff = 0.0f;
  CalcAnchorLengthStatistics(local_score, lcskpp_indices, &avg_diff, &std_diff);
  double min_diff = avg_diff - 3*std_diff;

  LOG_DEBUG_SPEC_NO_HEADER("Starting a new chain: %ld\n", chains.size());
//  fflush(stdout);
  for (int64_t i=0; i<lcskpp_indices.size(); i++) {
    double max_score_gap = 1.0, max_score_mismatch = 1.0;
    int64_t max_score_id = -1;
    int64_t next_id = i;
    while (next_id < lcskpp_indices.size()) {
      if (chains.back().size() == 0) {
        int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
        GetPositionsFromRegistry(registry_entries, lcskpp_indices, next_id, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
        int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));
        chains.back().push_back(lcskpp_indices[next_id]);
        chain_cov_bases.back() += cov_bases;
        LOG_DEBUG_SPEC_NO_HEADER("  (I) Adding ancho %ld (chains.back().size() == 0): start = [%d, %d], end = [%d, %d]\n", 0, qpos_start, rpos_start, qpos_end, rpos_end);
//        fflush(stdout);
//        next_id += 1;
        continue;
      }

      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
      GetPositionsFromRegistry2(registry_entries, chains.back().back(), &qpos_start, &rpos_start, &qpos_end, &rpos_end);
      int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));

      int32_t next_qpos_start = 0, next_rpos_start = 0, next_qpos_end = 0, next_rpos_end = 0;
      GetPositionsFromRegistry(registry_entries, lcskpp_indices, next_id, &next_qpos_start, &next_rpos_start, &next_qpos_end, &next_rpos_end);
      int32_t next_cov_bases = std::max((next_qpos_end - next_qpos_start), (next_rpos_end - next_rpos_start));

      LOG_DEBUG_SPEC_NO_HEADER("  Considering anchor %ld (i = %ld, lcskpp_indices[%ld] = %d, chains.back().back() = %d): start = [%d, %d], end = [%d, %d]. max_score_id = %ld\n", next_id, i, next_id, lcskpp_indices[next_id], chains.back().back(), next_qpos_start, next_rpos_start, next_qpos_end, next_rpos_end, max_score_id);

      double score_gap = 0.0, score_mismatch = 0.0;
      int64_t ret_val = CalcScore(qpos_start, rpos_start, next_qpos_end, next_rpos_end, indel_bandwidth_margin, chain_len_lookahead, &score_gap, &score_mismatch);

      if (ret_val > 0) {
        LOG_DEBUG_SPEC_NO_HEADER("    - Breaking. CalcScore returned a value > 0 (ret_val = %ld)\n", ret_val);
        break;
      }

      if (max_score_id == -1 || (score_gap >= max_score_gap && score_mismatch >= max_score_mismatch)) {
        LOG_DEBUG_SPEC_NO_HEADER("    - This will be the new max_score_id.\n", ret_val);
        max_score_gap = score_gap;
        max_score_mismatch = score_mismatch;
        max_score_id = next_id;
      }
      next_id += 1;
    }

    if (max_score_id == -1) { // If nothing was found, just add the current lcskpp_indices[i].
      std::vector<int32_t> new_chain;
      chains.push_back(new_chain);
      chains.back().reserve(lcskpp_indices.size() - i);
      chain_cov_bases.push_back(0);

      LOG_DEBUG_SPEC_NO_HEADER("  End of chain.\n\n");
      LOG_DEBUG_SPEC_NO_HEADER("  Creating a new chain (new chains.size() = %ld\n", chains.size());

      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
      GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
      int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));
      chains.back().push_back(lcskpp_indices[i]);
      chain_cov_bases.back() += cov_bases;

      LOG_DEBUG_SPEC_NO_HEADER("  (II) Adding anchor %ld (max_score_id == -1): start = [%d, %d], end = [%d, %d]. (Case where max_score_id == -1)\n", i, qpos_start, rpos_start, qpos_end, rpos_end);
//      fflush(stdout);
//      LOG_DEBUG_SPEC_NO_HEADER("Starting a new chain: %ld\n", chains.size());
//      fflush(stdout);

    } else {
//    } else if (chains.back().back() != lcskpp_indices[max_score_id]) {  // In this case, there was a max_score found, but we also check if the lcskpp index was already added to the array.
      int32_t max_qpos_start = 0, max_rpos_start = 0, max_qpos_end = 0, max_rpos_end = 0;
      GetPositionsFromRegistry(registry_entries, lcskpp_indices, max_score_id, &max_qpos_start, &max_rpos_start, &max_qpos_end, &max_rpos_end);
      int32_t max_cov_bases = std::max((max_qpos_end - max_qpos_start), (max_rpos_end - max_rpos_start));

      chains.back().push_back(lcskpp_indices[max_score_id]);
      chain_cov_bases.back() += max_cov_bases;
      i = max_score_id;

      LOG_DEBUG_SPEC_NO_HEADER("  (III) Adding anchor (max_score_id = %ld), start = [%d, %d], end = [%d, %d]\n", max_score_id, max_qpos_start, max_rpos_start, max_qpos_end, max_rpos_end);
//      fflush(stdout);
    }
  }

  ret_filtered_lcskpp_indices.size();
  if (ret_cluster_ids) { ret_cluster_ids->clear(); }



  LOG_DEBUG_SPEC_NO_HEADER("Filtering short chains.\n");
  for (int64_t i=0; i<chains.size(); i++) {
    // Remove chains with too low of a count.
    if (chains[i].size() == 0) { continue; }

    int32_t qpos_start1 = 0, rpos_start1 = 0, qpos_end1 = 0, rpos_end1 = 0;
    GetPositionsFromRegistry2(registry_entries, chains[i].front(), &qpos_start1, &rpos_start1, &qpos_end1, &rpos_end1);
    int32_t qpos_start2 = 0, rpos_start2 = 0, qpos_end2 = 0, rpos_end2 = 0;
    GetPositionsFromRegistry2(registry_entries, chains[i].back(), &qpos_start2, &rpos_start2, &qpos_end2, &rpos_end2);

//    printf ("chains[%ld].size() = %ld, chain_cov_bases[%ld] = %ld, qlen = %d, rlen = %d, qpos_start1 = %d, qpos_end2 = %d;   ", i, chains[i].size(), i, chain_cov_bases[i], (qpos_end2 - qpos_start1), (rpos_end2 - rpos_start1), qpos_start1, qpos_end2);
    if (chains[i].size() <= cluster_size_cutoff && chain_cov_bases[i] < min_covered_bases) {
      LOG_DEBUG_SPEC_NO_HEADER("Removed!");
      chains[i].clear();
    }
    LOG_DEBUG_SPEC_NO_HEADER("\n");
  }

  LOG_DEBUG_SPEC_NO_HEADER("cluster_size_cutoff = %ld\n", cluster_size_cutoff);
  LOG_DEBUG_SPEC_NO_HEADER("min_covered_bases = %ld\n", min_covered_bases);

  int64_t num_filtered_indices = 0;
  for (int64_t i=0; i<chains.size(); i++) { num_filtered_indices += chains[i].size(); }
  ret_filtered_lcskpp_indices.reserve(num_filtered_indices);
  if (ret_cluster_ids) { ret_cluster_ids->reserve(num_filtered_indices); }

  // Fill the return indices.
  for (int64_t i=0; i<chains.size(); i++) {
    // Remove chains with too low of a count.
    if (chains[i].size() <= cluster_size_cutoff && chain_cov_bases[i] < min_covered_bases) { continue; }

    int32_t qpos_start1 = 0, rpos_start1 = 0, qpos_end1 = 0, rpos_end1 = 0;
    GetPositionsFromRegistry2(registry_entries, chains[i].front(), &qpos_start1, &rpos_start1, &qpos_end1, &rpos_end1);
    int32_t qpos_start2 = 0, rpos_start2 = 0, qpos_end2 = 0, rpos_end2 = 0;
    GetPositionsFromRegistry2(registry_entries, chains[i].back(), &qpos_start2, &rpos_start2, &qpos_end2, &rpos_end2);

    LOG_DEBUG_SPEC_NO_HEADER("chains[%ld].size() = %ld, chain_cov_bases[%ld] = %ld, qlen = %d, rlen = %d, qpos_start1 = %d, qpos_end2 = %d;   ", i, chains[i].size(), i, chain_cov_bases[i], (qpos_end2 - qpos_start1), (rpos_end2 - rpos_start1), qpos_start1, qpos_end2);
    LOG_DEBUG_SPEC_NO_HEADER("\n");

    for (int64_t j=0; j<chains[i].size(); j++) {
      ret_filtered_lcskpp_indices.push_back(chains[i][j]);
      if (ret_cluster_ids) {
        ret_cluster_ids->push_back(chains.size() - i);
      }
    }
  }

//  fflush(stdout);
//  exit(1);

  return 0;

}
