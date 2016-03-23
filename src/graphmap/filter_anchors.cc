/*
 * filter.cc
 *
 *  Created on: Mar 15, 2016
 *      Author: isovic
 */

#include "filter_anchors.h"

//#define DEBUG_TEST1

int64_t CalcScore(int32_t qpos, int32_t rpos, int32_t next_qpos, int32_t next_rpos, double indel_bandwidth_margin, int32_t fwd_length) {
//  int32_t min_dist = (next_qpos <= qpos) ? (qpos - next_qpos) :
//                     (next_qpos >= (qpos + seed_length_q)) ? (next_qpos - qpos - seed_length_q) :
//                     -1;
  int32_t min_dist = (next_qpos <= qpos) ? (qpos - next_qpos) : (next_qpos - qpos);

  if (min_dist < 0 || min_dist > fwd_length) {
#ifdef DEBUG_TEST1
  printf ("  --[min_dist not satisfied!] min_dist = %ld, fwd_length = %ld\n", min_dist, fwd_length);
  fflush(stdout);
#endif
    return 1;
  }

  int32_t l1 = rpos - qpos, l2 = next_rpos - next_qpos;
  double a = abs(l2 - l1);
  double band = a * (sqrt(2.0) / 2.0);
  double dist = sqrt((double) (rpos - next_rpos) * (rpos - next_rpos) + (qpos - next_qpos) * (qpos - next_qpos));
  int32_t bandwidth = indel_bandwidth_margin * dist;
#ifdef DEBUG_TEST1
  printf ("  bandwidth = %ld, band = %f, dist = %f, a = %f\n", bandwidth, band, dist, a);
  fflush(stdout);
#endif
  if (band > bandwidth) {
#ifdef DEBUG_TEST1
  printf ("  --[bandwidth not satisfied!]\n", bandwidth, band, dist, a);
  fflush(stdout);
#endif
    return 1;
  }     // Do I need to check the diagonal, or the actual gap? 'a' is the vertical/horizontal distance, and 'band' is the diagonal between two lines.

  return -a;
}

// Used for packed 128-bit hits, mainly for Owler.
inline void GetPositionsFrom128bit(const std::vector<uint128_t> &hits, const std::vector<int> &lcskpp_indices, int64_t lcskpp_id, int32_t seed_len, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end) {
  *qpos_start = get128_qpos(hits[lcskpp_indices[lcskpp_id]]);
  *rpos_start = get128_rpos(hits[lcskpp_indices[lcskpp_id]]);
  *qpos_end = *qpos_start + seed_len;
  *rpos_end = *rpos_start + seed_len;
}

// Used for the vertices in GraphMap. The vertices are the ugly representation of several arrays which makes it more cache friendly.
inline void GetPositionsFromRegistry(const Vertices& registry_entries, const std::vector<int> &lcskpp_indices, int64_t lcskpp_id, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end) {
  int32_t vertex_id = lcskpp_indices[lcskpp_id];
  *qpos_start = registry_entries.query_starts[vertex_id];
  *rpos_start = registry_entries.reference_starts[vertex_id];
  *qpos_end = registry_entries.query_ends[vertex_id];
  *rpos_end = registry_entries.reference_ends[vertex_id];
}

// Used for the vertices in GraphMap. The vertices are the ugly representation of several arrays which makes it more cache friendly.
inline void GetPositionsFromRegistry2(const Vertices& registry_entries, int64_t vertex_id, int32_t *qpos_start, int32_t *rpos_start, int32_t *qpos_end, int32_t *rpos_end) {
  *qpos_start = registry_entries.query_starts[vertex_id];
  *rpos_start = registry_entries.reference_starts[vertex_id];
  *qpos_end = registry_entries.query_ends[vertex_id];
  *rpos_end = registry_entries.reference_ends[vertex_id];
}

int FilterAnchors(const SingleSequence* seq, ScoreRegistry* local_score, const ProgramParameters *parameters,
                  const std::vector<int> &lcskpp_indices, double indel_bandwidth_margin, int32_t max_dist, int32_t lookahead_dist_factor, int64_t min_covered_bases, int32_t cluster_size_cutoff,
                  std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {

//  int32_t bandwidth = 10;
  //  int32_t max_dist = 0.10f * seq->get_sequence_length();

//  double bandwidth_margin = 0.10f;
//  int32_t max_dist = 200;
//  int32_t lookahead_dist_factor = 2;

  const Vertices& registry_entries = local_score->get_registry_entries();

  std::vector<std::vector<int32_t> > chains;
  std::vector<int32_t> chain_cov_bases;

  {
    std::vector<int32_t> new_chain;
    new_chain.reserve(lcskpp_indices.size());
    chains.push_back(new_chain);
    chain_cov_bases.push_back(0);
  }
  for (int64_t i=0; i<lcskpp_indices.size(); i++) {
    int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
    GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
    int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));
//    int32_t qpos = get128_qpos(hits[lcskpp_indices[i]]);
//    int32_t rpos = get128_rpos(hits[lcskpp_indices[i]]);

#ifdef DEBUG_TEST1
    printf ("[i = %ld] Taking the next seed. qpos = %d, rpos = %d\n", i, qpos, rpos);
    if (chains.back().size() > 0) {
      printf ("Chain info:\n");
      printf ("  get128_qpos(hits[chains.back().front() = %d, get128_qpos(hits[chains.back().back()]) = %d\n", get128_qpos(hits[chains.back().front()]), get128_qpos(hits[chains.back().back()]));
      printf ("  get128_rpos(hits[chains.back().front() = %d, get128_rpos(hits[chains.back().back()]) = %d\n", get128_rpos(hits[chains.back().front()]), get128_rpos(hits[chains.back().back()]));
    }
#endif

//    int64_t chain_len_q = (chains.back().size() == 0) ? (fwd_length) :
//                          (get128_qpos(hits[chains.back().front()]) - get128_qpos(hits[chains.back().back()]) + seed_length);
//    int64_t chain_len_r = (chains.back().size() == 0) ? (fwd_length) :
//                          (get128_rpos(hits[chains.back().front()]) - get128_rpos(hits[chains.back().back()]) + seed_length);
//    int64_t chain_len = std::max(chain_len_q, chain_len_r);
//    int64_t chain_len = std::max((int64_t) max_dist, (int64_t) chains.back().size() * lookahead_dist_factor);

//    int32_t chain_qpos_start1 = 0, chain_rpos_start1 = 0, chain_qpos_end1 = 0, chain_rpos_end1 = 0;
//    int32_t chain_qpos_start2 = 0, chain_rpos_start2 = 0, chain_qpos_end2 = 0, chain_rpos_end2 = 0;
//    if (chains.back().size() > 0) {
//      GetPositionsFromRegistry2(registry_entries, chains.back().front(), &chain_qpos_start1, &chain_rpos_start1, &chain_qpos_end1, &chain_rpos_end1);
//      GetPositionsFromRegistry2(registry_entries, chains.back().back(), &chain_qpos_start2, &chain_rpos_start2, &chain_qpos_end2, &chain_rpos_end2);
//    }
//    int64_t chain_length = std::max((chain_qpos_end1 - chain_qpos_start2), (chain_rpos_end1 - chain_rpos_end2));
//    int64_t chain_len_lookahead = std::max((int64_t) max_dist, (int64_t) chain_length * lookahead_dist_factor);
    int64_t chain_len_lookahead = max_dist;

//    printf ("chain_len_q = %ld, chain_len_r = %ld, chain_len = %ld\n", chain_len_q, chain_len_r, chain_len);

#ifdef DEBUG_TEST1
    printf ("chain_len = %ld\n", chain_len);
#endif

    int64_t max_score = 1;
    int64_t max_score_id = -1;
    int64_t next_id = (i + 1);
    while (next_id < lcskpp_indices.size()) {
//      int32_t next_qpos = get128_qpos(hits[lcskpp_indices[next_id]]);
//      int32_t next_rpos = get128_rpos(hits[lcskpp_indices[next_id]]);
      int32_t next_qpos_start = 0, next_rpos_start = 0, next_qpos_end = 0, next_rpos_end = 0;
      GetPositionsFromRegistry(registry_entries, lcskpp_indices, next_id, &next_qpos_start, &next_rpos_start, &next_qpos_end, &next_rpos_end);
      int32_t next_cov_bases = std::max((next_qpos_end - next_qpos_start), (next_rpos_end - next_rpos_start));

      #ifdef DEBUG_TEST1
            printf ("  [next_id = %ld] Considering next seed for chaining. qpos = %d, rpos = %d\n", next_id, next_qpos, next_rpos);
      #endif
      int64_t score = CalcScore(qpos_start, rpos_start, next_qpos_end, next_rpos_end, indel_bandwidth_margin, chain_len_lookahead);
      #ifdef DEBUG_TEST1
            if (score > 0) { printf ("  Break.\n"); break; }
      #endif
      if (score > 0) { break; }
      if (max_score > 0 || score > max_score) {
        #ifdef DEBUG_TEST1
                printf ("  New max score found. Previously: score = %d, score_id = %d. Now: score = %d, score_id = %d\n", max_score, max_score_id, score, next_id);
        #endif
        max_score = score;
        max_score_id = next_id;
      }
      next_id += 1;
    }
    if (max_score_id == -1) {
      chains[chains.size()-1].push_back(lcskpp_indices[i]);
      chain_cov_bases.back() += cov_bases;

      #ifdef DEBUG_TEST1
            printf ("  -> Num chains: %ld, last chain size: %ld\n", chains.size(), chains.back().size());
      #endif
      std::vector<int32_t> new_chain;
      new_chain.reserve(lcskpp_indices.size() - i);
      chains.push_back(new_chain);
      chain_cov_bases.push_back(0);
//      i = max_score_id - 1;
      #ifdef DEBUG_TEST1
            printf ("  (I) Going for next i. i = %ld\n", i);
            printf ("\n");
      #endif
      continue;
    }
    int32_t max_qpos_start = 0, max_rpos_start = 0, max_qpos_end = 0, max_rpos_end = 0;
    GetPositionsFromRegistry(registry_entries, lcskpp_indices, max_score_id, &max_qpos_start, &max_rpos_start, &max_qpos_end, &max_rpos_end);
    int32_t max_cov_bases = std::max((max_qpos_end - max_qpos_start), (max_rpos_end - max_rpos_start));
    chains[chains.size()-1].push_back(lcskpp_indices[max_score_id]);
    chain_cov_bases.back() += max_cov_bases;
    i = max_score_id - 1;
    #ifdef DEBUG_TEST1
        printf ("  -> Num chains: %ld, last chain size: %ld\n", chains.size(), chains.back().size());
        printf ("  (II) Going for next i. i = %ld\n", i);
        printf ("\n");
    #endif

//    printf ("%d %d\n", qpos, rpos);
  }

#ifdef DEBUG_TEST1
  for (int64_t i=0; i<chains.size(); i++) {
    printf ("Chain %ld:\n", i);
    for (int64_t j=0; j<chains[i].size(); j++) {
      printf ("  [%ld] %d %d\n", j, get128_qpos(hits[chains[i][j]]), get128_rpos(hits[chains[i][j]]));
    }
    printf ("\n");
  }
#endif



//  std::vector<int64_t> chain_sizes;
//  chain_sizes.reserve(chains.size());
//  for (int64_t i=0; i<chains.size(); i++) {
//    chain_sizes.push_back(chains[i].size());
//  }
//  std::sort(chain_sizes.begin(), chain_sizes.end());
//  int64_t median_size = (chain_sizes.size() > 0) ? (chain_sizes[chain_sizes.size()/4]) : (0);
//  int64_t size_cutoff = std::max((int64_t) 3, median_size);

#ifdef DEBUG_TEST1
  printf ("Generating filtered LCSk anchors. Median of chain sizes: %ld, size cutoff: %ld.\n", median_size, size_cutoff);
#endif

  ret_filtered_lcskpp_indices.size();
  if (ret_cluster_ids) {
    ret_cluster_ids->clear();
  }

  int64_t num_filtered_indices = 0;
  for (int64_t i=0; i<chains.size(); i++) {
    num_filtered_indices += chains[i].size();
  }
  ret_filtered_lcskpp_indices.reserve(num_filtered_indices);
  if (ret_cluster_ids) {
    ret_cluster_ids->reserve(num_filtered_indices);
  }
//  for (int64_t i=(chains.size()-1); i>=0; i--) {
  for (int64_t i=0; i<chains.size(); i++) {
    if (chains[i].size() <= cluster_size_cutoff && chain_cov_bases[i] < min_covered_bases) { continue; }
//    if (chain_cov_bases[i] < min_covered_bases) { continue; }
//    for (int64_t j=(chains[i].size()-1); j>=0; j--) {
    for (int64_t j=0; j<chains[i].size(); j++) {
      ret_filtered_lcskpp_indices.push_back(chains[i][j]);
      if (ret_cluster_ids) {
        ret_cluster_ids->push_back(chains.size() - i);
      }
    }
  }

#ifdef DEBUG_TEST1
  printf ("\n");
  printf ("Chain sizes:\n");
  for (int64_t i=0; i<chains.size(); i++) {
    printf ("%ld ", chains[i].size());
  }
  printf ("\n");
  for (int64_t i=0; i<chain_sizes.size(); i++) {
    printf ("%ld ", chain_sizes[i]);
  }
  printf ("\n");
#endif

#ifdef DEBUG_TEST1
  printf ("ret_filtered_lcskpp_indices.size() = %ld\n", ret_filtered_lcskpp_indices.size());
#endif

  return 0;

}
