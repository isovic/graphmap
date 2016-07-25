/*
 * lcskpp_packed_filter.cc
 *
 *  Created on: Jul 25, 2016
 *      Author: isovic
 */

#include <string>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include "owler2/lcskpp_packed.h"

//#define DEBUG_VERBOSE_LCSK_FILTER_OWLER
//#define DEBUG_VERBOSE_CLUSTERS

const int32_t chain_dist_aab = 5;
const int32_t chain_dist_dbm = 15;

#define Dist(x1, y1, x2, y2) sqrt((double) (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))

// @parameter dist_aab Accept All Bandwidths (aab) if query and reference distsaces are smaller than this value. This value should be very small (~5bp). Used to prevent breaking chains for very small indels - the indel_bandwidth_margin is not used in this case.
// @parameter dist_dbm Double Bandwidth Margin (dbm) is a heuristic which doubles the indel_bandwidth_margin if both query and reference distances are all smaller than dist_dbm. For longer homopolimer runs (which are e.g. still shorter than 20bp), there might be valid indels which have a larger gap in one dimension.
int64_t CalcScorePacked(int32_t qpos, int32_t rpos, int32_t next_qpos, int32_t next_rpos, double indel_bandwidth_margin, int32_t fwd_length, int32_t dist_aab, int32_t dist_dbm, double *score_gap, double *score_dist) {
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

  int32_t x_dist = abs(qpos - next_qpos);
  int32_t y_dist = abs(rpos - next_rpos);

  int32_t bandwidth = indel_bandwidth_margin * dist;
  if (x_dist < dist_aab && y_dist < dist_aab) { bandwidth = dist; }
  else if (x_dist < dist_dbm && y_dist < dist_dbm) { bandwidth = indel_bandwidth_margin * 2.0f * dist; }

//  int32_t bandwidth = indel_bandwidth_margin * fwd_length;
//  printf ("  -> fwd_length = %d, dist = %f, mismatch = %f, band = %f\n", fwd_length, dist, mismatch, band);
//  fflush(stdout);

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
    // All numbers are too high, break the chain.
    return 3;
  }

  return 4;
}

int FilterAnchorsByChainingPacked(const uint128_t *packed_hits, int64_t n_hits, int32_t k, const SingleSequence* read, const ProgramParameters *parameters,
                  const std::vector<int> &lcskpp_indices, double indel_bandwidth_margin, int32_t max_dist, int32_t lookahead_dist_factor, int64_t min_covered_bases, int32_t cluster_size_cutoff,
                  std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {

  int64_t chain_len_lookahead = max_dist;

//  const Vertices& registry_entries = local_score->get_registry_entries();

  std::vector<std::vector<int32_t> > chains;
  std::vector<int32_t> chain_cov_bases;

  {
    std::vector<int32_t> new_chain;
    new_chain.reserve(lcskpp_indices.size());
    chains.push_back(new_chain);
    chain_cov_bases.push_back(0);
  }

#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
  LOG_DEBUG_SPEC_NO_HEADER("Starting a new chain: %ld\n", chains.size());
#endif

//  fflush(stdout);
  for (int64_t i=0; i<lcskpp_indices.size(); i++) {
    double max_score_gap = 1.0, max_score_mismatch = 1.0;
    int64_t max_score_id = -1;
    int64_t next_id = i;
    while (next_id < lcskpp_indices.size()) {
      if (chains.back().size() == 0) {
//        int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
//        GetPositionsFromRegistry(registry_entries, lcskpp_indices, next_id, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
        int32_t qpos_start = get_lcsk128_qpos(packed_hits[lcskpp_indices[next_id]]);
        int32_t rpos_start = get_lcsk128_rpos(packed_hits[lcskpp_indices[next_id]]);
        int32_t qpos_end = qpos_start + k, rpos_end = rpos_start + k;

        int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));
        chains.back().push_back(lcskpp_indices[next_id]);
        chain_cov_bases.back() += cov_bases;
#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
        LOG_DEBUG_SPEC_NO_HEADER("  (I) Adding ancho %ld (chains.back().size() == 0): start = [%d, %d], end = [%d, %d]\n", 0, qpos_start, rpos_start, qpos_end, rpos_end);
#endif
//        fflush(stdout);
//        next_id += 1;
        continue;
      }

//      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
//      GetPositionsFromRegistry2(registry_entries, chains.back().back(), &qpos_start, &rpos_start, &qpos_end, &rpos_end);
      int32_t qpos_start = get_lcsk128_qpos(packed_hits[chains.back().back()]);
      int32_t rpos_start = get_lcsk128_rpos(packed_hits[chains.back().back()]);
      int32_t qpos_end = qpos_start + k, rpos_end = rpos_start + k;
      int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));

//      int32_t next_qpos_start = 0, next_rpos_start = 0, next_qpos_end = 0, next_rpos_end = 0;
//      GetPositionsFromRegistry(registry_entries, lcskpp_indices, next_id, &next_qpos_start, &next_rpos_start, &next_qpos_end, &next_rpos_end);
      int32_t next_qpos_start = get_lcsk128_qpos(packed_hits[lcskpp_indices[next_id]]);
      int32_t next_rpos_start = get_lcsk128_rpos(packed_hits[lcskpp_indices[next_id]]);
      int32_t next_qpos_end = next_qpos_start + k, next_rpos_end = next_rpos_start + k;
      int32_t next_cov_bases = std::max((next_qpos_end - next_qpos_start), (next_rpos_end - next_rpos_start));

#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
      LOG_DEBUG_SPEC_NO_HEADER("  Considering anchor %ld (i = %ld, lcskpp_indices[%ld] = %d, chains.back().back() = %d): start = [%d, %d], end = [%d, %d]. max_score_id = %ld\n", next_id, i, next_id, lcskpp_indices[next_id], chains.back().back(), next_qpos_start, next_rpos_start, next_qpos_end, next_rpos_end, max_score_id);
#endif

      double score_gap = 0.0, score_mismatch = 0.0;
      int64_t ret_val = CalcScorePacked(qpos_start, rpos_start, next_qpos_end, next_rpos_end, indel_bandwidth_margin, chain_len_lookahead, chain_dist_aab, chain_dist_dbm, &score_gap, &score_mismatch);

      if (ret_val > 0) {
#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
        LOG_DEBUG_SPEC_NO_HEADER("    - Breaking. CalcScore returned a value > 0 (ret_val = %ld)\n", ret_val);
#endif
        break;
      }

      if (max_score_id == -1 || (score_gap >= max_score_gap && score_mismatch >= max_score_mismatch)) {
#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
        LOG_DEBUG_SPEC_NO_HEADER("    - This will be the new max_score_id.\n", ret_val);
#endif
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

#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
      LOG_DEBUG_SPEC_NO_HEADER("  End of chain.\n\n");
      LOG_DEBUG_SPEC_NO_HEADER("  Creating a new chain (new chains.size() = %ld\n", chains.size());
#endif

//      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
//      GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
      int32_t qpos_start = get_lcsk128_qpos(packed_hits[lcskpp_indices[i]]);
      int32_t rpos_start = get_lcsk128_rpos(packed_hits[lcskpp_indices[i]]);
      int32_t qpos_end = qpos_start + k, rpos_end = rpos_start + k;
      int32_t cov_bases = std::max((qpos_end - qpos_start), (rpos_end - rpos_start));
      chains.back().push_back(lcskpp_indices[i]);
      chain_cov_bases.back() += cov_bases;

#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
      LOG_DEBUG_SPEC_NO_HEADER("  (II) Adding anchor %ld (max_score_id == -1): start = [%d, %d], end = [%d, %d]. (Case where max_score_id == -1)\n", i, qpos_start, rpos_start, qpos_end, rpos_end);
#endif
//      fflush(stdout);
//      LOG_DEBUG_SPEC_NO_HEADER("Starting a new chain: %ld\n", chains.size());
//      fflush(stdout);

    } else {
//    } else if (chains.back().back() != lcskpp_indices[max_score_id]) {  // In this case, there was a max_score found, but we also check if the lcskpp index was already added to the array.
//      int32_t max_qpos_start = 0, max_rpos_start = 0, max_qpos_end = 0, max_rpos_end = 0;
//      GetPositionsFromRegistry(registry_entries, lcskpp_indices, max_score_id, &max_qpos_start, &max_rpos_start, &max_qpos_end, &max_rpos_end);
      int32_t max_qpos_start = get_lcsk128_qpos(packed_hits[lcskpp_indices[max_score_id]]);
      int32_t max_rpos_start = get_lcsk128_rpos(packed_hits[lcskpp_indices[max_score_id]]);
      int32_t max_qpos_end = max_qpos_start + k, max_rpos_end = max_rpos_start + k;
      int32_t max_cov_bases = std::max((max_qpos_end - max_qpos_start), (max_rpos_end - max_rpos_start));

      chains.back().push_back(lcskpp_indices[max_score_id]);
      chain_cov_bases.back() += max_cov_bases;
      i = max_score_id;

#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
      LOG_DEBUG_SPEC_NO_HEADER("  (III) Adding anchor (max_score_id = %ld), start = [%d, %d], end = [%d, %d]\n", max_score_id, max_qpos_start, max_rpos_start, max_qpos_end, max_rpos_end);
#endif
//      fflush(stdout);
    }
  }

  ret_filtered_lcskpp_indices.size();
  if (ret_cluster_ids) { ret_cluster_ids->clear(); }



#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
  LOG_DEBUG_SPEC_NO_HEADER("Filtering short chains.\n");
#endif

  for (int64_t i=0; i<chains.size(); i++) {
    // Remove chains with too low of a count.
    if (chains[i].size() == 0) { continue; }

//    int32_t qpos_start1 = 0, rpos_start1 = 0, qpos_end1 = 0, rpos_end1 = 0;
//    GetPositionsFromRegistry2(registry_entries, chains[i].front(), &qpos_start1, &rpos_start1, &qpos_end1, &rpos_end1);
    int32_t qpos_start1 = get_lcsk128_qpos(packed_hits[chains[i].front()]);
    int32_t rpos_start1 = get_lcsk128_rpos(packed_hits[chains[i].front()]);
    int32_t qpos_end1 = qpos_start1 + k, rpos_end1 = rpos_start1 + k;
//    int32_t qpos_start2 = 0, rpos_start2 = 0, qpos_end2 = 0, rpos_end2 = 0;
//    GetPositionsFromRegistry2(registry_entries, chains[i].back(), &qpos_start2, &rpos_start2, &qpos_end2, &rpos_end2);
    int32_t qpos_start2 = get_lcsk128_qpos(packed_hits[chains[i].back()]);
    int32_t rpos_start2 = get_lcsk128_rpos(packed_hits[chains[i].back()]);
    int32_t qpos_end2 = qpos_start2 + k, rpos_end2 = rpos_start2 + k;

//    printf ("chains[%ld].size() = %ld, chain_cov_bases[%ld] = %ld, qlen = %d, rlen = %d, qpos_start1 = %d, qpos_end2 = %d;   ", i, chains[i].size(), i, chain_cov_bases[i], (qpos_end2 - qpos_start1), (rpos_end2 - rpos_start1), qpos_start1, qpos_end2);
    if (chains[i].size() <= cluster_size_cutoff && chain_cov_bases[i] < min_covered_bases) {
#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
      LOG_DEBUG_SPEC_NO_HEADER("Removed!");
#endif
      chains[i].clear();
    }
#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
    LOG_DEBUG_SPEC_NO_HEADER("\n");
#endif
  }

#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
  LOG_DEBUG_SPEC_NO_HEADER("cluster_size_cutoff = %ld\n", cluster_size_cutoff);
  LOG_DEBUG_SPEC_NO_HEADER("min_covered_bases = %ld\n", min_covered_bases);
#endif

  int64_t num_filtered_indices = 0;
  for (int64_t i=0; i<chains.size(); i++) { num_filtered_indices += chains[i].size(); }
  ret_filtered_lcskpp_indices.reserve(num_filtered_indices);
  if (ret_cluster_ids) { ret_cluster_ids->reserve(num_filtered_indices); }

  // Fill the return indices.
  for (int64_t i=0; i<chains.size(); i++) {
    // Remove chains with too low of a count.
    if (chains[i].size() <= cluster_size_cutoff && chain_cov_bases[i] < min_covered_bases) { continue; }

//    int32_t qpos_start1 = 0, rpos_start1 = 0, qpos_end1 = 0, rpos_end1 = 0;
//    GetPositionsFromRegistry2(registry_entries, chains[i].front(), &qpos_start1, &rpos_start1, &qpos_end1, &rpos_end1);
    int32_t qpos_start1 = get_lcsk128_qpos(packed_hits[chains[i].front()]);
    int32_t rpos_start1 = get_lcsk128_rpos(packed_hits[chains[i].front()]);
    int32_t qpos_end1 = qpos_start1 + k, rpos_end1 = rpos_start1 + k;
//    int32_t qpos_start2 = 0, rpos_start2 = 0, qpos_end2 = 0, rpos_end2 = 0;
//    GetPositionsFromRegistry2(registry_entries, chains[i].back(), &qpos_start2, &rpos_start2, &qpos_end2, &rpos_end2);
    int32_t qpos_start2 = get_lcsk128_qpos(packed_hits[chains[i].back()]);
    int32_t rpos_start2 = get_lcsk128_rpos(packed_hits[chains[i].back()]);
    int32_t qpos_end2 = qpos_start2 + k, rpos_end2 = rpos_start2 + k;

#ifdef DEBUG_VERBOSE_LCSK_FILTER_OWLER
    LOG_DEBUG_SPEC_NO_HEADER("chains[%ld].size() = %ld, chain_cov_bases[%ld] = %ld, qlen = %d, rlen = %d, qpos_start1 = %d, qpos_end2 = %d;   ", i, chains[i].size(), i, chain_cov_bases[i], (qpos_end2 - qpos_start1), (rpos_end2 - rpos_start1), qpos_start1, qpos_end2);
    LOG_DEBUG_SPEC_NO_HEADER("\n");
#endif

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

int GenerateClustersPacked(const uint128_t *packed_hits, int64_t n_hits, int32_t k, int64_t min_num_anchors_in_cluster, int64_t min_cluster_length, int64_t min_cluster_covered_bases,
                     float min_cluster_coverage, std::vector<int> &lcskpp_indices,
                     const SingleSequence* read, const ProgramParameters* parameters, std::vector<ClusterAndIndices *> &ret_clusters,
                     std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {

  ret_clusters.clear();

#ifdef DEBUG_VERBOSE_CLUSTERS
  LOG_DEBUG_SPEC("Starting to generate clusters from lcskpp filtered anchors. lcskpp_indices.size() = %ld\n", lcskpp_indices.size());
#endif

  if (lcskpp_indices.size() == 0) { return 1; }

  ClusterAndIndices *new_cluster = NULL;

  new_cluster = new ClusterAndIndices;

  double indel_bandwidth_margin = 0.20f;
  double score_gap = 0.0, score_dist = 0.0;

  for (int64_t i=(lcskpp_indices.size() - 1); i >= 0; i--) {
    if (new_cluster->lcskpp_indices.size() == 0) {
      new_cluster->lcskpp_indices.push_back(lcskpp_indices[i]);
      new_cluster->coverage += k; // local_score->get_registry_entries().covered_bases_queries[lcskpp_indices[i]];
      new_cluster->num_anchors += 1;

#ifdef DEBUG_VERBOSE_CLUSTERS
      LOG_DEBUG_SPEC("  Adding a new lcskpp index to an empty cluster. Cluster ID: %ld.\n", ret_clusters.size());
#endif

    } else {
      int32_t qpos_start = get_lcsk128_qpos(packed_hits[lcskpp_indices[i]]);
      int32_t rpos_start = get_lcsk128_rpos(packed_hits[lcskpp_indices[i]]);
      int32_t qpos_end = qpos_start + k, rpos_end = rpos_start + k;

      int32_t prev_qpos_start = get_lcsk128_qpos(packed_hits[lcskpp_indices[(i + 1)]]);
      int32_t prev_rpos_start = get_lcsk128_rpos(packed_hits[lcskpp_indices[(i + 1)]]);
      int32_t prev_qpos_end = prev_qpos_start + k, prev_rpos_end = prev_rpos_start + k;

      int64_t ret_val = CalcScorePacked(qpos_start, rpos_start, prev_qpos_end, prev_rpos_end, indel_bandwidth_margin, -1, chain_dist_aab, chain_dist_dbm, &score_gap, &score_dist);

      if (ret_val > 0) {
#ifdef DEBUG_VERBOSE_CLUSTERS
        LOG_DEBUG_SPEC("  Closing the cluster. Cluster ID: %ld, ret_val = %ld.\n\n", ret_clusters.size(), ret_val);
#endif
        ret_clusters.push_back(new_cluster);
        new_cluster = new ClusterAndIndices;
      }

      new_cluster->lcskpp_indices.push_back(lcskpp_indices[i]);
      new_cluster->coverage += k; // local_score->get_registry_entries().covered_bases_queries[lcskpp_indices[i]];
      new_cluster->num_anchors += 1;
#ifdef DEBUG_VERBOSE_CLUSTERS
      LOG_DEBUG_SPEC("  Adding a new lcskpp index to Cluster ID: %ld, new_cluster->lcskpp_indices.size() = %ld (after adding).\n", ret_clusters.size(), new_cluster->lcskpp_indices.size());
#endif
    }
  }

  if (new_cluster->lcskpp_indices.size() != 0) {
#ifdef DEBUG_VERBOSE_CLUSTERS
    LOG_DEBUG_SPEC("  Closing the cluster. Cluster ID: %ld.\n\n", ret_clusters.size());
#endif
    ret_clusters.push_back(new_cluster);
  }

#ifdef DEBUG_VERBOSE_CLUSTERS
  LOG_DEBUG_SPEC_NEWLINE;
#endif

  /// Generate the final ret_clusters that will be returned and used further.
  for (int64_t i=0; i<ret_clusters.size(); i++) {
    ret_clusters[i]->query.start = get_lcsk128_qpos(packed_hits[ret_clusters[i]->lcskpp_indices.front()]); // registry_entries.query_starts[ret_clusters[i]->lcskpp_indices.front()];
    ret_clusters[i]->query.end = get_lcsk128_qpos(packed_hits[ret_clusters[i]->lcskpp_indices.back()]) + k - 1; // registry_entries.query_ends[ret_clusters[i]->lcskpp_indices.back()] - 1;
    ret_clusters[i]->ref.start = get_lcsk128_rpos(packed_hits[ret_clusters[i]->lcskpp_indices.front()]); // registry_entries.reference_starts[ret_clusters[i]->lcskpp_indices.front()];
    ret_clusters[i]->ref.end = get_lcsk128_rpos(packed_hits[ret_clusters[i]->lcskpp_indices.back()]) + k - 1; // registry_entries.reference_ends[ret_clusters[i]->lcskpp_indices.back()] - 1;

#ifdef DEBUG_VERBOSE_CLUSTERS
    LOG_DEBUG_SPEC("Anchors in Cluster %ld:\n", i);
    for (int64_t j=0; j<ret_clusters[i]->lcskpp_indices.size(); j++) {
      LOG_DEBUG_SPEC("  start = [%ld, %ld], end = [%ld, %ld]\n",
                     get_lcsk128_qpos(packed_hits[ret_clusters[i]->lcskpp_indices[j]]), get_lcsk128_rpos(packed_hits[ret_clusters[i]->lcskpp_indices[j]]),
                     get_lcsk128_qpos(packed_hits[ret_clusters[i]->lcskpp_indices[j]]) + k - 1, get_lcsk128_rpos(packed_hits[ret_clusters[i]->lcskpp_indices[j]]) + k - 1);
    }
    LOG_DEBUG_SPEC_NEWLINE;
#endif

    int64_t cluster_length = ret_clusters[i]->query.end - ret_clusters[i]->query.start + 1;
    int64_t covered_bases = ret_clusters[i]->coverage;
    float cluster_coverage = ((float) covered_bases) / ((float) cluster_length);

    // Filter clusters which are too small.
    if ((ret_clusters[i]->lcskpp_indices.size() <= min_num_anchors_in_cluster &&
        (cluster_length < min_cluster_length)) || covered_bases < min_cluster_covered_bases) {

#ifdef DEBUG_VERBOSE_CLUSTERS
      LOG_DEBUG_SPEC("  Filtering cluster %ld because lengths not satisfied. ret_clusters[i]->lcskpp_indices.size() = %ld / %ld, cluster_length = %ld / %ld, covered_bases = %ld / %ld\n",
                     i, ret_clusters[i]->lcskpp_indices.size(), min_num_anchors_in_cluster, cluster_length, min_cluster_length, covered_bases, min_cluster_covered_bases);
#endif
      ret_clusters[i]->lcskpp_indices.clear();
      continue;
    }

    ret_filtered_lcskpp_indices.insert(ret_filtered_lcskpp_indices.end(), ret_clusters[i]->lcskpp_indices.begin(), ret_clusters[i]->lcskpp_indices.end());

    /// Create indices for debugging purposes (so we can differentiate clusters).
    if (ret_cluster_ids) {
      std::vector<int32_t> cluster_indices(ret_clusters[i]->lcskpp_indices.size(), i);
      ret_cluster_ids->insert(ret_cluster_ids->end(), cluster_indices.begin(), cluster_indices.end());
    }

#ifdef DEBUG_VERBOSE_CLUSTERS
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[Cluster %ld] cluster_length = %ld, covered_bases = %ld\n", i, cluster_length, covered_bases), std::string(__FUNCTION__));
    for (int64_t j=0; j<ret_clusters[i]->lcskpp_indices.size(); j++) {
      int32_t qpos_start = get_lcsk128_qpos(packed_hits[ret_clusters[i]->lcskpp_indices[j]]);
      int32_t rpos_start = get_lcsk128_rpos(packed_hits[ret_clusters[i]->lcskpp_indices[j]]);
      int32_t qpos_end = qpos_start + k, rpos_end = rpos_start + k;
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("    [%ld] start: [%d, %d], end: [%d, %d]\n", j, qpos_start, rpos_start, qpos_end, rpos_end), "[]");
    }
#endif
  }

  for (int64_t i = (ret_clusters.size() - 1); i >= 0; i--) {
    if (ret_clusters[i]->lcskpp_indices.size() == 0) {
      ret_clusters.erase(ret_clusters.begin() + i);
    }
  }

#ifdef DEBUG_VERBOSE_CLUSTERS
  LOG_DEBUG_SPEC_NEWLINE;
#endif

  return 0;
}
