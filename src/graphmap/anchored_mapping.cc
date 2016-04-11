/*
 * experimental.cc
 *
 *  Created on: May 1, 2015
 *      Author: isovic
 */

#include "graphmap/graphmap.h"
#include "graphmap/filter_anchors.h"

/// These are some constants used for filtering shady anchors.
/// TODO: This can be omitted if dynamic programming was used to penalize the anchor distances.
/// int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * MIN_CLUSTER_COVERAGE_FACTOR;
#define MIN_CLUSTER_COVERAGE_FACTOR 0.05f
/// int64_t min_cluster_length = read->get_sequence_length() * MIN_CLUSTER_LENGTH_FACTOR;
#define MIN_CLUSTER_LENGTH_FACTOR 0.03f

int GenerateClusters(int64_t min_num_anchors_in_cluster, int64_t min_cluster_length, int64_t min_cluster_covered_bases, float min_cluster_coverage, std::vector<int> &lcskpp_indices,
                     ScoreRegistry* local_score, MappingData* mapping_data, const std::vector<Index *> indexes,
                     const SingleSequence* read, const ProgramParameters* parameters, std::vector<ClusterAndIndices *> &ret_clusters,
                     std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {
//  std::vector<ClusterAndIndices *> clusters;
  ret_clusters.clear();

  if (lcskpp_indices.size() == 0) { return 1; }

  ClusterAndIndices *new_cluster = NULL;

  new_cluster = new ClusterAndIndices;

//  printf ("lcskpp_indices.size() = %ld\n", lcskpp_indices.size());
//  printf ("local_score->get_registry_entries().size() = %ld\n", local_score->get_registry_entries().num_vertices);
//  printf ("lcskpp_indices.back() = %ld\n", lcskpp_indices.back());
//  printf ("lcskpp_indices.front() = %ld\n", lcskpp_indices.front());
//  fflush(stdout);

//  new_cluster->query.start = local_score->get_registry_entries().query_starts[lcskpp_indices.back()];
//  new_cluster->query.end = local_score->get_registry_entries().query_ends[lcskpp_indices.front()] - 1;      // The '- 1' is because anchor end coordinate points to one base after the last inclusive base. Clusters should demarcate the exact start and end positions (start is the first base, end is the last base).
//  new_cluster->ref.start = local_score->get_registry_entries().reference_starts[lcskpp_indices.back()];
//  new_cluster->ref.end = local_score->get_registry_entries().reference_ends[lcskpp_indices.front()] - 1;
//

  const Vertices& registry_entries = local_score->get_registry_entries();

  double indel_bandwidth_margin = 0.20f;
  double score_gap = 0.0, score_dist = 0.0;

  for (int64_t i=(lcskpp_indices.size() - 1); i >= 0; i--) {
    if (new_cluster->lcskpp_indices.size() == 0) {
      new_cluster->lcskpp_indices.push_back(lcskpp_indices[i]);
      new_cluster->coverage += local_score->get_registry_entries().covered_bases_queries[lcskpp_indices[i]];
      new_cluster->num_anchors += 1;

//      printf ("Tu sam 1! i = %ld / %ld, new_cluster->lcskpp_indices.size() = %ld\n", i, (lcskpp_indices.size() - 1), new_cluster->lcskpp_indices.size());
//      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
//      GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
//      printf ("    start: [%d, %d], end: [%d, %d], i = %ld\n", qpos_start, rpos_start, qpos_end, rpos_end, i);
//      fflush(stdout);

    } else {
      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
      GetPositionsFromRegistry(registry_entries, lcskpp_indices, i, &qpos_start, &rpos_start, &qpos_end, &rpos_end);
//      printf ("    start: [%d, %d], end: [%d, %d], i = %ld\n", qpos_start, rpos_start, qpos_end, rpos_end, i);
//      fflush(stdout);

      int32_t prev_qpos_start = 0, prev_rpos_start = 0, prev_qpos_end = 0, prev_rpos_end = 0;
      GetPositionsFromRegistry(registry_entries, lcskpp_indices, (i + 1), &prev_qpos_start, &prev_rpos_start, &prev_qpos_end, &prev_rpos_end);

      int64_t ret_val = CalcScore(qpos_start, rpos_start, prev_qpos_end, prev_rpos_end, indel_bandwidth_margin, -1, &score_gap, &score_dist);

      if (ret_val > 0) {
        ret_clusters.push_back(new_cluster);
        new_cluster = new ClusterAndIndices;
      }
      new_cluster->lcskpp_indices.push_back(lcskpp_indices[i]);
      new_cluster->coverage += local_score->get_registry_entries().covered_bases_queries[lcskpp_indices[i]];
      new_cluster->num_anchors += 1;
//      printf ("Tu sam 2! i = %ld / %ld, new_cluster->lcskpp_indices.size() = %ld\n", i, (lcskpp_indices.size() - 1), new_cluster->lcskpp_indices.size());
//      fflush(stdout);
    }
  }

  if (new_cluster->lcskpp_indices.size() != 0) {
    ret_clusters.push_back(new_cluster);
  }

  /// Generate the final ret_clusters that will be returned and used further.
  for (int64_t i=0; i<ret_clusters.size(); i++) {
    ret_clusters[i]->query.start = registry_entries.query_starts[ret_clusters[i]->lcskpp_indices.front()];
    ret_clusters[i]->query.end = registry_entries.query_ends[ret_clusters[i]->lcskpp_indices.back()] - 1;
    ret_clusters[i]->ref.start = registry_entries.reference_starts[ret_clusters[i]->lcskpp_indices.front()];
    ret_clusters[i]->ref.end = registry_entries.reference_ends[ret_clusters[i]->lcskpp_indices.back()] - 1;

    int64_t cluster_length = ret_clusters[i]->query.end - ret_clusters[i]->query.start + 1;
    int64_t covered_bases = ret_clusters[i]->coverage;
    float cluster_coverage = ((float) covered_bases) / ((float) cluster_length);

    // Filter clusters which are too small.
//    if ((min_num_anchors_in_cluster > 0 && min_cluster_length > 0 && min_cluster_coverage > 0) && ret_clusters[i]->lcskpp_indices.size() <= min_num_anchors_in_cluster && cluster_length < min_cluster_length) {
    if (ret_clusters[i]->lcskpp_indices.size() <= min_num_anchors_in_cluster &&
        cluster_length < min_cluster_length || covered_bases < min_cluster_covered_bases) {
      ret_clusters[i]->lcskpp_indices.clear();
      continue;
    }

    ret_filtered_lcskpp_indices.insert(ret_filtered_lcskpp_indices.end(), ret_clusters[i]->lcskpp_indices.begin(), ret_clusters[i]->lcskpp_indices.end());
//    printf ("ret_filtered_lcskpp_indices.size() = %ld\n", ret_filtered_lcskpp_indices.size());
//    fflush(stdout);

    /// Create indices for debugging purposes (so we can differentiate clusters).
    if (ret_cluster_ids) {
      std::vector<int32_t> cluster_indices(ret_clusters[i]->lcskpp_indices.size(), i);
      ret_cluster_ids->insert(ret_cluster_ids->end(), cluster_indices.begin(), cluster_indices.end());
    }

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[Cluster %ld] cluster_length = %ld, covered_bases = %ld\n", i, cluster_length, covered_bases), std::string(__FUNCTION__));
    for (int64_t j=0; j<ret_clusters[i]->lcskpp_indices.size(); j++) {
      int32_t qpos_start = 0, rpos_start = 0, qpos_end = 0, rpos_end = 0;
      GetPositionsFromRegistry2(registry_entries, ret_clusters[i]->lcskpp_indices[j], &qpos_start, &rpos_start, &qpos_end, &rpos_end);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("    [%ld] start: [%d, %d], end: [%d, %d]\n", j, qpos_start, rpos_start, qpos_end, rpos_end), "[]");
    }
  }

  for (int64_t i = (ret_clusters.size() - 1); i >= 0; i--) {
    if (ret_clusters[i]->lcskpp_indices.size() == 0) {
      ret_clusters.erase(ret_clusters.begin() + i);
    }
  }

  return 0;
}

int GenerateClustersDummy(int64_t min_cluster_length, float min_cluster_coverage, std::vector<int> &lcskpp_indices,
                     ScoreRegistry* local_score, MappingData* mapping_data, const std::vector<Index *> indexes,
                     const SingleSequence* read, const ProgramParameters* parameters, std::vector<ClusterAndIndices *> &ret_clusters,
                     std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {
//  std::vector<ClusterAndIndices *> clusters;
  ret_clusters.clear();

  if (lcskpp_indices.size() == 0) { return 1; }

  ClusterAndIndices *new_cluster = NULL;

  new_cluster = new ClusterAndIndices;

//  printf ("lcskpp_indices.size() = %ld\n", lcskpp_indices.size());
//  printf ("local_score->get_registry_entries().size() = %ld\n", local_score->get_registry_entries().num_vertices);
//  printf ("lcskpp_indices.back() = %ld\n", lcskpp_indices.back());
//  printf ("lcskpp_indices.front() = %ld\n", lcskpp_indices.front());
//  fflush(stdout);

  new_cluster->query.start = local_score->get_registry_entries().query_starts[lcskpp_indices.back()];
  new_cluster->query.end = local_score->get_registry_entries().query_ends[lcskpp_indices.front()] - 1;      // The '- 1' is because anchor end coordinate points to one base after the last inclusive base. Clusters should demarcate the exact start and end positions (start is the first base, end is the last base).
  new_cluster->ref.start = local_score->get_registry_entries().reference_starts[lcskpp_indices.back()];
  new_cluster->ref.end = local_score->get_registry_entries().reference_ends[lcskpp_indices.front()] - 1;

  for (int64_t i=(lcskpp_indices.size() - 1); i >= 0; i--) {
    new_cluster->lcskpp_indices.push_back(lcskpp_indices[i]);
    new_cluster->coverage += local_score->get_registry_entries().covered_bases_queries[lcskpp_indices[i]];
    new_cluster->num_anchors += 1;
  }

  ret_clusters.push_back(new_cluster);

  /// Generate the final ret_clusters that will be returned and used further.
  for (int64_t i=0; i<ret_clusters.size(); i++) {
    int64_t cluster_length = ret_clusters[i]->query.end - ret_clusters[i]->query.start + 1;
    int64_t covered_bases = ret_clusters[i]->coverage;
    float cluster_coverage = ((float) covered_bases) / ((float) cluster_length);

    ret_filtered_lcskpp_indices.insert(ret_filtered_lcskpp_indices.end(), ret_clusters[i]->lcskpp_indices.begin(), ret_clusters[i]->lcskpp_indices.end());

    /// Create indices for debugging purposes (so we can differentiate clusters).
    if (ret_cluster_ids) {
      std::vector<int32_t> cluster_indices(ret_clusters[i]->lcskpp_indices.size(), i);
      ret_cluster_ids->insert(ret_cluster_ids->end(), cluster_indices.begin(), cluster_indices.end());
    }

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[Cluster %ld] cluster_length = %ld, covered_bases = %ld\n", i, cluster_length, covered_bases), "L1-PostProcessRegionWithLCS_");
  }

  return 0;
}

int GraphMap::VerboseClustersToFile_(std::string out_file, const ScoreRegistry* local_score, const MappingData* mapping_data, const std::vector<Index *> &indexes, const SingleSequence* read, const ProgramParameters* parameters, const std::vector<ClusterAndIndices *> &clusters) {
  FILE *fp = fopen(out_file.c_str(), "a");
  if (fp == NULL) { return 1; }

  int64_t num_nonempty_clusters = 0;
  for (int64_t i=0; i<clusters.size(); i++) {
    if (clusters[i] && clusters[i]->lcskpp_indices.size() > 0) { num_nonempty_clusters += 1; }
  }

  // 1. Number of clusters, 2. Read ID, 3. Read len, 4. Target ID, 4. Target len, 5. Target reversed
  fprintf (fp, "region\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%d\n", local_score->get_scores_id(), num_nonempty_clusters, read->get_sequence_id(), read->get_sequence_length(),
           (local_score->get_region().reference_id % indexes[0]->get_num_sequences_forward()), indexes[0]->get_reference_lengths()[local_score->get_region().reference_id], (int32_t) local_score->get_region().reference_id >= indexes[0]->get_num_sequences_forward());
  int64_t current_cluster = 0;
  for (int64_t i=0; i<clusters.size(); i++) {
    if (clusters[i] && clusters[i]->lcskpp_indices.size() > 0) {
      current_cluster += 1;
      // Cluster line:
      //  cluster_id qstart qend rstart rend num_anchors num_covered_bases
      fprintf (fp, "cluster\t%ld\t%ld\t%ld\t%ld\t%ld\t%d\t%d\n", current_cluster, clusters[i]->query.start, clusters[i]->query.end, clusters[i]->ref.start, clusters[i]->ref.end, clusters[i]->num_anchors, clusters[i]->coverage);
    }
  }
  fprintf (fp, "#\n");
  fclose(fp);
  return 0;
}

int GraphMap::AnchoredPostProcessRegionWithLCS_(ScoreRegistry* local_score, MappingData* mapping_data, const std::vector<Index *> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || ((int64_t) read->get_sequence_id()) == parameters->debug_read), FormatString("Entering function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB] current_readid = %ld, current_local_score = %ld\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024), read->get_sequence_id(), local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS");
  int lcskpp_length = 0;
  std::vector<int> lcskpp_indices;
  CalcLCSFromLocalScoresCacheFriendly_(&(local_score->get_registry_entries()), false, 0, 0, &lcskpp_length, &lcskpp_indices);
  if (lcskpp_length == 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Current local scores: %ld, lcskpp_length == 0 || best_score == NULL\n", local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS");
    return 1;
  }

  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LOG_DEBUG_SPEC("After LCSk:\n", local_score->get_scores_id());
    for (int64_t i = 0; i < lcskpp_indices.size(); i++) {
      LOG_DEBUG_SPEC_NO_HEADER("[%ld] %s\n", i, local_score->get_registry_entries().VerboseToString(lcskpp_indices[i]).c_str());
    }
    LOG_DEBUG_SPEC_NEWLINE;
  }

//  std::reverse(first_filtered_lcskpp_indices.begin(), first_filtered_lcskpp_indices.end());

  /// Filter the LCSk anchors.
  std::vector<ClusterAndIndices *> clusters;
  std::vector<int> cluster_indices;
  std::vector<int32_t> cluster_ids;

  /// Filter the LCSk anchors, first pass. This pass filters outliers, but does not generate clusters.
  std::vector<int> first_filtered_lcskpp_indices;
  std::vector<int> second_filtered_lcskpp_indices;

  FilterAnchorsByChaining(read, local_score, parameters, lcskpp_indices, parameters->error_rate/2 + 0.01f, 200.0f, 0, 50, 2, first_filtered_lcskpp_indices, NULL);
  FilterAnchorsByDiff(read, local_score, parameters, first_filtered_lcskpp_indices, second_filtered_lcskpp_indices);
  GenerateClusters(2, 50, 50, 0.0, second_filtered_lcskpp_indices, local_score, mapping_data, indexes, read, parameters, clusters, cluster_indices, &cluster_ids);

  // Find the L1 parameters (median line and the confidence intervals).
  float l_diff = read->get_sequence_length() * parameters->error_rate;
  float maximum_allowed_deviation = l_diff * sqrt(2.0f) / 2.0f;
  float sigma_L2 = 0.0f, confidence_L1 = 0.0f;
  int64_t k = 0, l = 0;
  // Actuall L1 calculation.
  int ret_L1 = CalculateL1ParametersWithMaximumDeviation_(local_score, cluster_indices, maximum_allowed_deviation, &k, &l, &sigma_L2, &confidence_L1);
  // Sanity check.
  if (ret_L1) {
    LOG_DEBUG_SPEC("An error occured, L1 function (I) returned with %d!\n", ret_L1);
    if (ret_L1 == 1) {  LOG_DEBUG_SPEC("  lcskpp_indices.size() == 0\n"); }
    else if (ret_L1 == 2) { LOG_DEBUG_SPEC("  num_points_under_max_dev_threshold == 0\n"); }

    #ifndef RELEASE_VERSION
      if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
        LOG_DEBUG_SPEC("Writing all anchors to file scores-%ld.\n",  local_score->get_scores_id());
        VerboseLocalScoresToFile(FormatString("temp/local_scores/scores-%ld.csv", local_score->get_scores_id()), read, local_score, NULL, 0, 0, false);
        LOG_DEBUG_SPEC("Writing LCSk anchors to file LCS-%ld.\n",  local_score->get_scores_id());
        VerboseLocalScoresToFile(FormatString("temp/local_scores/LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, 0, 0, false);
        LOG_DEBUG_SPEC("Writing LCSk anchors to file LCSL1-%ld.\n",  local_score->get_scores_id());
        VerboseLocalScoresToFile(FormatString("temp/local_scores/LCSL1-%ld.csv", local_score->get_scores_id()), read, local_score, &second_filtered_lcskpp_indices, 0, 0, false);
        LOG_DEBUG_SPEC("Writing LCSk anchors to file double_LCS-%ld.\n",  local_score->get_scores_id());
        VerboseLocalScoresToFile(FormatString("temp/local_scores/double_LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &cluster_indices, 0, 0, false, &cluster_ids);
      }
    #endif

    return 1;
  }
  float allowed_L1_deviation = 3.0f * confidence_L1;

  // Count the number of covered bases, and find the first and last element of the LCSk.
  int64_t indexfirst = -1;
  int64_t indexlast = -1;

  int64_t covered_bases = 0;
  int64_t covered_bases_query = 0, covered_bases_reference = 0;
  int64_t num_covering_kmers = 0;
  LOG_DEBUG_SPEC("Counting the covered bases and finding the first and the last brick index.\n");
  for (uint64_t i = 0; i < cluster_indices.size(); i++) {
    covered_bases_query += local_score->get_registry_entries().covered_bases_queries[cluster_indices[i]];
    covered_bases_reference += local_score->get_registry_entries().covered_bases_references[cluster_indices[i]];
    num_covering_kmers += local_score->get_registry_entries().num_kmers[cluster_indices[i]];
  }
  covered_bases = std::max(covered_bases_query, covered_bases_reference);

  if (cluster_indices.size() > 0) {
    indexfirst = cluster_indices.back();
    indexlast = cluster_indices.front();
  }

  // There are no valid graph paths! All scores were dismissed because of high deviation.
  // This is most likely a false positive.
  if (indexfirst == -1 || indexlast == -1) {
    LOG_DEBUG_SPEC("An error occured, indexfirst = %ld, indexlast = %ld\n", indexfirst, indexlast);
    return 1;
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

  #ifndef RELEASE_VERSION
    LOG_DEBUG_SPEC("Clusters:\n");
  #endif

  for (int64_t i=0; i<clusters.size(); i++) {
    if (clusters[i]) {
        Cluster mapping_cluster;
        mapping_cluster.query = clusters[i]->query;
        mapping_cluster.ref = clusters[i]->ref;
        mapping_info.clusters.push_back(mapping_cluster);
      #ifndef RELEASE_VERSION
        int64_t reference_start = indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id];
        int64_t region_start = local_score->get_region().start;
        LOG_DEBUG_SPEC("start(%ld, %ld), end(%ld, %ld)\tstart(%ld, %ld), end(%ld, %ld)\n",
                       mapping_cluster.query.start, mapping_cluster.ref.start, mapping_cluster.query.end, mapping_cluster.ref.end,
                       mapping_cluster.query.start, mapping_cluster.ref.start - reference_start, mapping_cluster.query.end, mapping_cluster.ref.end - reference_start);
      #endif
    }
  }



  L1Results l1_info;
  l1_info.l1_l = l;
  l1_info.l1_k = 1.0f;
  l1_info.l1_lmin = l - l_diff;
  l1_info.l1_lmax = l + l_diff;
  l1_info.l1_confidence_abs = confidence_L1;
  l1_info.l1_std = sigma_L2;
  l1_info.l1_rough_start = l1_info.l1_k * 0 + l1_info.l1_lmin;
  l1_info.l1_rough_end = l1_info.l1_k * read->get_sequence_length() + l1_info.l1_lmax;
  if (l1_info.l1_rough_start < indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id])
    l1_info.l1_rough_start = indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id];
  if (l1_info.l1_rough_end >= (indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id] + indexes[0]->get_reference_lengths()[local_score->get_region().reference_id]))
    l1_info.l1_rough_end = (indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id] + indexes[0]->get_reference_lengths()[local_score->get_region().reference_id]) - 1;

  mapping_info.is_mapped = true;

  PathGraphEntry *new_entry = new PathGraphEntry(indexes[0], read, parameters, (Region &) local_score->get_region(), &mapping_info, &l1_info);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
  mapping_data->intermediate_mappings.push_back(new_entry);

  #ifndef RELEASE_VERSION
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LOG_DEBUG_SPEC("Writing all anchors to file scores-%ld.\n",  local_score->get_scores_id());
      VerboseLocalScoresToFile(FormatString("temp/local_scores/scores-%ld.csv", local_score->get_scores_id()), read, local_score, NULL, 0, 0, false);
      LOG_DEBUG_SPEC("Writing LCSk anchors to file LCS-%ld.\n",  local_score->get_scores_id());
      VerboseLocalScoresToFile(FormatString("temp/local_scores/LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, 0, 0, false, NULL);
      LOG_DEBUG_SPEC("Writing cluster anchors to file LCSL1-%ld.\n",  local_score->get_scores_id());
      VerboseLocalScoresToFile(FormatString("temp/local_scores/LCSL1-%ld.csv", local_score->get_scores_id()), read, local_score, &second_filtered_lcskpp_indices, 0, 0, false);
      LOG_DEBUG_SPEC("Writing cluster anchors (again) to file double_LCS-%ld.\n",  local_score->get_scores_id());
      VerboseLocalScoresToFile(FormatString("temp/local_scores/double_LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &cluster_indices, l, 3.0f * confidence_L1, false, &cluster_ids);

      LOG_DEBUG_SPEC("LCSk clusters:\n");
      for (int64_t i=0; i<clusters.size(); i++) {
        LOG_DEBUG_SPEC_NO_HEADER("[%ld] num_anchors: %ld, length: %ld, coverage: %ld, query.start = %ld, query.end = %ld\n", i, clusters[i]->num_anchors, (clusters[i]->query.end - clusters[i]->query.start), clusters[i]->coverage, clusters[i]->query.start, clusters[i]->query.end);
      }
      LOG_DEBUG_SPEC_NEWLINE;
      LOG_DEBUG_SPEC("mapping_info.clusters:\n");
      for (int64_t i=0; i<mapping_info.clusters.size(); i++) {
        LOG_DEBUG_SPEC_NO_HEADER("[%ld] query.start = %ld, query.end = %ld, ref.start = %ld, ref.end = %ld\n", i, mapping_info.clusters[i].query.start, mapping_info.clusters[i].query.end, mapping_info.clusters[i].ref.start, mapping_info.clusters[i].ref.end);
      }
      LOG_DEBUG_SPEC_NEWLINE;
      VerboseClustersToFile_(FormatString("temp/clusters/clusters-read-%ld.csv", read->get_sequence_id()), local_score, mapping_data, indexes, read, parameters, clusters);
    }
    LOG_DEBUG_SPEC("Exiting function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024));
  #endif

  for (int64_t i=0; i<clusters.size(); i++) {
    if (clusters[i]) {
      delete clusters[i];
    }
  }

  return 0;
}
