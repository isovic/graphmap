/*
 * experimental.cc
 *
 *  Created on: May 1, 2015
 *      Author: isovic
 */

#include "graphmap/graphmap.h"
#include "graphmap/filter_anchors.h"

int GraphMap::AnchoredPostProcessRegionWithLCS_(ScoreRegistry* local_score, MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  LOG_DEBUG_SPEC("Entering function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB] current_readid = %ld, current_local_score = %ld\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024), read->get_sequence_id(), local_score->get_scores_id());
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

  // Parameters for anchor filtering.
  double indel_bandwidth_margin = parameters->error_rate/2 + 0.01f;
  int32_t max_dist = 200;
  int64_t min_covered_bases = 50; // TODO: need to experiment with this: std::min((int64_t) (read->get_sequence_length() * 0.10f), (int64_t) 50);
  int64_t cluster_size_cutoff = 2; // TODO: need to experiment with this. 1; // 2
  FilterAnchorsByChaining(read, local_score, parameters, lcskpp_indices, indel_bandwidth_margin, max_dist, 0, min_covered_bases, cluster_size_cutoff, first_filtered_lcskpp_indices, NULL);
  FilterAnchorsByDiff(read, local_score, parameters, first_filtered_lcskpp_indices, second_filtered_lcskpp_indices);
//  FilterAnchorsByChaining(read, local_score, parameters, lcskpp_indices, parameters->error_rate/2 + 0.01f, 200.0f, 0, 50, 2, first_filtered_lcskpp_indices, NULL);
//  FilterAnchorsByDiff(read, local_score, parameters, first_filtered_lcskpp_indices, second_filtered_lcskpp_indices);

  // Parameters for cluster filtering.
  int64_t min_num_anchors_in_cluster = 2; // TODO: need to experiment with this. 1; // 2;
  int64_t min_cluster_length = min_covered_bases; // 50;
  int64_t min_cluster_covered_bases = min_cluster_length; // 50;
  float min_cluster_coverage = 0.0f;
  GenerateClusters(min_num_anchors_in_cluster, min_cluster_length, min_cluster_covered_bases, min_cluster_coverage,
                   second_filtered_lcskpp_indices, local_score, mapping_data, read, parameters, clusters, cluster_indices, &cluster_ids);
//  GenerateClusters(2, 50, 50, 0.0, second_filtered_lcskpp_indices, local_score, mapping_data, indexes, read, parameters, clusters, cluster_indices, &cluster_ids);

  // Find the L1 parameters (median line and the confidence intervals).
  float l_diff = read->get_sequence_length() * parameters->error_rate;
  float maximum_allowed_deviation = l_diff * sqrt(2.0f) / 2.0f;
  float sigma_L2 = 0.0f, confidence_L1 = 0.0f;
  int64_t k = 0, l = 0;
  // Actual L1 calculation.
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
        mapping_cluster.coverage = clusters[i]->coverage;
        mapping_cluster.num_anchors = clusters[i]->num_anchors;
        mapping_cluster.valid = true;
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
//      VerboseClustersToFile_(FormatString("temp/clusters/clusters-read-%ld.csv", read->get_sequence_id()), local_score, mapping_data, indexes, read, parameters, clusters);
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
