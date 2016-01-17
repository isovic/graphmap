/*
 * experimental.cc
 *
 *  Created on: May 1, 2015
 *      Author: isovic
 */

#include "graphmap/graphmap.h"

/// These are some constants used for filtering shady anchors.
/// TODO: This can be omitted if dynamic programming was used to penalize the anchor distances.
/// int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * MIN_CLUSTER_COVERAGE_FACTOR;
#define MIN_CLUSTER_COVERAGE_FACTOR 0.05f
/// int64_t min_cluster_length = read->get_sequence_length() * MIN_CLUSTER_LENGTH_FACTOR;
#define MIN_CLUSTER_LENGTH_FACTOR 0.03f



struct ClusterAndIndices {
  Range query;
  Range ref;
  int32_t num_anchors = 0;
  int32_t coverage = 0;
  std::vector<int> lcskpp_indices;
};

bool CheckDistanceTooBig(const Vertices& registry_entries, int64_t index_last, int64_t index_current, const ProgramParameters* parameters) {
  int64_t distance_query = registry_entries.query_ends[index_current] - registry_entries.query_starts[index_last];
  int64_t distance_ref = registry_entries.reference_ends[index_current] - registry_entries.reference_starts[index_last];
  float max_length = ((float) std::max(distance_query, distance_ref));
  float min_length = ((float) std::min(distance_query, distance_ref));
  if ((min_length == 0 && max_length != 0) || (min_length > 0 && (max_length / min_length - 1.0f) > parameters->error_rate / 2.0f)) {
    return true;
  }

  return false;
}

bool CheckDistanceStep(const Vertices& registry_entries, int64_t index_first, int64_t index_last, int64_t index_current, float max_diff) { // , const ProgramParameters* parameters) {
  int64_t seed_length = 12;
  int64_t distance_query = registry_entries.query_ends[index_current] - registry_entries.query_starts[index_first];
  int64_t distance_query_before = registry_entries.query_ends[index_last] - registry_entries.query_starts[index_first];

  int64_t distance_ref = registry_entries.reference_ends[index_current] - registry_entries.reference_starts[index_first];
  int64_t distance_ref_before = registry_entries.reference_ends[index_last] - registry_entries.reference_starts[index_first];

  float diff_query = ((float) distance_query) / ((float) distance_query_before);
  float diff_ref = ((float) distance_ref) / ((float) distance_ref_before);

  if (diff_query > max_diff || diff_ref > max_diff)
    return true;

  return false;
}

int FilterAnchorBreakpoints(int64_t min_cluster_length, float min_cluster_coverage, std::vector<int> &lcskpp_indices, ScoreRegistry* local_score, MappingData* mapping_data, const std::vector<Index *> indexes, const SingleSequence* read, const ProgramParameters* parameters, std::vector<ClusterAndIndices *> &ret_clusters, std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {
  std::vector<ClusterAndIndices *> clusters;

  ClusterAndIndices *new_cluster = NULL;
  int64_t last_nonskipped_i = lcskpp_indices.size() + 1;
  for (int64_t i=(lcskpp_indices.size() - 1); i >= 0; i--) {
//  for (int64_t i=0; i<(lcskpp_indices.size()); i++) {
    /// Skip anchors which might be too erroneous.
    int64_t current_lcskp_index = lcskpp_indices.at(i);
    if (CheckDistanceTooBig(local_score->get_registry_entries(), current_lcskp_index, current_lcskp_index, parameters) == true) {
      continue;
    }

    if (new_cluster == NULL || last_nonskipped_i > lcskpp_indices.size()) {

    } else {
      /// This is going to work, because last_nonskipped_i will be set the second iteration of the loop. The value of i starts counting from int64_t i=(lcskpp_indices.size() - 1).
      int64_t previous_lcskp_index = lcskpp_indices.at(last_nonskipped_i);

      bool wrong_to_previous1 = CheckDistanceTooBig(local_score->get_registry_entries(), previous_lcskp_index, current_lcskp_index, parameters);
      bool wrong_to_previous2 = (new_cluster->lcskpp_indices.size() < 2) ? false :
                                (CheckDistanceTooBig(local_score->get_registry_entries(), new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-2], current_lcskp_index, parameters));
      bool wrong_by_distance = false; // CheckDistanceStep(local_score->get_registry_entries(), new_cluster->lcskpp_indices.front(), previous_lcskp_index, current_lcskp_index, 1.5f);
      if ((wrong_to_previous1 == true && wrong_to_previous2 == true) || wrong_by_distance == true) {
        /// In this case, the new point is a general outlier to the previous LCSk, because it doesn't fit neither to the previous point, nor to the point before that.
        /// Check if the currently collected cluster satisfies the conditions for it to be taken into account.
        if (new_cluster != NULL) {
          int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;

          if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {

            clusters.push_back(new_cluster);
            /// TODO: The lines with a + sign was added to the master branch to dirty-fix a bug with anchored alignments. The bug was that an anchor might have overlapped a previous anchor, because of wrong endpoints
            /// This should have been solved on the dev branch earlier than this, but for some reasons I couldn't have merged the branches, and introduced a divergence.
            /// However, it may prove that this fix was actually useful, so I left it here just in case.
            //  +           if (clusters.size() > 0 && (new_cluster->query.start <= clusters.back()->query.end || new_cluster->ref.start <= clusters.back()->ref.end)) {
            //  +               delete new_cluster;
            //  +               new_cluster = NULL;
            //  +           } else {
            //  +                 clusters.push_back(new_cluster);
            //  +                 new_cluster = NULL;
            //  +           }
            //  +

          } else {
//            printf ("Deleted [Cluster %ld]\n", ret_clusters.size());
//            printf ("(new_cluster->query.end - new_cluster->query.start + 1) = %ld\n", (new_cluster->query.end - new_cluster->query.start + 1));
//            printf ("min_cluster_length = %ld\n", min_cluster_length);
//            printf ("new_cluster->coverage = %ld\n", new_cluster->coverage);
//            printf ("min_covered_bases = %ld\n", min_covered_bases);
//            printf ("min_cluster_coverage = %f\n", min_cluster_coverage);
//            fflush(stdout);
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[Deleted Cluster %ld] (a)\n", clusters.size()), std::string(__FUNCTION__));
//            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read,
//                                         FormatString("(new_cluster->query.end - new_cluster->query.start + 1) = %ld\n", (new_cluster->query.end - new_cluster->query.start + 1)), std::string(__FUNCTION__));
//            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read,
//                                         FormatString("min_cluster_length = %ld\n", min_cluster_length), std::string(__FUNCTION__));
//            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read,
//                                         FormatString("new_cluster->coverage = %ld\n", new_cluster->coverage), std::string(__FUNCTION__));
//            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read,
//                                         FormatString("min_covered_bases = %ld\n", min_covered_bases), std::string(__FUNCTION__));
//            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read,
//                                         FormatString("min_cluster_coverage = %f\n", min_cluster_coverage), std::string(__FUNCTION__));

            delete new_cluster;
          }
          new_cluster = NULL;
        }
      } else if (wrong_to_previous1 == true && wrong_to_previous2 == false) {
        /// In this case, the previous point was an outlier, because the new point fits better to the one before the previous one. Overwrite the previous entry in new_cluster.
        new_cluster->query.end = local_score->get_registry_entries().query_ends[current_lcskp_index] - 1;
        new_cluster->ref.end = local_score->get_registry_entries().reference_ends[current_lcskp_index] - 1;
        new_cluster->coverage -= local_score->get_registry_entries().covered_bases_queries[previous_lcskp_index];
        new_cluster->coverage += local_score->get_registry_entries().covered_bases_queries[current_lcskp_index];
        new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-1] = current_lcskp_index;

        if (new_cluster->lcskpp_indices.size() == 1) {
          new_cluster->query.start = local_score->get_registry_entries().query_starts[current_lcskp_index];
          new_cluster->ref.start = local_score->get_registry_entries().reference_starts[current_lcskp_index];
        }
        last_nonskipped_i = i;

        continue;
      }
    }

    /// TODO: This line with a + sign was added to the master branch to dirty-fix a bug with anchored alignments. The bug was that an anchor might have overlapped a previous anchor, because of wrong endpoints
    /// This should have been solved on the dev branch earlier than this, but for some reasons I couldn't have merged the branches, and introduced a divergence.
    /// However, it may prove that this fix was actually useful, so I left it here just in case.
    ///  +  if (clusters.size() == 0 || (clusters.size() > 0 && (local_score->get_registry_entries().query_starts[current_lcskp_index] > clusters.back()->query.end && local_score->get_registry_entries().reference_starts[current_lcskp_index] > c kp_index] > clusters.back()->query.end && local_score->get_registry_entries().reference_starts[current_lcskp_index] > clusters.back()->ref.end))) {
    if (new_cluster == NULL) {
      new_cluster = new ClusterAndIndices;
      new_cluster->query.start = local_score->get_registry_entries().query_starts[current_lcskp_index];
      new_cluster->ref.start = local_score->get_registry_entries().reference_starts[current_lcskp_index];
    }
    new_cluster->query.end = local_score->get_registry_entries().query_ends[current_lcskp_index] - 1;
    new_cluster->ref.end = local_score->get_registry_entries().reference_ends[current_lcskp_index] - 1;
    new_cluster->num_anchors += 1;
    new_cluster->coverage += local_score->get_registry_entries().covered_bases_queries[current_lcskp_index];
    new_cluster->lcskpp_indices.push_back(current_lcskp_index);

//    printf ("i = %ld\n", i);
//    printf ("new_cluster->query.start = %ld\n", new_cluster->query.start);
//    printf ("new_cluster->query.end = %ld\n", new_cluster->query.end);
//    printf ("\n");
//    fflush(stdout);

    last_nonskipped_i = i;
    ///  +  }

  }
  /// Push the last cluster.
  if (new_cluster != NULL) {
    int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;
    if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {
      if (clusters.size() > 0 && (new_cluster->query.start <= clusters.back()->query.end || new_cluster->ref.start <= clusters.back()->ref.end)) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[Deleted Cluster %ld] (a)\n", clusters.size()), std::string(__FUNCTION__));

        delete new_cluster;
        new_cluster = NULL;
      } else {
        clusters.push_back(new_cluster);
        new_cluster = NULL;
      }
    } else {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[Deleted Cluster %ld] (a)\n", clusters.size()), std::string(__FUNCTION__));

      delete new_cluster;
      new_cluster = NULL;
    }
    new_cluster = NULL;
  }

  ret_filtered_lcskpp_indices.clear();
  if (ret_cluster_ids) {
    ret_cluster_ids->clear();
  }

  //////////////////////
  /// Filter shady anchors by thresholding the coverage by some statistics on the clusters.
  /// Concretely, for each cluster the percentage of covered bases is calculated. Using the median and the standard deviation
  /// we set the threshold for the minimum number of covered bases in a cluster.
  /// If below this threshold, the cluster is dismissed.
  std::vector<float> cluster_coverages;
  for (int64_t i=0; i<clusters.size(); i++) {
    int64_t cluster_length = clusters[i]->query.end - clusters[i]->query.start + 1;
    int64_t covered_bases = clusters[i]->coverage;
    cluster_coverages.push_back(((float) covered_bases) / ((float) cluster_length));
  }
  std::sort(cluster_coverages.begin(), cluster_coverages.end());
  float median = 0.0f;
  int64_t num_clusters = cluster_coverages.size();
  if (num_clusters > 0) {
    median = ((num_clusters % 2) == 0) ? ((cluster_coverages[num_clusters/2] + cluster_coverages[num_clusters/2+1]) / 2.0f) :
                                         cluster_coverages[(int64_t) floor(((float) num_clusters) / 2.0f)];
  }
  float mean = 0.0f;
  for (int64_t i=0; i<cluster_coverages.size(); i++) {
    mean += cluster_coverages[i];
  }
  if (cluster_coverages.size() > 0)
    mean /= cluster_coverages.size();
  float std = 0.0f, std_med = 0.0f;
  for (int64_t i=0; i<cluster_coverages.size(); i++) {
    std += (cluster_coverages[i] - mean) * (cluster_coverages[i] - mean);
    std_med += (cluster_coverages[i] - median) * (cluster_coverages[i] - median);
  }
  if (cluster_coverages.size() > 1) {
    std /= (cluster_coverages.size() - 1);
    std_med /= (cluster_coverages.size() - 1);
  }
  std = sqrt(std);
  std_med = sqrt(std_med);
  float min_relative_cluster_coverage = median - std_med - 0.001;  /// The 0.001 is subtracted to avoid the numerical error.

  if (min_relative_cluster_coverage <= 0.0f)
    min_relative_cluster_coverage = 0.001f;



  /// Generate the final ret_clusters that will be returned and used further.
  for (int64_t i=0; i<clusters.size(); i++) {
    int64_t cluster_length = clusters[i]->query.end - clusters[i]->query.start + 1;
    int64_t covered_bases = clusters[i]->coverage;
    float cluster_coverage = ((float) covered_bases) / ((float) cluster_length);

    if (cluster_coverage >= min_relative_cluster_coverage) {
      ret_filtered_lcskpp_indices.insert(ret_filtered_lcskpp_indices.end(), clusters[i]->lcskpp_indices.begin(), clusters[i]->lcskpp_indices.end());

      /// Create indices for debugging purposes (so we can differentiate clusters).
      if (ret_cluster_ids) {
        std::vector<int32_t> cluster_indices(clusters[i]->lcskpp_indices.size(), i);
        ret_cluster_ids->insert(ret_cluster_ids->end(), cluster_indices.begin(), cluster_indices.end());
      }
      ret_clusters.push_back(clusters[i]);

    } else {
      delete clusters[i];
    }

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[Cluster %ld] cluster_length = %ld, covered_bases = %ld\n", i, cluster_length, covered_bases), "L1-PostProcessRegionWithLCS_");
  }

  return 0;
}

int GraphMap::AnchoredPostProcessRegionWithLCS_(ScoreRegistry* local_score, MappingData* mapping_data, const std::vector<Index *> indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || ((int64_t) read->get_sequence_id()) == parameters->debug_read), FormatString("Entering function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB] current_readid = %ld, current_local_score = %ld\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024), read->get_sequence_id(), local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS");
  int lcskpp_length = 0;
  std::vector<int> lcskpp_indices;
  CalcLCSFromLocalScoresCacheFriendly_(&(local_score->get_registry_entries()), false, 0, 0, &lcskpp_length, &lcskpp_indices);
  if (lcskpp_length == 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Current local scores: %ld, lcskpp_length == 0 || best_score == NULL\n", local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS");
    return 1;
  }

  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("After LCSk:\n", local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
    for (int64_t i = 0; i < lcskpp_indices.size(); i++) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld] %s\n", i, local_score->get_registry_entries().VerboseToString(lcskpp_indices[i]).c_str()), "[]");
    }
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
  }

  /// Filter the LCSk anchors.
  std::vector<ClusterAndIndices *> clusters;
  std::vector<int> cluster_indices;
  std::vector<int32_t> cluster_ids;
  FilterAnchorBreakpoints(read->get_sequence_length() * MIN_CLUSTER_LENGTH_FACTOR, MIN_CLUSTER_COVERAGE_FACTOR, lcskpp_indices, local_score, mapping_data, indexes, read, parameters, clusters, cluster_indices, &cluster_ids);

  // Find the L1 parameters (median line and the confidence intervals).
  float l_diff = read->get_sequence_length() * parameters->error_rate;
  float maximum_allowed_deviation = l_diff * sqrt(2.0f) / 2.0f;
  float sigma_L2 = 0.0f, confidence_L1 = 0.0f;
  int64_t k = 0, l = 0;
  // Actuall L1 calculation.
  int ret_L1 = CalculateL1ParametersWithMaximumDeviation_(local_score, cluster_indices, maximum_allowed_deviation, &k, &l, &sigma_L2, &confidence_L1);
  // Sanity check.
  if (ret_L1) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("An error occured, L1 function (I) returned with %d!\n", ret_L1), "L1-PostProcessRegionWithLCS_");

    #ifndef RELEASE_VERSION
      if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing all anchors to file scores-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
        VerboseLocalScoresToFile(FormatString("temp/local_scores/scores-%ld.csv", local_score->get_scores_id()), read, local_score, NULL, 0, 0, false);
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing LCSk anchors to file LCS-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
        VerboseLocalScoresToFile(FormatString("temp/local_scores/LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, 0, 0, false);
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing LCSk anchors to file LCS-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
        VerboseLocalScoresToFile(FormatString("temp/local_scores/LCSL1-%ld.csv", local_score->get_scores_id()), read, local_score, &cluster_indices, 0, 0, false);
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing LCSk anchors to file LCS-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
        VerboseLocalScoresToFile(FormatString("temp/local_scores/double_LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, 0, 0, false);
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
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Counting the covered bases and finding the first and the last brick index.\n"), "PostProcessRegionWithLCS_-DoubleLCSk");
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
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("An error occured, indexfirst = %ld, indexlast = %ld\n", indexfirst, indexlast), "L1-PostProcessRegionWithLCS_");
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
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Clusters:\n"), "ExperimentalPostProcessRegionWithLCS_");
#endif

  for (int64_t i=0; i<clusters.size(); i++) {
    if (clusters[i]) {
//      int64_t cluster_length = clusters[i]->query.end - clusters[i]->query.start + 1;
//      int64_t covered_bases = clusters[i]->coverage;
//      if (cluster_length >= min_cluster_length && covered_bases >= min_covered_bases) {
        Cluster mapping_cluster;
        mapping_cluster.query = clusters[i]->query;
        mapping_cluster.ref = clusters[i]->ref;
        mapping_info.clusters.push_back(mapping_cluster);
#ifndef RELEASE_VERSION
      int64_t reference_start = indexes[0]->get_reference_starting_pos()[local_score->get_region().reference_id];
      int64_t region_start = local_score->get_region().start;
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("start(%ld, %ld), end(%ld, %ld)\tstart(%ld, %ld), end(%ld, %ld)\n", mapping_cluster.query.start, mapping_cluster.ref.start, mapping_cluster.query.end, mapping_cluster.ref.end,
                                                                                                                                    mapping_cluster.query.start, mapping_cluster.ref.start - reference_start, mapping_cluster.query.end, mapping_cluster.ref.end - reference_start), "");
#endif
//      }
//      delete clusters[i];
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

//  CheckMinimumMappingConditions_(&mapping_info, &l1_info, index, read, parameters);

  mapping_info.is_mapped = true;

  PathGraphEntry *new_entry = new PathGraphEntry(indexes[0], read, parameters, (Region &) local_score->get_region(), &mapping_info, &l1_info);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
  mapping_data->intermediate_mappings.push_back(new_entry);

#ifndef RELEASE_VERSION
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing all anchors to file scores-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
    VerboseLocalScoresToFile(FormatString("temp/local_scores/scores-%ld.csv", local_score->get_scores_id()), read, local_score, NULL, 0, 0, false);

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing LCSk anchors to file LCS-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
    VerboseLocalScoresToFile(FormatString("temp/local_scores/LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, 0, 0, false);

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing cluster anchors to file LCSL1-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
    VerboseLocalScoresToFile(FormatString("temp/local_scores/LCSL1-%ld.csv", local_score->get_scores_id()), read, local_score, &cluster_indices, 0, 0, false, &cluster_ids);

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Writing cluster anchors (again) to file double_LCS-%ld.\n",  local_score->get_scores_id()), "ExperimentalPostProcessRegionWithLCS_");
    VerboseLocalScoresToFile(FormatString("temp/local_scores/double_LCS-%ld.csv", local_score->get_scores_id()), read, local_score, &lcskpp_indices, l, 3.0f * confidence_L1, true);

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("LCSk clusters:\n"), "ExperimentalPostProcessRegionWithLCS_");
    for (int64_t i=0; i<clusters.size(); i++) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld] num_anchors: %ld, length: %ld, coverage: %ld, query.start = %ld, query.end = %ld\n", i, clusters[i]->num_anchors, (clusters[i]->query.end - clusters[i]->query.start), clusters[i]->coverage, clusters[i]->query.start, clusters[i]->query.end), "[]");
    }
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("mapping_info.clusters:\n"), "ExperimentalPostProcessRegionWithLCS_");
    for (int64_t i=0; i<mapping_info.clusters.size(); i++) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld] query.start = %ld, query.end = %ld, ref.start = %ld, ref.end = %ld\n", i, mapping_info.clusters[i].query.start, mapping_info.clusters[i].query.end, mapping_info.clusters[i].ref.start, mapping_info.clusters[i].ref.end), "[]");
    }
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || read->get_sequence_id() == parameters->debug_read), FormatString("Exiting function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "PostProcessRegionWithLCS_");
#endif

  for (int64_t i=0; i<clusters.size(); i++) {
    if (clusters[i]) {
      delete clusters[i];
    }
  }

  return 0;
}

int AnchoredAlignment(bool is_linear, bool end_to_end, AlignmentFunctionType AlignmentFunctionNW, AlignmentFunctionType AlignmentFunctionSHW, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                           int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part, int64_t *ret_nonclipped_left_part,
                           int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part, int64_t *ret_nonclipped_right_part,
                           SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                           int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op) {

  if (best_path->get_mapping_data().clusters.size() <= 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "No valid anchors exist!", "LocalRealignmentCircular");
    return ALIGNMENT_WRONG_CLUSTER_SIZE;
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Entering anchored alignment.\n", "AnchoredAlignment");
  int8_t *ref_data = (int8_t *) index->get_data();
  int64_t region_length_joined = 0, start_offset = 0, position_of_ref_end = 0;

  int64_t absolute_reference_id = best_path->get_region_data().reference_id;
  int64_t reference_id = best_path->get_region_data().reference_id;
  int64_t reference_start = index->get_reference_starting_pos()[absolute_reference_id];
  int64_t reference_length = index->get_reference_lengths()[absolute_reference_id];
  SeqOrientation orientation = (best_path->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);


  if (is_linear == false) {
    if (best_path->get_region_data().is_split == false) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Called the function for handling the circular part of the genome, but alignment is not split. best_path->region.is_split == false.\n\n"), "LocalRealignmentCircular");
      return ALIGNMENT_NOT_CIRCULAR;
    }

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Concatenating regions for circular alignment.\n", "AnchoredAlignment");
    ConcatenateSplitRegion(index, (Region &) best_path->get_region_data(), &ref_data, &region_length_joined, &start_offset, &position_of_ref_end);
    reference_start = 0;
    reference_length = region_length_joined;
  }




  int64_t edit_distance = 0;
  std::vector<unsigned char> alignment;

  int64_t clip_count_front = best_path->get_mapping_data().clusters.front().query.start;
  int64_t clip_count_back = read->get_sequence_length() - (best_path->get_mapping_data().clusters.back().query.end) - 1;  /////I

  int64_t alignment_position_start = best_path->get_mapping_data().clusters.front().ref.start; // - clip_count_front;
  int64_t alignment_position_end = best_path->get_mapping_data().clusters.back().ref.end; // + clip_count_back;  /////I

  int64_t query_start = best_path->get_mapping_data().clusters.front().query.start; // - clip_count_front;
  int64_t query_end = best_path->get_mapping_data().clusters.back().query.end; // + clip_count_back;  /////I

  /// Aligning the begining of the read (in front of the first anchor).
  if (clip_count_front > 0) {
    /// Check if we need to extend the alignment to the left boundary. Also, even if the user specified it, if we are to close to the boundary, just clip it.
    if (end_to_end == false || ((alignment_position_start - clip_count_front*2) < reference_start)) {
      std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
      alignment.insert(alignment.begin(), insertions_front.begin(), insertions_front.end());

    } else {
      if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                            "Aligning the begining of the read (overhang).\n", "LocalRealignmentLinear");
      }

      /// Reversing the sequences to make the semiglobal alignment of the trailing and leading parts.
      int8_t *reversed_query_front = reverse_data(read->get_data(), clip_count_front);
      int8_t *reversed_ref_front = reverse_data(ref_data + (alignment_position_start - 1) - (clip_count_front*2 - 1), clip_count_front*2);

      int64_t leftover_left_start = 0, leftover_left_end = 0, leftover_left_edit_distance = 0;
      std::vector<unsigned char> leftover_left_alignment;
      int ret_code_right = AlignmentFunctionSHW(reversed_query_front, (clip_count_front),
                                       (int8_t *) (reversed_ref_front), (clip_count_front*2),
                                       -1, parameters.match_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                       &leftover_left_start, &leftover_left_end,
                                       &leftover_left_edit_distance, leftover_left_alignment);
      /// Check if the return code is ok. Otherwise, just clip the front.
      /// Added on 15.11.2015.: check if the edit distance of the front part is too high. EDlib will automatically return an error code, but SeqAn won't.
      /// An example is when the entire front part does not match (e.g. alignment of a read to a part of the reference consisted of N bases).
      if (ret_code_right != 0 || leftover_left_edit_distance > clip_count_front/2) {
        // TODO: This is a nasty hack. EDlib used to crash when query and target are extremely small, e.g. query = "C" and target = "TC".
        // In this manner we just ignore the leading part, and clip it.
        std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
        alignment.insert(alignment.begin(), insertions_front.begin(), insertions_front.end());

      } else {
        if (leftover_left_alignment.size() == 0) {
          std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
          alignment.insert(alignment.begin(), insertions_front.begin(), insertions_front.end());
        } else {
          unsigned char *reversed_alignment = reverse_data(&(leftover_left_alignment[0]), leftover_left_alignment.size());
          alignment.insert(alignment.begin(), reversed_alignment, reversed_alignment + leftover_left_alignment.size());
          alignment_position_start -= leftover_left_end + 1;
          if (reversed_alignment)
            free(reversed_alignment);
        }

        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
          std::string alignment_as_string = "";

          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    FormatString("End of the beginning part of the read: %ld\n", (alignment_position_start - 1) - (clip_count_front*2 - 1)), "[]");
          alignment_as_string = PrintAlignmentToString((const unsigned char *) reversed_query_front, clip_count_front,
                                                       (const unsigned char *) (reversed_ref_front), clip_count_front*2,
                                                       (unsigned char *) &(leftover_left_alignment[0]), leftover_left_alignment.size(),
                                                       (0), MYERS_MODE_SHW);
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    FormatString("Aligning the beginning of the read:\n%s\n", alignment_as_string.c_str()), "[]");
        }

        if (reversed_query_front)
          free(reversed_query_front);
        if (reversed_ref_front)
          free(reversed_ref_front);
      }

    }
  }

  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
    for (int64_t i=0; i<best_path->get_mapping_data().clusters.size(); i++) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("[anchor %d] [%ld, %ld]-[%ld, %ld]\n", i, best_path->get_mapping_data().clusters[i].query.start, best_path->get_mapping_data().clusters[i].ref.start,
                                                       best_path->get_mapping_data().clusters[i].query.end, best_path->get_mapping_data().clusters[i].ref.end), "[]");
    }
  }

  for (int64_t i=0; i<best_path->get_mapping_data().clusters.size(); i++) {
    if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          "Aligning an anchor.\n", "LocalRealignmentLinear");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("[anchor %d] [%ld, %ld]-[%ld, %ld]\n", i, best_path->get_mapping_data().clusters[i].query.start, best_path->get_mapping_data().clusters[i].ref.start,
                                                       best_path->get_mapping_data().clusters[i].query.end, best_path->get_mapping_data().clusters[i].ref.end), "[]");
    }

    ///////////////////////////
    /// Align the anchor.
    int64_t query_start = best_path->get_mapping_data().clusters[i].query.start;
    int64_t query_end = best_path->get_mapping_data().clusters[i].query.end;
    int64_t ref_start = best_path->get_mapping_data().clusters[i].ref.start;
    int64_t ref_end = best_path->get_mapping_data().clusters[i].ref.end;
    int64_t query_alignment_length = query_end - query_start + 1;
    int64_t ref_alignment_length = ref_end - ref_start + 1;

    int64_t anchor_alignment_position_start = 0, anchor_alignment_position_end = 0, anchor_edit_distance = 0;
    std::vector<unsigned char> anchor_alignment;
    int ret_code1 = AlignmentFunctionNW(read->get_data() + query_start, (query_alignment_length),
                                        (int8_t *) (ref_data + ref_start), (ref_alignment_length),
                                     -1, parameters.match_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                     &anchor_alignment_position_start, &anchor_alignment_position_end,
                                     &anchor_edit_distance, anchor_alignment);
    if (ret_code1 != 0 || anchor_alignment.size() == 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("Alignment returned with error! ret_code1 = %d\n", ret_code1), "LocalRealignmentLinear");
      return ret_code1;
    }

    if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
      std::string alignment_as_string = "";
      alignment_as_string = PrintAlignmentToString((const unsigned char *) (read->get_data() + query_start), query_alignment_length,
                                                   (const unsigned char *) (ref_data + ref_start), (ref_alignment_length),
                                                   (unsigned char *) &(anchor_alignment[0]), anchor_alignment.size(),
                                                   (0), MYERS_MODE_NW);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                FormatString("Aligned anchor %d:\n%s\n", i, alignment_as_string.c_str()), "[]");
    }

    edit_distance += anchor_edit_distance;
    /// Check for a special case when previous global alignment ended with deletions or insertions, and the new one starts with deletions or insertions.
    /// Switching from deletions to insertions is basically a mismatch streak.
    if (alignment.size() > 0 && anchor_alignment.size() > 0 && ((alignment.back() == EDLIB_D && anchor_alignment[0] == EDLIB_I) || (alignment.back() == EDLIB_I && anchor_alignment[0] == EDLIB_D))) {
      int64_t num_trailing_indels = 0;
      int64_t num_leading_indels = 0;
      int64_t current_op1 = alignment.size() - 1;
      while (current_op1 >= 0) {
        if ((current_op1 + 1) < alignment.size() && alignment[current_op1] != alignment[current_op1+1])
          break;
        num_trailing_indels += 1;
        current_op1 -= 1;
      }
      int64_t current_op2 = 0;
      while (current_op2 < anchor_alignment.size()) {
        if (current_op2 > 0 && anchor_alignment[current_op2] != anchor_alignment[current_op2-1])
          break;
        num_leading_indels += 1;
        current_op2 += 1;
      }

      int64_t min_count = std::min(num_trailing_indels, num_leading_indels);

      for (current_op1 = 0; current_op1 < min_count; current_op1++) {
        if ((ref_data + ref_start + anchor_alignment_position_start - current_op1 - 1) == (read->get_data() + (query_start) - current_op1))
          alignment[alignment.size() - current_op1 - 1] = EDLIB_EQUAL;
        else
          alignment[alignment.size() - current_op1 - 1] = EDLIB_X;
      }

      alignment.insert(alignment.end(), anchor_alignment.begin() + min_count, anchor_alignment.end());

    } else {
      alignment.insert(alignment.end(), anchor_alignment.begin(), anchor_alignment.end());
    }


    ///////////////////////////
    ///////////////////////////
    /// Align in between the anchors.
    if ((i + 1) < best_path->get_mapping_data().clusters.size()) {
      int64_t next_query_start = best_path->get_mapping_data().clusters[i+1].query.start;
      int64_t next_query_end = best_path->get_mapping_data().clusters[i+1].query.end;
      int64_t next_ref_start = best_path->get_mapping_data().clusters[i+1].ref.start;
      int64_t next_ref_end = best_path->get_mapping_data().clusters[i+1].ref.end;
      int64_t inbetween_query_length = (next_query_start - (query_end)) - 1;  /////I
      int64_t inbetween_ref_length = (next_ref_start - (ref_end)) - 1;  /////I

      /// Check if there is actually any distance between the queries, or between the references.
      /// If there is no difference, that means there is a clean insertion/deletion.
      if (inbetween_query_length <= 0 && inbetween_ref_length > 0) {
        std::vector<unsigned char> deletions_inbetween(inbetween_ref_length, EDLIB_D);
        alignment.insert(alignment.end(), deletions_inbetween.begin(), deletions_inbetween.end());

      } else if (inbetween_ref_length <= 0 && inbetween_query_length > 0) {
        std::vector<unsigned char> insertions_inbetween(inbetween_query_length, EDLIB_I);
        alignment.insert(alignment.end(), insertions_inbetween.begin(), insertions_inbetween.end());

      } else if (inbetween_query_length < 0 && inbetween_ref_length < 0) {
        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    "Problem aligning in between anchors!\n", "LocalRealignmentLinear");
        }
        return ALIGNMENT_DISTANCE_BETWEEN_ANCHORS_PROBLEM;

      } else if (inbetween_query_length > 0 && inbetween_ref_length > 0) {
        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                              "Aligning in between anchors.\n", "LocalRealignmentLinear");
        }

        int64_t between_alignment_position_start = 0, between_alignment_position_end = 0, between_anchor_edit_distance = 0;
        std::vector<unsigned char> between_anchor_alignment;
        int ret_code2 = AlignmentFunctionNW(read->get_data() + (query_end) + 1, inbetween_query_length,
                                            (int8_t *) (ref_data + ref_end) + 1, inbetween_ref_length,
                                         -1, parameters.match_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                         &between_alignment_position_start, &between_alignment_position_end,
                                         &between_anchor_edit_distance, between_anchor_alignment);

        if (ret_code2 != 0 || between_anchor_alignment.size() == 0) {
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                              FormatString("Alignment returned with error! ret_code2 = %d\n", ret_code2), "LocalRealignmentLinear");
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                              FormatString("inbetween_query_length = %ld\ninbetween_ref_length = %ld\nnext_ref_start = %ld\nref_end = %ld\n",
                                                           inbetween_query_length, inbetween_ref_length, next_ref_start, ref_end), "[]");
          return ret_code2;
        }
        edit_distance += between_anchor_edit_distance;

        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
          std::string alignment_as_string = "";
          alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + (query_end) + 1, inbetween_query_length,
                                                       (const unsigned char *) (ref_data + ref_end) + 1, inbetween_ref_length,
                                                       (unsigned char *) &(between_anchor_alignment[0]), between_anchor_alignment.size(),
                                                       (0), MYERS_MODE_NW);
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    FormatString("Aligning in between anchors %d and %d:\n%s\n", i, (i+1), alignment_as_string.c_str()), "[]");
        }

        /// Check for a special case when previous global alignment ended with deletions or insertions, and the new one starts with deletions or insertions.
        /// Switching from deletions to insertions is basically a mismatch streak.
        if (alignment.size() > 0 && between_anchor_alignment.size() > 0 && ((alignment.back() == EDLIB_D && between_anchor_alignment[0] == EDLIB_I) || (alignment.back() == EDLIB_I && between_anchor_alignment[0] == EDLIB_D))) {
          int64_t num_trailing_indels = 0;
          int64_t num_leading_indels = 0;
          int64_t current_op1 = alignment.size() - 1;
          while (current_op1 >= 0) {
            if ((current_op1 + 1) < alignment.size() && alignment[current_op1] != alignment[current_op1+1])
              break;
            num_trailing_indels += 1;
            current_op1 -= 1;
          }
          int64_t current_op2 = 0;
          while (current_op2 < between_anchor_alignment.size()) {
            if (current_op2 > 0 && between_anchor_alignment[current_op2] != between_anchor_alignment[current_op2-1])
              break;
            num_leading_indels += 1;
            current_op2 += 1;
          }

          int64_t min_count = std::min(num_trailing_indels, num_leading_indels);
          for (current_op1 = 0; current_op1 < min_count; current_op1++) {
            if ((ref_data + ref_end + between_alignment_position_start - current_op1 - 1) == (read->get_data() + (query_end) - current_op1))
              alignment[alignment.size() - current_op1 - 1] = EDLIB_EQUAL;
            else
              alignment[alignment.size() - current_op1 - 1] = EDLIB_X;
          }

          alignment.insert(alignment.end(), between_anchor_alignment.begin() + min_count, between_anchor_alignment.end());

        } else {
          alignment.insert(alignment.end(), between_anchor_alignment.begin(), between_anchor_alignment.end());
        }
      }
    }
  }



  /// Aligning the end of the read.
  if (clip_count_back > 0) {
    /// Handle the clipping at the end, or extend alignment to the end of the sequence.
    if (end_to_end == false || (alignment_position_end + 1 + clip_count_back * 2) >= (reference_start + reference_length)) {
        std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
        alignment.insert(alignment.end(), insertions_back.begin(), insertions_back.end());

    } else {
      if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                            "Aligning the end of the read (overhang).\n", "LocalRealignmentLinear");
      }

      int64_t leftover_right_start = 0, leftover_right_end = 0, leftover_right_edit_distance = 0;
      std::vector<unsigned char> leftover_right_alignment;
      int ret_code_right = AlignmentFunctionSHW(read->get_data() + query_end + 1, (clip_count_back),
                                       (int8_t *) (ref_data + alignment_position_end + 1), (clip_count_back*2),
                                       -1, parameters.match_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                       &leftover_right_start, &leftover_right_end,
                                       &leftover_right_edit_distance, leftover_right_alignment);
      /// Check if the return code is ok. Otherwise, just clip the back.
      /// Added on 15.11.2015.: check if the edit distance of the back part is too high. EDlib will automatically return an error code, but SeqAn won't.
      /// An example is when the entire back part does not match (e.g. alignment of a read to a part of the reference consisted of N bases).
      if (ret_code_right != 0 || leftover_right_edit_distance > clip_count_back/2) {
        // TODO: This is a nasty hack. EDlib used to crash when query and target are extremely small, e.g. query = "C" and target = "TC".
        // In this manner we just ignore the trailing part, and clip it.
        std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
        alignment.insert(alignment.end(), insertions_back.begin(), insertions_back.end());

      } else {
        if (leftover_right_alignment.size() == 0) {
          std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
          alignment.insert(alignment.end(), insertions_back.begin(), insertions_back.end());
        } else {

          if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
            std::string alignment_as_string = "";
            alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + query_end + 1, clip_count_back,
                                                         (const unsigned char *) (ref_data + alignment_position_end + 1), clip_count_back*2,
                                                         (unsigned char *) &(leftover_right_alignment[0]), leftover_right_alignment.size(),
                                                         (0), MYERS_MODE_SHW);
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                      FormatString("Aligning the end of the read:\n%s\n", alignment_as_string.c_str()), "[]");
          }

          /// Check for a special case when previous global alignment ended with deletions or insertions, and the new one starts with deletions or insertions.
          /// Switching from deletions to insertions is basically a mismatch streak.
          if (alignment.size() > 0 && leftover_right_alignment.size() > 0 && ((alignment.back() == EDLIB_D && leftover_right_alignment[0] == EDLIB_I) || (alignment.back() == EDLIB_I && leftover_right_alignment[0] == EDLIB_D))) {
            int64_t num_trailing_indels = 0;
            int64_t num_leading_indels = 0;
            int64_t current_op1 = alignment.size() - 1;
            while (current_op1 >= 0) {
              if ((current_op1 + 1) < alignment.size() && alignment[current_op1] != alignment[current_op1+1])
                break;
              num_trailing_indels += 1;
              current_op1 -= 1;
            }
            int64_t current_op2 = 0;
            while (current_op2 < leftover_right_alignment.size()) {
              if (current_op2 > 0 && leftover_right_alignment[current_op2] != leftover_right_alignment[current_op2-1])
                break;
              num_leading_indels += 1;
              current_op2 += 1;
            }

            int64_t min_count = std::min(num_trailing_indels, num_leading_indels);
            for (current_op1 = 0; current_op1 < min_count; current_op1++) {
              if ((ref_data + alignment_position_end + 1 + leftover_right_start - current_op1 - 1) == (read->get_data() + (query_end + 1) - current_op1))
                alignment[alignment.size() - current_op1 - 1] = EDLIB_EQUAL;
              else
                alignment[alignment.size() - current_op1 - 1] = EDLIB_X;
            }

            alignment.insert(alignment.end(), leftover_right_alignment.begin() + min_count, leftover_right_alignment.end());
          } else {
            alignment.insert(alignment.end(), leftover_right_alignment.begin(), leftover_right_alignment.end());
          }
          alignment_position_end += leftover_right_end + 1;
        }

      }
    }
  }



  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                        FormatString("alignment_position_start = %ld\nalignment_position_end = %ld\n",
                                        alignment_position_start, alignment_position_end), "LocalRealignmentLinear");
  }

  ConvertInsertionsToClipping((unsigned char *) &(alignment[0]), alignment.size());

  CountAlignmentOperations(alignment, read->get_data(), ref_data, reference_id, alignment_position_start, orientation,
                           parameters.evalue_match, parameters.evalue_mismatch, parameters.evalue_gap_open, parameters.evalue_gap_extend,
                           ret_eq_op, ret_x_op, ret_i_op, ret_d_op, ret_AS_left_part, ret_nonclipped_left_part);

  if (parameters.verbose_level > 5 && read->get_sequence_id() == parameters.debug_read) {
    std::string alignment_as_string = "";
    alignment_as_string = PrintAlignmentToString((const unsigned char *) (read->get_data()), read->get_sequence_length(),
                                               (const unsigned char *) (ref_data + alignment_position_start), (alignment_position_end - alignment_position_start + 1),
                                               (unsigned char *) &(alignment[0]), alignment.size(),
                                               (0), MYERS_MODE_NW);
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                             FormatString("Alignment:\n%s\n\nalignment_position_start = %ld\n\n", alignment_as_string.c_str(), alignment_position_start), "AnchoredAlignment");
  }



  int64_t best_aligning_position = 0;

  if (is_linear == true) {
    *ret_cigar_left_part = AlignmentToCigar((unsigned char *) &(alignment[0]), alignment.size());
    *ret_cigar_right_part = "";

    if (orientation == kForward) {
      index->RawPositionConverterWithRefId(alignment_position_start, absolute_reference_id, 0, NULL, &best_aligning_position, NULL);
    } else {
      index->RawPositionConverterWithRefId(alignment_position_end, absolute_reference_id, 0, NULL, &best_aligning_position, NULL);
      reference_id -= index->get_num_sequences_forward();
    }




    *ret_alignment_position_left_part = best_aligning_position;
    *ret_alignment_position_right_part = 0;
    *ret_orientation = orientation;
    *ret_reference_id = reference_id;
    *ret_position_ambiguity = 0;

    if (parameters.verbose_level > 5 && read->get_sequence_id() == parameters.debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Alignment array:\n"), "[]");
      for (int i1=0; i1<alignment.size(); i1++) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("%d", alignment[i1]), "[]");
      }
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("\n"), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("CIGAR string:\n%s\n", ret_cigar_left_part->c_str()), "AnchoredAlignment");

    }

    if (CheckAlignmentSane(alignment, read, index, reference_id, best_aligning_position)) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("Alignment is insane!\n"), "LocalRealignmentLinear");
      return ALIGNMENT_NOT_SANE;
    }




  } else {
    *ret_AS_right_part = *ret_AS_left_part;
    *ret_nonclipped_right_part = *ret_nonclipped_left_part;



    int64_t best_aligning_position = 0;

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("best_aligning_position_start = %ld\n", alignment_position_start), "LocalRealignmentCircular");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("best_aligning_position_end = %ld\n", alignment_position_end), "LocalRealignmentCircular");

    int64_t left_alignment_length = 0, left_alignment_start = 0, left_alignment_end = 0;
    int64_t right_alignment_length = 0, right_alignment_start = 0, right_alignment_end = 0;
    unsigned char *left_alignment = NULL;
    unsigned char *right_alignment = NULL;
    std::vector<unsigned char> alignment_right_part;
    if (ClipCircularAlignment(alignment_position_start, alignment_position_end, (unsigned char *) &(alignment[0]), alignment.size(),
                          (int64_t) (read->get_sequence_length()), (int64_t) (index->get_reference_starting_pos()[absolute_reference_id]),
                          (int64_t) (index->get_reference_lengths()[absolute_reference_id]),
                          start_offset, position_of_ref_end,
                          &left_alignment, &left_alignment_length, &left_alignment_start, &left_alignment_end,
                          &right_alignment, &right_alignment_length, &right_alignment_start, &right_alignment_end) != 0) {
      alignment.clear();
      alignment.assign(left_alignment, (left_alignment + left_alignment_length));
      if (left_alignment)
        free(left_alignment);

      alignment_right_part.clear();
      alignment_right_part.assign(right_alignment, (right_alignment + right_alignment_length));
      if (right_alignment)
        free(right_alignment);
    }

    int64_t best_aligning_position_left_part = 0;
    if (alignment.size() > 0) {
      *ret_cigar_left_part = AlignmentToCigar((unsigned char *) &(alignment[0]), alignment.size());
  //    *ret_AS_left_part = RescoreAlignment((unsigned char *) &(alignment_left_part[0]), alignment_left_part.size(), parameters.match_score, parameters.mismatch_penalty, parameters.gap_open_penalty, parameters.gap_extend_penalty);

      if (orientation == kForward) {
        index->RawPositionConverterWithRefId(left_alignment_start, absolute_reference_id, 0, NULL, &best_aligning_position_left_part, NULL);
      } else {
        index->RawPositionConverterWithRefId(left_alignment_end, absolute_reference_id, 0, NULL, &best_aligning_position_left_part, NULL);
      }
    } else {
      *ret_cigar_left_part = "";
    }

    int64_t best_aligning_position_right_part = 0;
    if (alignment_right_part.size() > 0) {
      *ret_cigar_right_part = AlignmentToCigar((unsigned char *) &(alignment_right_part[0]), alignment_right_part.size());
  //    *ret_AS_right_part = RescoreAlignment((unsigned char *) &(alignment_right_part[0]), alignment_right_part.size(), parameters.match_score, parameters.mismatch_penalty, parameters.gap_open_penalty, parameters.gap_extend_penalty);

      if (orientation == kForward) {
        index->RawPositionConverterWithRefId(right_alignment_start, absolute_reference_id, 0, NULL, &best_aligning_position_right_part, NULL);
      } else {
        index->RawPositionConverterWithRefId(right_alignment_end, absolute_reference_id, 0, NULL, &best_aligning_position_right_part, NULL);
      }
    } else {
      *ret_cigar_right_part = "";
    }

    *ret_alignment_position_left_part = best_aligning_position_left_part;
    *ret_alignment_position_right_part = best_aligning_position_right_part;
    *ret_orientation = orientation;
    *ret_reference_id = reference_id;
    *ret_position_ambiguity = 0;

    if (ref_data)
      delete[] ref_data;

    if (CheckAlignmentSane(alignment, read, index, reference_id, best_aligning_position)) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("Alignment is insane!\n"), "LocalRealignmentLinear");
      return ALIGNMENT_NOT_SANE;
    }
  }

  alignment.clear();

  return ((int) edit_distance);
}



//int GraphMap::FilterAnchorBreakpoints(const std::vector<int> &lcskpp_indices, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length, int64_t min_cluster_length, float min_cluster_coverage, OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters, std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids) {
////  int64_t min_cluster_length = 0;
////  int64_t min_covered_bases = std::max(30.0f, read->get_sequence_length() * 0.02f);
//
//  std::vector<ClusterAndIndices *> clusters;
//  ClusterAndIndices *new_cluster = NULL;
//  int64_t last_nonskipped_i = lcskpp_indices.size() + 1;
//  for (int64_t i=(lcskpp_indices.size() - 1); i >= 0; i--) {
//    /// Skip anchors which might be too erroneous.
//    int64_t current_lcskp_index = lcskpp_indices.at(i) + ref_hits_start;
//    if (CheckDistanceTooBig(owler_data, current_lcskp_index, current_lcskp_index, parameters->error_rate / 2.0f) == true)
//      continue;
//
//    if (last_nonskipped_i > lcskpp_indices.size()) {
//
//    } else {
//      /// This is going to work, because last_nonskipped_i will be set the second iteration of the loop. The value of i starts counting from int64_t i=(lcskpp_indices.size() - 1).
//      int64_t previous_lcskp_index = lcskpp_indices.at(last_nonskipped_i) + ref_hits_start;
//
//      bool wrong_to_previous1 = CheckDistanceTooBig(owler_data, previous_lcskp_index, current_lcskp_index, parameters->error_rate / 2.0f);
//      bool wrong_to_previous2 = (new_cluster->lcskpp_indices.size() < 2) ? false :
//                                (CheckDistanceTooBig(owler_data, new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-2] + ref_hits_start, current_lcskp_index, parameters->error_rate / 2.0f));
//      bool wrong_by_distance = CheckDistanceStep(owler_data, new_cluster->lcskpp_indices.front(), previous_lcskp_index, current_lcskp_index, 1.5f);
//
//      if ((wrong_to_previous1 == true && wrong_to_previous2 == true) || wrong_by_distance == true) {
//        /// In this case, the new point is a general outlier to the previous LCSk, because it doesn't fit nesither to the previous point, nor to the point before that.
//        if (new_cluster != NULL) {
//          int64_t cov_bases_read = 0, cov_bases_ref = 0;
//          CalcCoveredBases(owler_data->seed_hits2, seed_length, new_cluster->lcskpp_indices, ref_hits_start, ref_hits_end, &cov_bases_read, &cov_bases_ref);
//          new_cluster->coverage = std::max(cov_bases_read, cov_bases_ref);
//
//          int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;
//
//          if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {
//            clusters.push_back(new_cluster);
//          } else {
//            delete new_cluster;
//          }
//          new_cluster = NULL;
//        }
//      } else if (wrong_to_previous1 == true && wrong_to_previous2 == false) {
//        /// In this case, the previous point was an outlier, because the new point fits better to the one before the previous one. Overwrite the previous entry in new_cluster.
//
//        new_cluster->query.end = owler_data->seed_hits2[current_lcskp_index].query_pos + 12 - 1;
//        new_cluster->ref.end = owler_data->seed_hits2[current_lcskp_index].ref_pos + 12 - 1;
//
//        /// This should not change, as we remove 12 bases and add 12 bases.
////          new_cluster->coverage -= local_score->get_registry_entries().covered_bases_queries[previous_lcskp_index];
////          new_cluster->coverage += local_score->get_registry_entries().covered_bases_queries[current_lcskp_index];
//        new_cluster->lcskpp_indices[new_cluster->lcskpp_indices.size()-1] = current_lcskp_index - ref_hits_start;
//
//        if (new_cluster->lcskpp_indices.size() == 1) {
//          new_cluster->query.start = owler_data->seed_hits2[current_lcskp_index].query_pos;
//          new_cluster->ref.start = owler_data->seed_hits2[current_lcskp_index].ref_pos;
//        }
//        last_nonskipped_i = i;
//
//        continue;
//      }
//    }
//
//    if (new_cluster == NULL) {
//      new_cluster = new ClusterAndIndices;
//      new_cluster->query.start = owler_data->seed_hits2[current_lcskp_index].query_pos;
//      new_cluster->ref.start = owler_data->seed_hits2[current_lcskp_index].ref_pos;
//    }
//
//    new_cluster->query.end = owler_data->seed_hits2[current_lcskp_index].query_pos + 12 - 1;
//    new_cluster->ref.end = owler_data->seed_hits2[current_lcskp_index].ref_pos + 12 - 1;
//    new_cluster->num_anchors += 1;
////    new_cluster->coverage += 12;
//    new_cluster->lcskpp_indices.push_back(current_lcskp_index - ref_hits_start);
//
//    last_nonskipped_i = i;
//  }
//  if (new_cluster != NULL) {
//    int64_t cov_bases_read = 0, cov_bases_ref = 0;
//    CalcCoveredBases(owler_data->seed_hits2, seed_length, new_cluster->lcskpp_indices, ref_hits_start, ref_hits_end, &cov_bases_read, &cov_bases_ref);
//    new_cluster->coverage = std::max(cov_bases_read, cov_bases_ref);
//
//    int64_t min_covered_bases = (new_cluster->query.end - new_cluster->query.start + 1) * min_cluster_coverage;
//
//    if ((new_cluster->query.end - new_cluster->query.start + 1) >= min_cluster_length && new_cluster->coverage >= min_covered_bases) {
//      clusters.push_back(new_cluster);
//    } else {
//      delete new_cluster;
//    }
//    new_cluster = NULL;
//  }
//
//  if (ret_cluster_ids) {
//    ret_cluster_ids->clear();
//  }
//
//  int ret_val = 0;
//  /// Check if the leftover clusters are linear and that only outlier anchors are filtered. This is important for overlapping.
//  for (int64_t i=1; i<clusters.size(); i++) {
//    int64_t current_lcskp_index = clusters[i]->lcskpp_indices.front() + ref_hits_start;
//    int64_t previous_lcskp_index = clusters[i-1]->lcskpp_indices.back() + ref_hits_start;
//
//    bool wrong_to_previous1 = CheckDistanceTooBig(owler_data, previous_lcskp_index, current_lcskp_index, parameters->error_rate);
//    if (wrong_to_previous1 == true) {
//      ret_val += 1;
//    }
//  }
//
//  ret_filtered_lcskpp_indices.clear();
////  std::vector<int> cluster_indices;
//  for (int64_t i=0; i<clusters.size(); i++) {
////    int64_t cluster_length = clusters[i]->query.end - clusters[i]->query.start + 1;
////    if (cluster_length >= min_cluster_length && clusters[i]->coverage >= min_covered_bases) {
//    ret_filtered_lcskpp_indices.insert(ret_filtered_lcskpp_indices.end(), clusters[i]->lcskpp_indices.begin(), clusters[i]->lcskpp_indices.end());
//
//    /// Create indices for debugging purposes (so we can differentiate clusters).
//    if (ret_cluster_ids) {
//      std::vector<int32_t> cluster_indices(clusters[i]->lcskpp_indices.size(), i);
//      ret_cluster_ids->insert(ret_cluster_ids->end(), cluster_indices.begin(), cluster_indices.end());
//    }
////    }
//
//    if (clusters[i])
//      delete clusters[i];
//  }
//
////  int num_clusters = clusters.size();
//
//  clusters.clear();
//
////  return num_clusters;
//  return ret_val;
//}

int AnchoredAlignmentMex(bool is_linear, bool end_to_end, AlignmentFunctionTypeMex AlignmentFunctionNW, AlignmentFunctionTypeMex AlignmentFunctionSHW, const SingleSequence *read, const Index *index, const ProgramParameters &parameters, const PathGraphEntry *best_path,
                           int64_t *ret_alignment_position_left_part, std::string *ret_cigar_left_part, int64_t *ret_AS_left_part, int64_t *ret_nonclipped_left_part,
                           int64_t *ret_alignment_position_right_part, std::string *ret_cigar_right_part, int64_t *ret_AS_right_part, int64_t *ret_nonclipped_right_part,
                           SeqOrientation *ret_orientation, int64_t *ret_reference_id, int64_t *ret_position_ambiguity,
                           int64_t *ret_eq_op, int64_t *ret_x_op, int64_t *ret_i_op, int64_t *ret_d_op) {

  if (best_path->get_mapping_data().clusters.size() <= 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "No valid anchors exist!", "LocalRealignmentCircular");
    return ALIGNMENT_WRONG_CLUSTER_SIZE;
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Entering anchored alignment.\n", "AnchoredAlignment");
  int8_t *ref_data = (int8_t *) index->get_data();
  int64_t region_length_joined = 0, start_offset = 0, position_of_ref_end = 0;

  int64_t absolute_reference_id = best_path->get_region_data().reference_id;
  int64_t reference_id = best_path->get_region_data().reference_id;
  int64_t reference_start = index->get_reference_starting_pos()[absolute_reference_id];
  int64_t reference_length = index->get_reference_lengths()[absolute_reference_id];
  SeqOrientation orientation = (best_path->get_region_data().reference_id >= index->get_num_sequences_forward()) ? (kReverse) : (kForward);


  if (is_linear == false) {
    if (best_path->get_region_data().is_split == false) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Called the function for handling the circular part of the genome, but alignment is not split. best_path->region.is_split == false.\n\n"), "LocalRealignmentCircular");
      return ALIGNMENT_NOT_CIRCULAR;
    }

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, "Concatenating regions for circular alignment.\n", "AnchoredAlignment");
    ConcatenateSplitRegion(index, (Region &) best_path->get_region_data(), &ref_data, &region_length_joined, &start_offset, &position_of_ref_end);
    reference_start = 0;
    reference_length = region_length_joined;
  }




  int64_t edit_distance = 0;
  std::vector<unsigned char> alignment;

  int64_t clip_count_front = best_path->get_mapping_data().clusters.front().query.start;
  int64_t clip_count_back = read->get_sequence_length() - (best_path->get_mapping_data().clusters.back().query.end) - 1;  /////I

  int64_t alignment_position_start = best_path->get_mapping_data().clusters.front().ref.start; // - clip_count_front;
  int64_t alignment_position_end = best_path->get_mapping_data().clusters.back().ref.end; // + clip_count_back;  /////I

  int64_t query_start = best_path->get_mapping_data().clusters.front().query.start; // - clip_count_front;
  int64_t query_end = best_path->get_mapping_data().clusters.back().query.end; // + clip_count_back;  /////I

  if (clip_count_front > 0) {
    /// Check if we need to extend the alignment to the left boundary. Also, even if the user specified it, if we are to close to the boundary, just clip it.
    if (end_to_end == false || ((alignment_position_start - clip_count_front*2) < reference_start)) {
      std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
      alignment.insert(alignment.begin(), insertions_front.begin(), insertions_front.end());

    } else {
      if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                            "Aligning the begining of the read (overhang).\n", "LocalRealignmentLinear");
      }

      /// Reversing the sequences to make the semiglobal alignment of the trailing and leading parts.
      int8_t *reversed_query_front = reverse_data(read->get_data(), clip_count_front);
      int8_t *reversed_ref_front = reverse_data(ref_data + (alignment_position_start - 1) - (clip_count_front*2 - 1), clip_count_front*2);

      int64_t leftover_left_start = 0, leftover_left_end = 0, leftover_left_edit_distance = 0;
      std::vector<unsigned char> leftover_left_alignment;
      int ret_code_right = AlignmentFunctionSHW(reversed_query_front, (clip_count_front),
                                       (int8_t *) (reversed_ref_front), (clip_count_front*2),
                                       -1, parameters.match_score, parameters.mex_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                       &leftover_left_start, &leftover_left_end,
                                       &leftover_left_edit_distance, leftover_left_alignment);
      if (ret_code_right != 0) {
        // TODO: This is a nasty hack. EDlib used to crash when query and target are extremely small, e.g. query = "C" and target = "TC".
        // In this manner we just ignore the leading part, and clip it.
        std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
        alignment.insert(alignment.begin(), insertions_front.begin(), insertions_front.end());

      } else {
        if (leftover_left_alignment.size() == 0) {
          std::vector<unsigned char> insertions_front(clip_count_front, EDLIB_I);
          alignment.insert(alignment.begin(), insertions_front.begin(), insertions_front.end());
        } else {
          unsigned char *reversed_alignment = reverse_data(&(leftover_left_alignment[0]), leftover_left_alignment.size());
          alignment.insert(alignment.begin(), reversed_alignment, reversed_alignment + leftover_left_alignment.size());
          alignment_position_start -= leftover_left_end + 1;
          if (reversed_alignment)
            free(reversed_alignment);
        }

        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
          std::string alignment_as_string = "";

          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    FormatString("End of the beginning part of the read: %ld\n", (alignment_position_start - 1) - (clip_count_front*2 - 1)), "[]");
          alignment_as_string = PrintAlignmentToString((const unsigned char *) reversed_query_front, clip_count_front,
                                                       (const unsigned char *) (reversed_ref_front), clip_count_front*2,
                                                       (unsigned char *) &(leftover_left_alignment[0]), leftover_left_alignment.size(),
                                                       (0), MYERS_MODE_SHW);
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    FormatString("Aligning the beginning of the read:\n%s\n", alignment_as_string.c_str()), "[]");
        }

        if (reversed_query_front)
          free(reversed_query_front);
        if (reversed_ref_front)
          free(reversed_ref_front);
      }

    }
  }

  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
    for (int64_t i=0; i<best_path->get_mapping_data().clusters.size(); i++) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("[anchor %d] [%ld, %ld]-[%ld, %ld]\n", i, best_path->get_mapping_data().clusters[i].query.start, best_path->get_mapping_data().clusters[i].ref.start,
                                                       best_path->get_mapping_data().clusters[i].query.end, best_path->get_mapping_data().clusters[i].ref.end), "[]");
    }
  }

  for (int64_t i=0; i<best_path->get_mapping_data().clusters.size(); i++) {
    if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          "Aligning an anchor.\n", "LocalRealignmentLinear");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("[anchor %d] [%ld, %ld]-[%ld, %ld]\n", i, best_path->get_mapping_data().clusters[i].query.start, best_path->get_mapping_data().clusters[i].ref.start,
                                                       best_path->get_mapping_data().clusters[i].query.end, best_path->get_mapping_data().clusters[i].ref.end), "[]");
    }

    ///////////////////////////
    /// Align the anchor.
    int64_t query_start = best_path->get_mapping_data().clusters[i].query.start;
    int64_t query_end = best_path->get_mapping_data().clusters[i].query.end;
    int64_t ref_start = best_path->get_mapping_data().clusters[i].ref.start;
    int64_t ref_end = best_path->get_mapping_data().clusters[i].ref.end;
    int64_t query_alignment_length = query_end - query_start + 1;
    int64_t ref_alignment_length = ref_end - ref_start + 1;

    int64_t anchor_alignment_position_start = 0, anchor_alignment_position_end = 0, anchor_edit_distance = 0;
    std::vector<unsigned char> anchor_alignment;
    int ret_code1 = AlignmentFunctionNW(read->get_data() + query_start, (query_alignment_length),
                                        (int8_t *) (ref_data + ref_start), (ref_alignment_length),
                                     -1, parameters.match_score, parameters.mex_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                     &anchor_alignment_position_start, &anchor_alignment_position_end,
                                     &anchor_edit_distance, anchor_alignment);
    if (ret_code1 != 0 || anchor_alignment.size() == 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("Alignment returned with error! ret_code1 = %d\n", ret_code1), "LocalRealignmentLinear");
      return ret_code1;
    }

    if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
      std::string alignment_as_string = "";
      alignment_as_string = PrintAlignmentToString((const unsigned char *) (read->get_data() + query_start), query_alignment_length,
                                                   (const unsigned char *) (ref_data + ref_start), (ref_alignment_length),
                                                   (unsigned char *) &(anchor_alignment[0]), anchor_alignment.size(),
                                                   (0), MYERS_MODE_NW);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                FormatString("Aligned anchor %d:\n%s\n", i, alignment_as_string.c_str()), "[]");
    }

    edit_distance += anchor_edit_distance;
    /// Check for a special case when previous global alignment ended with deletions or insertions, and the new one starts with deletions or insertions.
    /// Switching from deletions to insertions is basically a mismatch streak.
    if (alignment.size() > 0 && anchor_alignment.size() > 0 && ((alignment.back() == EDLIB_D && anchor_alignment[0] == EDLIB_I) || (alignment.back() == EDLIB_I && anchor_alignment[0] == EDLIB_D))) {
      int64_t num_trailing_indels = 0;
      int64_t num_leading_indels = 0;
      int64_t current_op1 = alignment.size() - 1;
      while (current_op1 >= 0) {
        if ((current_op1 + 1) < alignment.size() && alignment[current_op1] != alignment[current_op1+1])
          break;
        num_trailing_indels += 1;
        current_op1 -= 1;
      }
      int64_t current_op2 = 0;
      while (current_op2 < anchor_alignment.size()) {
        if (current_op2 > 0 && anchor_alignment[current_op2] != anchor_alignment[current_op2-1])
          break;
        num_leading_indels += 1;
        current_op2 += 1;
      }

      int64_t min_count = std::min(num_trailing_indels, num_leading_indels);

      for (current_op1 = 0; current_op1 < min_count; current_op1++) {
        if ((ref_data + ref_start + anchor_alignment_position_start - current_op1 - 1) == (read->get_data() + (query_start) - current_op1))
          alignment[alignment.size() - current_op1 - 1] = EDLIB_EQUAL;
        else
          alignment[alignment.size() - current_op1 - 1] = EDLIB_X;
      }

      alignment.insert(alignment.end(), anchor_alignment.begin() + min_count, anchor_alignment.end());

    } else {
      alignment.insert(alignment.end(), anchor_alignment.begin(), anchor_alignment.end());
    }


    ///////////////////////////
    ///////////////////////////
    /// Align in between the anchors.
    if ((i + 1) < best_path->get_mapping_data().clusters.size()) {
      int64_t next_query_start = best_path->get_mapping_data().clusters[i+1].query.start;
      int64_t next_query_end = best_path->get_mapping_data().clusters[i+1].query.end;
      int64_t next_ref_start = best_path->get_mapping_data().clusters[i+1].ref.start;
      int64_t next_ref_end = best_path->get_mapping_data().clusters[i+1].ref.end;
      int64_t inbetween_query_length = (next_query_start - (query_end)) - 1;  /////I
      int64_t inbetween_ref_length = (next_ref_start - (ref_end)) - 1;  /////I

      /// Check if there is actually any distance between the queries, or between the references.
      /// If there is no difference, that means there is a clean insertion/deletion.
      if (inbetween_query_length <= 0 && inbetween_ref_length > 0) {
        std::vector<unsigned char> deletions_inbetween(inbetween_ref_length, EDLIB_D);
        alignment.insert(alignment.end(), deletions_inbetween.begin(), deletions_inbetween.end());

      } else if (inbetween_ref_length <= 0 && inbetween_query_length > 0) {
        std::vector<unsigned char> insertions_inbetween(inbetween_query_length, EDLIB_I);
        alignment.insert(alignment.end(), insertions_inbetween.begin(), insertions_inbetween.end());

      } else if (inbetween_query_length < 0 && inbetween_ref_length < 0) {
        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    "Problem aligning in between anchors!\n", "LocalRealignmentLinear");
        }
        return ALIGNMENT_DISTANCE_BETWEEN_ANCHORS_PROBLEM;

      } else if (inbetween_query_length > 0 && inbetween_ref_length > 0) {
        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                              "Aligning in between anchors.\n", "LocalRealignmentLinear");
        }

        int64_t between_alignment_position_start = 0, between_alignment_position_end = 0, between_anchor_edit_distance = 0;
        std::vector<unsigned char> between_anchor_alignment;
        int ret_code2 = AlignmentFunctionNW(read->get_data() + (query_end) + 1, inbetween_query_length,
                                            (int8_t *) (ref_data + ref_end) + 1, inbetween_ref_length,
                                         -1, parameters.match_score, parameters.mex_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                         &between_alignment_position_start, &between_alignment_position_end,
                                         &between_anchor_edit_distance, between_anchor_alignment);

//        int64_t between_alignment_position_start = 0, between_alignment_position_end = 0, between_anchor_edit_distance = 0;
//        std::vector<unsigned char> between_anchor_alignment;
//        int ret_code2 = AlignmentFunctionNW(read->get_data() + (query_end) + 1 - parameters.k_graph, inbetween_query_length + parameters.k_graph*2,
//                                            (int8_t *) (ref_data + ref_end) + 1 - parameters.k_graph, inbetween_ref_length + parameters.k_graph*2,
//                                         -1, parameters.match_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
//                                         &between_alignment_position_start, &between_alignment_position_end,
//                                         &between_anchor_edit_distance, between_anchor_alignment);
//
//        if (ret_code2 != 0 || between_anchor_alignment.size() == 0) {
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
//                                              FormatString("Alignment returned with error! ret_code2 = %d\n", ret_code2), "LocalRealignmentLinear");
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
//                                              FormatString("inbetween_query_length = %ld\ninbetween_ref_length = %ld\nnext_ref_start = %ld\nref_end = %ld\n",
//                                                           inbetween_query_length, inbetween_ref_length, next_ref_start, ref_end), "[]");
//          return ret_code2;
//        }
//        edit_distance += between_anchor_edit_distance;
//
//        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
//          std::string alignment_as_string = "";
//          alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + (query_end) + 1 - parameters.k_graph, inbetween_query_length + parameters.k_graph*2,
//                                                       (const unsigned char *) (ref_data + ref_end) + 1 - parameters.k_graph, inbetween_ref_length + parameters.k_graph*2,
//                                                       (unsigned char *) &(between_anchor_alignment[0]), between_anchor_alignment.size(),
//                                                       (0), MYERS_MODE_NW);
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
//                                                    FormatString("(temp) Aligning in between anchors %d and %d:\n%s\n", i, (i+1), alignment_as_string.c_str()), "[]");
//        }
//
//        between_anchor_alignment = std::vector<unsigned char>(between_anchor_alignment.begin() + parameters.k_graph, between_anchor_alignment.end() - parameters.k_graph);
//
//        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
//          std::string alignment_as_string = "";
//          alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + (query_end) + 1, inbetween_query_length,
//                                                       (const unsigned char *) (ref_data + ref_end) + 1, inbetween_ref_length,
//                                                       (unsigned char *) &(between_anchor_alignment[0]), between_anchor_alignment.size(),
//                                                       (0), MYERS_MODE_NW);
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
//                                                    FormatString("Aligning in between anchors %d and %d:\n%s\n", i, (i+1), alignment_as_string.c_str()), "[]");
//        }

        if (ret_code2 != 0 || between_anchor_alignment.size() == 0) {
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                              FormatString("Alignment returned with error! ret_code2 = %d\n", ret_code2), "LocalRealignmentLinear");
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                              FormatString("inbetween_query_length = %ld\ninbetween_ref_length = %ld\nnext_ref_start = %ld\nref_end = %ld\n",
                                                           inbetween_query_length, inbetween_ref_length, next_ref_start, ref_end), "[]");
          return ret_code2;
        }
        edit_distance += between_anchor_edit_distance;

        if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
          std::string alignment_as_string = "";
          alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + (query_end) + 1, inbetween_query_length,
                                                       (const unsigned char *) (ref_data + ref_end) + 1, inbetween_ref_length,
                                                       (unsigned char *) &(between_anchor_alignment[0]), between_anchor_alignment.size(),
                                                       (0), MYERS_MODE_NW);
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                    FormatString("Aligning in between anchors %d and %d:\n%s\n", i, (i+1), alignment_as_string.c_str()), "[]");
        }

        /// Check for a special case when previous global alignment ended with deletions or insertions, and the new one starts with deletions or insertions.
        /// Switching from deletions to insertions is basically a mismatch streak.
        if (alignment.size() > 0 && between_anchor_alignment.size() > 0 && ((alignment.back() == EDLIB_D && between_anchor_alignment[0] == EDLIB_I) || (alignment.back() == EDLIB_I && between_anchor_alignment[0] == EDLIB_D))) {
          int64_t num_trailing_indels = 0;
          int64_t num_leading_indels = 0;
          int64_t current_op1 = alignment.size() - 1;
          while (current_op1 >= 0) {
            if ((current_op1 + 1) < alignment.size() && alignment[current_op1] != alignment[current_op1+1])
              break;
            num_trailing_indels += 1;
            current_op1 -= 1;
          }
          int64_t current_op2 = 0;
          while (current_op2 < between_anchor_alignment.size()) {
            if (current_op2 > 0 && between_anchor_alignment[current_op2] != between_anchor_alignment[current_op2-1])
              break;
            num_leading_indels += 1;
            current_op2 += 1;
          }

          int64_t min_count = std::min(num_trailing_indels, num_leading_indels);
          for (current_op1 = 0; current_op1 < min_count; current_op1++) {
            if ((ref_data + ref_end + between_alignment_position_start - current_op1 - 1) == (read->get_data() + (query_end) - current_op1))
              alignment[alignment.size() - current_op1 - 1] = EDLIB_EQUAL;
            else
              alignment[alignment.size() - current_op1 - 1] = EDLIB_X;
          }

          alignment.insert(alignment.end(), between_anchor_alignment.begin() + min_count, between_anchor_alignment.end());

        } else {
          alignment.insert(alignment.end(), between_anchor_alignment.begin(), between_anchor_alignment.end());
        }
      }
    }
  }



  if (clip_count_back > 0) {
    /// Handle the clipping at the end, or extend alignment to the end of the sequence.
    if (end_to_end == false || (alignment_position_end + 1 + clip_count_back * 2) >= (reference_start + reference_length)) {
        std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
        alignment.insert(alignment.end(), insertions_back.begin(), insertions_back.end());

    } else {
      if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                            "Aligning the end of the read (overhang).\n", "LocalRealignmentLinear");
      }

      int64_t leftover_right_start = 0, leftover_right_end = 0, leftover_right_edit_distance = 0;
      std::vector<unsigned char> leftover_right_alignment;
      int ret_code_right = AlignmentFunctionSHW(read->get_data() + query_end + 1, (clip_count_back),
                                       (int8_t *) (ref_data + alignment_position_end + 1), (clip_count_back*2),
                                       -1, parameters.match_score, parameters.mex_score, -parameters.mismatch_penalty, -parameters.gap_open_penalty, -parameters.gap_extend_penalty,
                                       &leftover_right_start, &leftover_right_end,
                                       &leftover_right_edit_distance, leftover_right_alignment);
      if (ret_code_right != 0) {
        // TODO: This is a nasty hack. EDlib used to crash when query and target are extremely small, e.g. query = "C" and target = "TC".
        // In this manner we just ignore the trailing part, and clip it.
        std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
        alignment.insert(alignment.end(), insertions_back.begin(), insertions_back.end());

      } else {
        if (leftover_right_alignment.size() == 0) {
          std::vector<unsigned char> insertions_back(clip_count_back, EDLIB_I);
          alignment.insert(alignment.end(), insertions_back.begin(), insertions_back.end());
        } else {

          if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
            std::string alignment_as_string = "";
            alignment_as_string = PrintAlignmentToString((const unsigned char *) read->get_data() + query_end + 1, clip_count_back,
                                                         (const unsigned char *) (ref_data + alignment_position_end + 1), clip_count_back*2,
                                                         (unsigned char *) &(leftover_right_alignment[0]), leftover_right_alignment.size(),
                                                         (0), MYERS_MODE_SHW);
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                                      FormatString("Aligning the end of the read:\n%s\n", alignment_as_string.c_str()), "[]");
          }

          /// Check for a special case when previous global alignment ended with deletions or insertions, and the new one starts with deletions or insertions.
          /// Switching from deletions to insertions is basically a mismatch streak.
          if (alignment.size() > 0 && leftover_right_alignment.size() > 0 && ((alignment.back() == EDLIB_D && leftover_right_alignment[0] == EDLIB_I) || (alignment.back() == EDLIB_I && leftover_right_alignment[0] == EDLIB_D))) {
            int64_t num_trailing_indels = 0;
            int64_t num_leading_indels = 0;
            int64_t current_op1 = alignment.size() - 1;
            while (current_op1 >= 0) {
              if ((current_op1 + 1) < alignment.size() && alignment[current_op1] != alignment[current_op1+1])
                break;
              num_trailing_indels += 1;
              current_op1 -= 1;
            }
            int64_t current_op2 = 0;
            while (current_op2 < leftover_right_alignment.size()) {
              if (current_op2 > 0 && leftover_right_alignment[current_op2] != leftover_right_alignment[current_op2-1])
                break;
              num_leading_indels += 1;
              current_op2 += 1;
            }

            int64_t min_count = std::min(num_trailing_indels, num_leading_indels);
            for (current_op1 = 0; current_op1 < min_count; current_op1++) {
              if ((ref_data + alignment_position_end + 1 + leftover_right_start - current_op1 - 1) == (read->get_data() + (query_end + 1) - current_op1))
                alignment[alignment.size() - current_op1 - 1] = EDLIB_EQUAL;
              else
                alignment[alignment.size() - current_op1 - 1] = EDLIB_X;
            }

            alignment.insert(alignment.end(), leftover_right_alignment.begin() + min_count, leftover_right_alignment.end());
          } else {
            alignment.insert(alignment.end(), leftover_right_alignment.begin(), leftover_right_alignment.end());
          }
          alignment_position_end += leftover_right_end + 1;
        }

      }
    }
  }



  if (parameters.verbose_level > 5 && ((int64_t) read->get_sequence_id()) == parameters.debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                        FormatString("alignment_position_start = %ld\nalignment_position_end = %ld\n",
                                        alignment_position_start, alignment_position_end), "LocalRealignmentLinear");
  }

  ConvertInsertionsToClipping((unsigned char *) &(alignment[0]), alignment.size());

  CountAlignmentOperations(alignment, read->get_data(), ref_data, reference_id, alignment_position_start, orientation,
                           parameters.evalue_match, parameters.evalue_mismatch, parameters.evalue_gap_open, parameters.evalue_gap_extend,
                           ret_eq_op, ret_x_op, ret_i_op, ret_d_op, ret_AS_left_part, ret_nonclipped_left_part);

#ifndef RELEASE_VERSION
  if (parameters.verbose_level > 5 && read->get_sequence_id() == parameters.debug_read) {
    std::string alignment_as_string = "";
    alignment_as_string = PrintAlignmentToString((const unsigned char *) (read->get_data()), read->get_sequence_length(),
                                               (const unsigned char *) (ref_data + alignment_position_start), (alignment_position_end - alignment_position_start + 1),
                                               (unsigned char *) &(alignment[0]), alignment.size(),
                                               (0), MYERS_MODE_NW);
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                             FormatString("Alignment:\n%s\n\nalignment_position_start = %ld\n\n", alignment_as_string.c_str(), alignment_position_start), "AnchoredAlignment");
  }
#endif



  int64_t best_aligning_position = 0;

  if (is_linear == true) {
    *ret_cigar_left_part = AlignmentToCigar((unsigned char *) &(alignment[0]), alignment.size());
    *ret_cigar_right_part = "";

    if (orientation == kForward) {
      index->RawPositionConverterWithRefId(alignment_position_start, absolute_reference_id, 0, NULL, &best_aligning_position, NULL);
    } else {
      index->RawPositionConverterWithRefId(alignment_position_end, absolute_reference_id, 0, NULL, &best_aligning_position, NULL);
      reference_id -= index->get_num_sequences_forward();
    }




    *ret_alignment_position_left_part = best_aligning_position;
    *ret_alignment_position_right_part = 0;
    *ret_orientation = orientation;
    *ret_reference_id = reference_id;
    *ret_position_ambiguity = 0;

    if (parameters.verbose_level > 5 && read->get_sequence_id() == parameters.debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("Alignment array:\n"), "[]");
      for (int i1=0; i1<alignment.size(); i1++) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("%d", alignment[i1]), "[]");
      }
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("\n"), "[]");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("CIGAR string:\n%s\n", ret_cigar_left_part->c_str()), "AnchoredAlignment");
    }

    if (CheckAlignmentSane(alignment, read, index, reference_id, best_aligning_position)) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("Alignment is insane!\n"), "LocalRealignmentLinear");
      return ALIGNMENT_NOT_SANE;
    }




  } else {
    *ret_AS_right_part = *ret_AS_left_part;
    *ret_nonclipped_right_part = *ret_nonclipped_left_part;



    int64_t best_aligning_position = 0;

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("best_aligning_position_start = %ld\n", alignment_position_start), "LocalRealignmentCircular");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read, FormatString("best_aligning_position_end = %ld\n", alignment_position_end), "LocalRealignmentCircular");

    int64_t left_alignment_length = 0, left_alignment_start = 0, left_alignment_end = 0;
    int64_t right_alignment_length = 0, right_alignment_start = 0, right_alignment_end = 0;
    unsigned char *left_alignment = NULL;
    unsigned char *right_alignment = NULL;
    std::vector<unsigned char> alignment_right_part;
    if (ClipCircularAlignment(alignment_position_start, alignment_position_end, (unsigned char *) &(alignment[0]), alignment.size(),
                          (int64_t) (read->get_sequence_length()), (int64_t) (index->get_reference_starting_pos()[absolute_reference_id]),
                          (int64_t) (index->get_reference_lengths()[absolute_reference_id]),
                          start_offset, position_of_ref_end,
                          &left_alignment, &left_alignment_length, &left_alignment_start, &left_alignment_end,
                          &right_alignment, &right_alignment_length, &right_alignment_start, &right_alignment_end) != 0) {
      alignment.clear();
      alignment.assign(left_alignment, (left_alignment + left_alignment_length));
      if (left_alignment)
        free(left_alignment);

      alignment_right_part.clear();
      alignment_right_part.assign(right_alignment, (right_alignment + right_alignment_length));
      if (right_alignment)
        free(right_alignment);
    }

    int64_t best_aligning_position_left_part = 0;
    if (alignment.size() > 0) {
      *ret_cigar_left_part = AlignmentToCigar((unsigned char *) &(alignment[0]), alignment.size());
  //    *ret_AS_left_part = RescoreAlignment((unsigned char *) &(alignment_left_part[0]), alignment_left_part.size(), parameters.match_score, parameters.mismatch_penalty, parameters.gap_open_penalty, parameters.gap_extend_penalty);

      if (orientation == kForward) {
        index->RawPositionConverterWithRefId(left_alignment_start, absolute_reference_id, 0, NULL, &best_aligning_position_left_part, NULL);
      } else {
        index->RawPositionConverterWithRefId(left_alignment_end, absolute_reference_id, 0, NULL, &best_aligning_position_left_part, NULL);
      }
    } else {
      *ret_cigar_left_part = "";
    }

    int64_t best_aligning_position_right_part = 0;
    if (alignment_right_part.size() > 0) {
      *ret_cigar_right_part = AlignmentToCigar((unsigned char *) &(alignment_right_part[0]), alignment_right_part.size());
  //    *ret_AS_right_part = RescoreAlignment((unsigned char *) &(alignment_right_part[0]), alignment_right_part.size(), parameters.match_score, parameters.mismatch_penalty, parameters.gap_open_penalty, parameters.gap_extend_penalty);

      if (orientation == kForward) {
        index->RawPositionConverterWithRefId(right_alignment_start, absolute_reference_id, 0, NULL, &best_aligning_position_right_part, NULL);
      } else {
        index->RawPositionConverterWithRefId(right_alignment_end, absolute_reference_id, 0, NULL, &best_aligning_position_right_part, NULL);
      }
    } else {
      *ret_cigar_right_part = "";
    }

    *ret_alignment_position_left_part = best_aligning_position_left_part;
    *ret_alignment_position_right_part = best_aligning_position_right_part;
    *ret_orientation = orientation;
    *ret_reference_id = reference_id;
    *ret_position_ambiguity = 0;

    if (ref_data)
      delete[] ref_data;

    if (CheckAlignmentSane(alignment, read, index, reference_id, best_aligning_position)) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, ((int64_t) read->get_sequence_id()) == parameters.debug_read,
                                          FormatString("Alignment is insane!\n"), "LocalRealignmentLinear");
      return ALIGNMENT_NOT_SANE;
    }
  }

  alignment.clear();

  return ((int) edit_distance);
}
