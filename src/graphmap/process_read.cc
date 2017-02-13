/*
 * process_read.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: isovic
 */

#include <ctime>
#include <limits>
#include <algorithm>
#include "graphmap/graphmap.h"

#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include "alignment/alignment.h"
#include "utility/tictoc.h"

int GraphMap::ProcessRead(MappingData *mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, std::shared_ptr<is::Transcriptome> transcriptome, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params) {
  if (indexes.size() == 0 || indexes[0] == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "No reference indexes are specified."));
  }

//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Entered function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "ProcessRead");
  LOG_DEBUG_SPEC_NEWLINE;
  LOG_DEBUG_SPEC("Entered function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024));

  // If the read length is too short, call it unmapped.
  if (read->get_sequence_length() < parameters->min_read_len) {
    std::stringstream ss;
    ss << "Unmapped_5__readlength_too_short" << "__readlength=" << read->get_sequence_length() << "__limit=" << parameters->min_read_len << "__num_region_iterations=" << mapping_data->num_region_iterations;
    mapping_data->unmapped_reason += ss.str();
    return 0;
  }

  ////////////////////////////////////
  ///// Perform Region Selection /////
  ////////////////////////////////////
  clock_t begin_clock = clock();
  int64_t bin_size = (parameters->overlapper == true) ? -1 : read->get_sequence_length() / 3;

//  RegionSelection_(bin_size, mapping_data, indexes, read, parameters);
//  RegionSelectionNoBins_(bin_size, mapping_data, indexes, read, parameters);
  RegionSelectionNoCopy_(bin_size, mapping_data, indexes, read, parameters);
//  RegionSelectionNoCopyWithMap_(bin_size, mapping_data, indexes, read, parameters);
//  RegionSelectionNoCopyWithDensehash_(bin_size, mapping_data, indexes, read, parameters);

//  int64_t bin_size = (parameters->alignment_approach == "overlapper") ? -1 : 100000;
//  RegionSelectionNoCopy_(bin_size, mapping_data, indexes, read, parameters);

  clock_t end_clock = clock();
  double elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_region_selection = elapsed_secs;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n+++++++++++++++++ Region selection elapsed time: %f sec.\n\n", mapping_data->time_region_selection), "ProcessRead");

  /// Sanity check.
  if (mapping_data->bins.size() <= 0 || (mapping_data->bins.size() > 0 && mapping_data->bins.front().bin_value == 0)) {
    std::stringstream ss;
    if (mapping_data->bins.size() <= 0) {
      ss << "Unmapped_2.1_bins.size()=" << mapping_data->bins.size() << "_which_is_lte_0" << "__readlength=" << read->get_sequence_length() << "__num_region_iterations=" << mapping_data->num_region_iterations;
    } else {
      ss << "Unmapped_2.2_best_bin_has_count_0" << "__readlength=" << read->get_sequence_length() << "__num_region_iterations=" << mapping_data->num_region_iterations;
    }
    mapping_data->unmapped_reason += ss.str();
    return 0;
  }

  begin_clock = clock();

  /////////////////////////////////////////////
  ///// Create a hash index from the read /////
  /////////////////////////////////////////////
  // Create the index for the current read. This index is used in graph construction.
  std::vector<std::string> graph_shapes = {std::string(parameters->k_graph, '1')};
  auto index_read = is::createMinimizerIndex(graph_shapes, 1.0);
  SequenceFile sf_read;
  sf_read.AddSequence((SingleSequence *) read, false);
  index_read->Create(sf_read, 0.0f, false, false, 1, 1);

  //////////////////////////////
  ///// Initialize stuff.  /////
  /////////////////////////////

  // Initialize the vertices of the graph.
  mapping_data->vertices.Resize(read->get_sequence_length());

  // Initialize the iteration counter and the value after which the counter should be reset.
  mapping_data->iteration = 0;

//  float threshold_step = 0.10f;
  float bin_value_threshold = mapping_data->bins.front().bin_value;
  bin_value_threshold = std::floor(mapping_data->bins.front().bin_value * (1.0f - parameters->bin_threshold_step));
  bin_value_threshold = std::max(bin_value_threshold, 2.0f);
  int64_t num_regions_within_threshold = CountBinsWithinThreshold_(mapping_data, bin_value_threshold);

  float min_allowed_bin_value = 0.0f;
  min_allowed_bin_value = std::min(std::floor(parameters->min_bin_percent * mapping_data->bins.front().bin_value), 2.0);

//  if (parameters->parsimonious_mode) {
//    min_allowed_bin_value = (0.50f * mapping_data->bins.front().bin_value);
//
//  } else {
//    min_allowed_bin_value = (parameters->min_bin_percent * mapping_data->bins.front().bin_value);
////    min_allowed_bin_value = 0.0f;
//  }

  int64_t max_num_regions = parameters->max_num_regions;
  if (max_num_regions < 0) { max_num_regions = mapping_data->bins.size(); }

// This was commented out on 14.04.2016.
//  if (parameters->` == true) {
//    max_num_regions = mapping_data->bins.size();
//    min_allowed_bin_value = std::min(0.10f * read->get_sequence_length(), 100.0f);
//  }

  mapping_data->num_region_iterations = 0;

  LOG_DEBUG_SPEC("min_allowed_bin_value = %f, bin_value_threshold = %f\n", min_allowed_bin_value, bin_value_threshold);

  // Just verbose.
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Top 100 scoring bins:\n"), "ProcessRead");
    for (int64_t i = 0; i < mapping_data->bins.size() && i < 100; i++) {
      if (mapping_data->bins[i].bin_value == 0) { break; }
      Region region = CalcRegionFromBin_(i, mapping_data, read, parameters);
      ScoreRegistry local_score(region, i);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read,
                                   FormatString("[i = %ld] ref_id = %ld, ref_id_fwd = %ld, location_start = %ld (%lX), location_end = %ld, is_reverse = %d, vote = %ld, region_index = %ld\n",
                                                i,mapping_data->bins[i].reference_id, (mapping_data->bins[i].reference_id % indexes[0]->get_num_sequences_forward()), region.start, region.start, region.end, (int) (region.start >= indexes[0]->get_data_length_forward()), region.region_votes, region.region_index), "ProcessRead");
    }
  }

#ifndef RELEASE_VERSION
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {

    std::string cluster_path = FormatString("temp/clusters/clusters-read-%ld.csv", read->get_sequence_id());
    FILE *fp_cluster_path = fopen(cluster_path.c_str(), "w");
    if (fp_cluster_path) {
      fclose(fp_cluster_path);
    }
  }
#endif

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n\n", "ProcessRead");

  int64_t num_regions_processed = 0;
  // Process regions one by one.
  for (int64_t i = 0; i < mapping_data->bins.size() && i < max_num_regions; i++) {
    // If the ret_check value is zero, then just continue as normal.
    int ret_check = 0;

    if (parameters->overlapper == false) {
      ret_check = CheckRegionSearchFinished_(i, min_allowed_bin_value, parameters->bin_threshold_step, &bin_value_threshold, mapping_data, read, parameters);
    } else {
      if (mapping_data->bins[i].bin_value < 1 || mapping_data->bins[i].bin_value < min_allowed_bin_value) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("CheckRegionSearchFinished because mapping_data->bins[i].bin_value < 1 || mapping_data->bins[i].bin_value < min_allowed_bin_value. mapping_data->bins[i].bin_value = %f, min_allowed_bin_value = %f\n", mapping_data->bins[i].bin_value, min_allowed_bin_value), "ProcessRead");
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Last region processed: i = %ld.\n", i), "ProcessRead");
        break;
      }
    }

    // Region search needs to stop.
    if (ret_check < 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("CheckRegionSearchFinished returned with value to break! ret_check = %d\n", ret_check), "ProcessRead");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Last region processed: i = %ld.\n", i), "ProcessRead");
      break;
      // Another iteration needs to be performed.
    } else if (ret_check > 0) {
      mapping_data->num_region_iterations += 1;
    }

    num_regions_processed += 1;

    Region region = CalcRegionFromBin_(i, mapping_data, read, parameters);
    ScoreRegistry local_score(region, i);

    bool is_reverse = (region.start >= indexes[0]->get_data_length_forward());
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[i = %ld] location_start = %ld, location_end = %ld, is_reverse = %d, vote = %ld, region_index = %ld\n", i, region.start, region.end, (int) (region.start >= indexes[0]->get_data_length_forward()), region.region_votes, region.region_index), "ProcessRead");

    // Perform the GraphMap on a single region.
    GraphMap_(&local_score, index_read, mapping_data, indexes, read, parameters);

    // Just verbose.
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Local scores (raw, before LCSk):\n"), "ProcessRead");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("%s", local_score.VerboseToString().c_str()), "ProcessRead");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Running PostProcessRegionWithLCS_. j = %ld / %ld, local_score.size() = %ld\n", i, mapping_data->bins.size(), local_score.get_registry_entries().num_vertices), "ProcessRead");
    }

    if (parameters->alignment_algorithm == "sg" || parameters->alignment_algorithm == "sggotoh") {
      int ret_value_lcs = SemiglobalPostProcessRegionWithLCS_(&local_score, mapping_data, indexes, read, parameters);
    } else {
      int ret_value_lcs = AnchoredPostProcessRegionWithLCS_(&local_score, mapping_data, indexes, read, parameters);
    }

    local_score.Clear();

    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("-----\n\n"), "ProcessRead");
    }
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Last region processed: num_regions_processed = %ld.\n", num_regions_processed), "ProcessRead");

//  if (index_read)
//    delete index_read;
//  index_read = NULL;
  mapping_data->vertices.Clear();

  end_clock = clock();
  elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_mapping = elapsed_secs;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n+++++++++++++++++ Read mapping elapsed time: %f sec.\n\n", elapsed_secs), "ProcessRead");

  begin_clock = clock();

  GenerateAlignments_(mapping_data, indexes[0], transcriptome, read, parameters, evalue_params);

  end_clock = clock();
  elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_alignment = elapsed_secs;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n+++++++++++++++++ GenerateAlignments elapsed time: %f sec.\n\n", elapsed_secs), "ProcessRead");

  // Just verbose.
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || read->get_sequence_id() == parameters->debug_read), FormatString("Exiting function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "ProcessRead");
  }

  return 0;
}

void GraphMap::VerboseRegions_(const ProgramParameters* parameters, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const std::vector<Region>& regions) {
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Top 100 scoring bins:\n"), "ProcessRead");
    for (int64_t i = 0; i < regions.size() && i < 100; i++) {
      auto& region = regions[i];
      ScoreRegistry local_score(region, i);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read,
                 FormatString("[i = %ld / %ld] ref_id = %ld, ref_id_fwd = %ld, location_start = %ld, location_end = %ld, dist = %ld, is_reverse = %d, vote = %ld, region_index = %ld\n",
                 i, regions.size(), region.reference_id, (region.reference_id % indexes[0]->get_num_sequences_forward()), region.start, region.end, (region.end - region.start - 1),
                 (int) (region.start >= indexes[0]->get_data_length_forward()), region.region_votes, region.region_index), "ProcessRead");
    }
  }
}

void GraphMap::OpenDebugClustersFile_(const ProgramParameters* parameters, const SingleSequence* read) {
  #ifndef RELEASE_VERSION
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {

      std::string cluster_path = FormatString("temp/clusters/clusters-read-%ld.csv", read->get_sequence_id());
      FILE *fp_cluster_path = fopen(cluster_path.c_str(), "w");
      if (fp_cluster_path) {
        fclose(fp_cluster_path);
      }
    }
  #endif
}

int GraphMap::ProcessRead2(MappingData *mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, std::shared_ptr<is::Transcriptome> transcriptome, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params) {
  if (indexes.size() == 0 || indexes[0] == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "No reference indexes are specified."));
  }

  LOG_DEBUG_SPEC_NEWLINE;
  LOG_DEBUG_SPEC("Entered function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024));

  // If the read length is too short, call it unmapped.
  if (read->get_sequence_length() < parameters->min_read_len) {
    std::stringstream ss;
    ss << "Unmapped_5__readlength_too_short" << "__readlength=" << read->get_sequence_length() << "__limit=" << parameters->min_read_len << "__num_region_iterations=" << mapping_data->num_region_iterations;
    mapping_data->unmapped_reason += ss.str();
    return 2;
  }

  ////////////////////////////////////
  ///// Perform Region Selection /////
  ////////////////////////////////////
  std::vector<Region> regions;

  TicToc tt_region;
  tt_region.start();

  RegionSelectionWithSort_(0, mapping_data, indexes, read, parameters, regions);

  tt_region.stop();

  mapping_data->time_region_selection = tt_region.get_secs();
  LOG_DEBUG_SPEC("\n+++++++++++++++++ Region selection (sort version) elapsed time: %f sec.\n\n", mapping_data->time_region_selection);

  /// Sanity check.
  if (regions.size() == 0) {
    mapping_data->unmapped_reason += "Unmapped_2.3_no_regions_were_selected._";
    return 1;
  }



  /////////////////////////////////////////////
  ///// Start the mapping process         /////
  /////////////////////////////////////////////
  TicToc tt_mapping;
  tt_mapping.start();

  /////////////////////////////////////////////
  ///// Create a hash index from the read /////
  /////////////////////////////////////////////
  // Temporary SequenceFile for the current read.
  SequenceFile sf_read;
  sf_read.AddSequence((SingleSequence *) read, false);

  // Create an index of the current read. This index is used in graph construction.
  std::vector<std::string> graph_shapes = {std::string(parameters->k_graph, '1')};
  auto index_read = is::createMinimizerIndex(graph_shapes, 1.0);
  index_read->Create(sf_read, 0.0f, false, false, 1, 1);

  //////////////////////////////
  ///// Initialize stuff.  /////
  /////////////////////////////

  // Initialize the vertices of the graph.
  mapping_data->vertices.Resize(read->get_sequence_length());

  // Initialize the iteration counter and the value after which the counter should be reset.
  mapping_data->iteration = 0;

//  float bin_value_threshold = std::max(std::floor(regions.front().region_votes * (1.0f - parameters->bin_threshold_step)), 2.0);
  int64_t min_allowed_bin_value = 5.0f;
  LOG_DEBUG_SPEC("regions.size() = %ld\n", regions.size());

  int64_t min_num_votes = std::max((int64_t) std::floor(regions.front().region_votes * (1.0f - parameters->bin_threshold_step)), min_allowed_bin_value);
  int64_t max_num_regions = (max_num_regions < 0) ? (parameters->max_num_regions) : (regions.size());


  // Just verbose.
  LOG_DEBUG_SPEC("min_allowed_bin_value = %f, bin_value_threshold = %f\n", min_allowed_bin_value);
  LOG_DEBUG_SPEC("regions.front().region_votes = %ld, regions.back().region_votes = %ld\n", regions.front().region_votes, regions.back().region_votes);
  LOG_DEBUG_SPEC("parameters->bin_threshold_step = %f\n", parameters->bin_threshold_step);
  VerboseRegions_(parameters, indexes, read, regions);
  OpenDebugClustersFile_(parameters, read);
  LOG_DEBUG_SPEC("\n\n");

  //////////////////////////////
  ///// Perform mapping.   /////
  /////////////////////////////
  mapping_data->num_region_iterations = 0;

  int64_t num_regions_processed = 0;
  // Process regions one by one.
  for (int64_t i=0; i<regions.size() && i < max_num_regions; i++) {
    auto& region = regions[i];

    if (region.region_votes < min_num_votes) {
      break;
    }

    mapping_data->num_region_iterations += 1;
    num_regions_processed += 1;

    ScoreRegistry local_score(region, i);

    bool is_reverse = (region.start >= indexes[0]->get_data_length_forward());
    LOG_DEBUG_SPEC("[i = %ld] location_start = %ld, location_end = %ld, is_reverse = %d, vote = %ld, region_index = %ld\n",
                   i, region.start, region.end, (int) (region.start >= indexes[0]->get_data_length_forward()), region.region_votes, region.region_index);

    // Perform the GraphMap on a single region.
    GraphMap_(&local_score, index_read, mapping_data, indexes, read, parameters);

    // Just verbose.
    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LOG_DEBUG_SPEC("Local scores (raw, before LCSk):\n");
      LOG_DEBUG_SPEC("%s", local_score.VerboseToString().c_str());
      LOG_DEBUG_SPEC("Running PostProcessRegionWithLCS_. j = %ld / %ld, local_score.size() = %ld\n", i, regions.size(), local_score.get_registry_entries().num_vertices);
    }

    if (parameters->alignment_algorithm == "sg" || parameters->alignment_algorithm == "sggotoh") {
      int ret_value_lcs = SemiglobalPostProcessRegionWithLCS_(&local_score, mapping_data, indexes, read, parameters);
    } else {
      int ret_value_lcs = AnchoredPostProcessRegionWithLCS_(&local_score, mapping_data, indexes, read, parameters);
    }

    local_score.Clear();

    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LOG_DEBUG_SPEC("-----\n\n");
    }
  }

  LOG_DEBUG_SPEC("Last region processed: num_regions_processed = %ld.\n", num_regions_processed);

  mapping_data->vertices.Clear();

  tt_mapping.stop();
  mapping_data->time_mapping = tt_mapping.get_secs();
  LOG_DEBUG_SPEC("\n+++++++++++++++++ Read mapping elapsed time: %f sec.\n\n", mapping_data->time_mapping);

  TicToc tt_alignment;
  tt_alignment.start();
  GenerateAlignments_(mapping_data, indexes[0], transcriptome, read, parameters, evalue_params);
  tt_alignment.stop();

  mapping_data->time_alignment = tt_alignment.get_secs();
  LOG_DEBUG_SPEC("\n+++++++++++++++++ GenerateAlignments elapsed time: %f sec.\n\n", mapping_data->time_alignment);

  // Just verbose.
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, ((parameters->num_threads == 1) || read->get_sequence_id() == parameters->debug_read), FormatString("Exiting function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "ProcessRead");
  }

  return 0;
}



int64_t GraphMap::CountBinsWithinThreshold_(const MappingData *mapping_data, float threshold) {
  int64_t num_regions_within_threshold = 0;
  for (int64_t i = 0; i < mapping_data->bins.size(); i++) {
    if (mapping_data->bins[i].bin_value > threshold && mapping_data->bins[i].bin_value > 0)
      num_regions_within_threshold += 1;
  }
  return num_regions_within_threshold;
}

Region GraphMap::CalcRegionFromBin_(int64_t sorted_bins_index, const MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters) {
  Region ret_region;

  int64_t extension_size = read->get_sequence_length();

  int64_t bin_ref_id = mapping_data->bins[sorted_bins_index].reference_id;
  int64_t bin_id = mapping_data->bins[sorted_bins_index].bin_id;
  int64_t bin_value = mapping_data->bins[sorted_bins_index].bin_value;
  int64_t reference_start = indexes_[0]->get_reference_starting_pos()[bin_ref_id];
  int64_t reference_length = indexes_[0]->get_reference_lengths()[bin_ref_id];

  std::string rname = indexes_[0]->get_headers()[bin_ref_id % indexes_[0]->get_num_sequences_forward()];

  int64_t location_start = (mapping_data->bin_size > 0) ? bin_id * mapping_data->bin_size : 0;
  int64_t location_end = (mapping_data->bin_size > 0) ? location_start + mapping_data->bin_size : reference_length;

  bool is_split = false;
  int64_t split_start = 0;
  int64_t split_end = 0;

  if (parameters->is_reference_circular == false) {
    location_start = (location_start > extension_size) ? (location_start - extension_size) : (0);
    location_end = ((location_end + extension_size) < reference_length) ? (location_end + extension_size) : (reference_length - 1);
    location_start += reference_start;
    location_end += reference_start;
    is_split = false;

  } else {

    // Check if this is the beginning of the reference sequence and there is an overhang to the left.
    if (location_start < extension_size) {
      location_start = (location_start > extension_size) ? (location_start - extension_size) : (0);
      location_end = ((location_end + extension_size) < reference_length) ? (location_end + extension_size) : (reference_length - 1);
      is_split = true;
      split_start = reference_start + reference_length - extension_size;
//      split_start = std::max((reference_start + reference_length - extension_size), ((int64_t) 0));
      split_end = reference_start + reference_length - 1;
      location_start += reference_start;
      location_end += reference_start;

      // A situation can occur when a read is longer than both fwd+rev reference. An example of such a read is in the read 9488 from
      // the Mikheyev & Tin dataset of the Lambda genome. The problem is that the read is absurdly long, around 110kbp, while the
      // reference is only 48502bp long. The alignment might make sense in the linear way, but in circular it would mean that
      // we circled the entire genome more than 2x, which makes no sense whatsoever. The resulting mpileup had 48502 extra bases
      // with coverage 1x, which equals this read. By setting the is_split = false, the read will be handled like the reference
      // is linear, instead of the circular case.
      if (split_start < 0)
        is_split = false;

      // Check if this is the end of the reference sequence and there is an overhang to the right.
    } else if ((location_end + extension_size) >= reference_length) {
      location_start = (location_start > extension_size) ? (location_start - extension_size) : (0);
      location_end = ((location_end + extension_size) < reference_length) ? (location_end + extension_size) : (reference_length - 1);
      is_split = true;
      split_start = reference_start + 0;
      split_end = reference_start + extension_size - 1;
      location_start += reference_start;
      location_end += reference_start;

      // Exact same situation as the previous if statement - if the read is longer than 2x the reference, something is wrong with the data.
      // Handle this case as linear alignment.
      if (split_end >= (reference_start + reference_length))
        is_split = false;

      // Otherwise, business as usual.
    } else {
      location_start = (location_start > extension_size) ? (location_start - extension_size) : (0);
      location_end = ((location_end + extension_size) < reference_length) ? (location_end + extension_size) : (reference_length - 1);
      location_start += reference_start;
      location_end += reference_start;
      is_split = false;

    }
  }

  ret_region.start = location_start;
  ret_region.end = location_end;
  ret_region.reference_id = bin_ref_id;
  ret_region.region_index = bin_id;
  ret_region.rname = rname;
  ret_region.region_votes = bin_value;
  ret_region.is_split = is_split;
  ret_region.split_start = split_start;
  ret_region.split_end = split_end;

  return ret_region;
}

// If EvaluateMappings_ returns 0, that means that the mappings have passed all tests, and we can stop the search.
// Otherwise, another iteration needs to be performed. Before allowing that, we need to check whether some conditions have
// been satisfied, such as that the bin_value is not too small (smaller than the threshold).
// Function returns 0 if the iteration has not finished yet, otherwise it returns value < 0 if the search needs to be stopped.
// If it returns value > 0, then the search needs to continue, but the ret code specifies why (return from EvaluateMappings_ function).
int GraphMap::CheckRegionSearchFinished_(int64_t current_region, float min_allowed_bin_value, float threshold_step, float *bin_value_threshold, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters) {
  if (mapping_data->bins[current_region].bin_value < (*bin_value_threshold)) {
    LOG_DEBUG_SPEC("mapping_data->bins[current_region].bin_value = %f, *bin_value_threshold = %f\n", mapping_data->bins[current_region].bin_value, *bin_value_threshold);
    // Check if there are "ok" mappings, should we continue iterating.
    int ret_evaluate_mappings = EvaluateMappings_(mapping_data, read, parameters);

    if (ret_evaluate_mappings == 0) {
      LOG_DEBUG_SPEC("\nret_evaluate_mappings = %d\n", ret_evaluate_mappings);
      LOG_DEBUG_SPEC("Stopping because all basic thresholds have been satisfied in the region batch. A potentially good mapping is already found.\n\n\n");

      return -1;

    } else {
      LOG_DEBUG_SPEC("\nret_evaluate_mappings = %d\n", ret_evaluate_mappings);

      if (mapping_data->bins[current_region].bin_value < 2) {
        LOG_DEBUG_SPEC("Stopping because bins_[i].bin_value < 1.\n\n\n");
        return -2;
      }

      // We need to go into another iteration.
      *bin_value_threshold -= std::floor(mapping_data->bins.front().bin_value * (threshold_step));
      *bin_value_threshold = std::max(*bin_value_threshold, 2.0f);
      int64_t num_regions_within_threshold = CountBinsWithinThreshold_(mapping_data, *bin_value_threshold);

      if (*bin_value_threshold < min_allowed_bin_value) {
        LOG_DEBUG_SPEC("Stopping because bin_value_threshold < (%f * bins_.front().bin_value). *bin_value_threshold = %f, min_allowed_bin_value = %f\n\n\n", threshold_step, *bin_value_threshold, min_allowed_bin_value);
        return -3;
      }

      if (num_regions_within_threshold > parameters->max_num_regions) {
        LOG_DEBUG_SPEC("Stopping because num_regions_within_threshold > parameters->max_num_regions. (num_regions_within_threshold = %ld, parameters->max_num_regions = %ld)\n\n\n", num_regions_within_threshold, parameters->max_num_regions);
        return -4;
      }
    }

    return ret_evaluate_mappings;
  }

  return 0;
}

int GraphMap::CheckRegionSearchFinished2_(int64_t current_region, float min_allowed_bin_value, float threshold_step,
                                          float *bin_value_threshold, MappingData *mapping_data, const std::vector<Region> &regions,
                                          const SingleSequence *read, const ProgramParameters *parameters) {
  float region_votes = regions[current_region].region_votes;

  LOG_DEBUG_SPEC("*bin_value_threshold = %f, region_votes = %f\n", *bin_value_threshold, region_votes);

  if (region_votes < (*bin_value_threshold)) {
    LOG_DEBUG_SPEC("region_votes = %f, *bin_value_threshold = %f\n", region_votes, *bin_value_threshold);
    // Check if there are "ok" mappings, should we continue iterating.
    int ret_evaluate_mappings = EvaluateMappings_(mapping_data, read, parameters);

    if (ret_evaluate_mappings == 0) {
      LOG_DEBUG_SPEC("\nret_evaluate_mappings = %d\n", ret_evaluate_mappings);
      LOG_DEBUG_SPEC("Stopping because all basic thresholds have been satisfied in the region batch. A potentially good mapping is already found.\n\n\n");

      return -1;

    } else {
      LOG_DEBUG_SPEC("\nret_evaluate_mappings = %d\n", ret_evaluate_mappings);

      if (regions[current_region].region_votes < 2) {
        LOG_DEBUG_SPEC("Stopping because region_votes < 1.\n\n\n");
        return -2;
      }

      // We need to go into another iteration.
      *bin_value_threshold -= std::floor(region_votes * (threshold_step));
      *bin_value_threshold = std::max(*bin_value_threshold, 2.0f);
//      int64_t num_regions_within_threshold = CountBinsWithinThreshold_(mapping_data, *bin_value_threshold);

      if (*bin_value_threshold < min_allowed_bin_value) {
        LOG_DEBUG_SPEC("Stopping because bin_value_threshold < (%f * bins_.front().bin_value). *bin_value_threshold = %f, min_allowed_bin_value = %f\n\n\n", threshold_step, *bin_value_threshold, min_allowed_bin_value);
        return -3;
      }

//      if (num_regions_within_threshold > parameters->max_num_regions) {
//        LOG_DEBUG_SPEC("Stopping because num_regions_within_threshold > parameters->max_num_regions. (num_regions_within_threshold = %ld, parameters->max_num_regions = %ld)\n\n\n", num_regions_within_threshold, parameters->max_num_regions);
//        return -4;
//      }
    }

    return ret_evaluate_mappings;
  }

  return 0;
}

// Returns 0 if a mapping seems to satisfy some basic thresholds. Otherwise, returns value > 0.
int GraphMap::EvaluateMappings_(MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters) {
  if (mapping_data->intermediate_mappings.size() == 0)
    return 1;

  std::sort(mapping_data->intermediate_mappings.begin(), mapping_data->intermediate_mappings.end(), path_graph_entry_greater_than_key());

  auto *best_entry = (mapping_data->intermediate_mappings.at(0));

  CollectFinalMappingsAndMapQ_(false, mapping_data, read, parameters);

  int64_t num_n_bases = ((int64_t) read->CalcNumberNBases());  // The number of 'N' bases in the read.
  int64_t num_valid_bases = ((int64_t) read->get_sequence_length()) - num_n_bases;  // The number of bases which are not 'N' in the read.
  int64_t threshold = (num_valid_bases * 0.10f);
  if (num_n_bases > 0)
    threshold = (num_valid_bases * 0.07f);

  // Threshold the alignment scores.
  if (best_entry->get_mapping_data().cov_bases_max < threshold)
    return 2;

  LOG_DEBUG_SPEC("best_entry->get_fpfilter_cov_bases() = %f < 0.50f\n", best_entry->get_fpfilter_cov_bases());
  LOG_DEBUG_SPEC("best_entry->get_fpfilter_query_len() = %f < 0.75f\n", best_entry->get_fpfilter_query_len());
  LOG_DEBUG_SPEC("best_entry->get_fpfilter_std() = %f < 0.20f\n", best_entry->get_fpfilter_std());

  if (best_entry->get_fpfilter_cov_bases() < 0.50f ||      // 0.20f
      best_entry->get_fpfilter_query_len() < 0.55f ||      // 0.20f
      best_entry->get_fpfilter_std() < 0.20f)              // 0.20f
    return 3;

  if (best_entry->CalcDistanceRatioSuppress() > parameters->error_rate)
    return 4;

  return 0;
}

int GraphMap::GenerateAlignments_(MappingData *mapping_data, std::shared_ptr<is::MinimizerIndex> index, std::shared_ptr<is::Transcriptome> transcriptome, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params) {
  if (mapping_data->intermediate_mappings.size() == 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("mapping_data->intermediate_mappings.size() == 0\n"), "GenerateAlignments_");
    return 1;
  }

  EvaluateMappings_(mapping_data, read, parameters);
  CollectFinalMappingsAndMapQ_(true, mapping_data, read, parameters);



  if (parameters->mapq_threshold >= 0 && mapping_data->mapping_quality < parameters->mapq_threshold) {
    mapping_data->unmapped_reason += FormatString("__Unmapped_because_mapq_too_low__mapq<%d", parameters->mapq_threshold);
    return 0;
  }

  for (int64_t i = 0; i < mapping_data->final_mapping_ptrs.size(); i++) {
    auto region_data = mapping_data->final_mapping_ptrs.at(i);

    /// If the region wasn't mapped for some reason, no need to perform alignment.
    if (region_data->get_mapping_data().is_mapped == false) {
      region_data->SetAligned(false);
      region_data->get_mapping_metadata().unmapped_reason += FormatString("__Unmapped_final_region_%d_because__region_data.is_mapped==false", i);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().is_mapped == false\n"), "GenerateAlignments_");
      continue;
    }

    /// Align the region and measure the time for execution.
    clock_t begin_clock = clock();
    int ret_aln = AlignRegion(read, index, transcriptome, parameters, evalue_params, true, region_data);
    clock_t end_clock = clock();
    double elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n+++++++++++++++++ Alignment elapsed time: %f sec.\n\n", elapsed_secs), "GenerateAlignments_");

    /// Set the number of different regions as an estimate of the mapping quality.
    for (int32_t j=0; j<region_data->get_alignments().size(); j++) {
      region_data->get_alignments()[j].num_secondary_alns = mapping_data->num_similar_mappings;
      region_data->get_alignments()[j].mapping_quality = mapping_data->mapping_quality;
    }

    /// Set the timing statistics.
    region_data->get_mapping_metadata().time_region_selection = mapping_data->time_region_selection;
    region_data->get_mapping_metadata().time_mapping = mapping_data->time_mapping;
    region_data->get_mapping_metadata().time_alignment = elapsed_secs;
    region_data->get_mapping_metadata().time_region_seed_lookup = mapping_data->time_region_seed_lookup;
    region_data->get_mapping_metadata().time_region_hitsort = mapping_data->time_region_hitsort;
    region_data->get_mapping_metadata().time_region_conversion = mapping_data->time_region_conversion;
    region_data->get_mapping_metadata().time_region_alloc = mapping_data->time_region_alloc;
    region_data->get_mapping_metadata().time_region_counting = mapping_data->time_region_counting;
    mapping_data->time_alignment = elapsed_secs;

    /// Check if something went wrong and the read is unmapped.
    if (ret_aln != 0) {
      region_data->get_mapping_metadata().unmapped_reason += FormatString("__AlignRegion_returned_with_error__ret_value=%d", ret_aln);
      region_data->SetAligned(false);

      /// Keep the output if alignment is insane for debug purposes.
      if (ret_aln != ALIGNMENT_NOT_SANE) {
        region_data->get_mapping_metadata().unmapped_reason += FormatString("__(Alignment_is_sane_but_there_is_another_problem)");
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Alignment sane, but ret_aln != 0 (ret_aln = %d)!\n", ret_aln), "GenerateAlignments_");
        continue;
      } else {
        region_data->get_mapping_metadata().unmapped_reason += FormatString("__(Alignment_is_insane)");
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Alignment is insane (ret_aln = %d)!\n", ret_aln), "GenerateAlignments_");
      }
    }

    /// Check if the alignment E-value satisfies the threshold.
    if (parameters->evalue_threshold >= ((double) 0.0)) {
      for (int32_t j=0; j<region_data->get_alignments().size(); j++) {
        if (region_data->get_alignments()[j].is_aligned == true && region_data->get_alignments()[j].evalue > parameters->evalue_threshold) {
          region_data->get_alignments()[j].is_aligned = false;
          region_data->get_mapping_metadata().unmapped_reason += FormatString("_evalue=%E>%E", region_data->get_alignments()[j].evalue, parameters->evalue_threshold);
          LOG_DEBUG_SPEC("Alignment %d/%d has E-value %E > %E. Marking alignment as unaligned.\n", region_data->get_alignments()[j].evalue, parameters->evalue_threshold);
        }
      }
    }
  }

  /// Check if after everything, there are no valid alignments. If so, concatenate all unmapped_reasons into the top one.
  if (mapping_data->IsAligned() == false) {
    for (int64_t i = 0; i < mapping_data->final_mapping_ptrs.size(); i++) {
      auto region_data = mapping_data->final_mapping_ptrs.at(i);
      for (int32_t j=0; j<region_data->get_alignments().size(); j++) {
        mapping_data->unmapped_reason += region_data->get_mapping_metadata().unmapped_reason;
      }
    }
  }

  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n\n\n"), "[]");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Intermediate mappings:\n%s", mapping_data->VerboseIntermediateMappingsToString(index, read).c_str()), "GenerateAlignments_");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n\n\n"), "[]");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Final mappings:\n%s", mapping_data->VerboseFinalMappingsToString(index, read).c_str()), "GenerateAlignments_");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n\n\n"), "[]");
  }

  return 0;
}

int GraphMap::CollectAlignments(const SingleSequence *read, const ProgramParameters *parameters, MappingData *mapping_data, std::string &ret_aln_lines) {
  std::stringstream ss;

  int64_t num_mapped_alignments = 0;
  int64_t num_unmapped_alignments = 0;
  for (int64_t i = 0; i < mapping_data->final_mapping_ptrs.size(); i++) {
    if (mapping_data->final_mapping_ptrs.at(i)->IsAligned() == false) {
      num_unmapped_alignments += 1;
      continue;
    }
    if (ss.tellp() > 0)
      ss << "\n";

    if (parameters->outfmt == "sam") {
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateSAM((num_mapped_alignments == 0), parameters->verbose_sam_output);  // TODO: Don't make the first alignment primary by default, but the best one.
    } else if (parameters->outfmt == "afg") {
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateAFG();
    } else if (parameters->outfmt == "m5") {
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateM5((num_mapped_alignments == 0), parameters->verbose_sam_output);
    } else if (parameters->outfmt == "mhap") {
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateMHAP((num_mapped_alignments == 0), parameters->verbose_sam_output);
    } else {  // Default to SAM output if the specified format is unknown.
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateSAM((num_mapped_alignments == 0), parameters->verbose_sam_output);
    }
    num_mapped_alignments += 1;
  }

  // If all reported alignments have been declared unmapped for some reason, output one of them to be consistent
  // and provide an alignment line for each read.
  // Also, there will always be at leas one final_mapping_ptrs_ entry, because we have handled the size() == 0 case
  // a little bit up.
  if (num_unmapped_alignments == mapping_data->final_mapping_ptrs.size() || mapping_data->unmapped_reason.size() > 0) {
    std::stringstream temp_ss;
    temp_ss << "__num_region_iterations=" << mapping_data->num_region_iterations;

    if (num_unmapped_alignments == 0)
      mapping_data->unmapped_reason += std::string("__Unmapped_3.1__no_valid_graph_paths.") + temp_ss.str();
    else
      mapping_data->unmapped_reason += std::string("__Unmapped_3.2__no_valid_graph_paths.") + temp_ss.str();
  }

  if (mapping_data->unmapped_reason.size() > 0) {
    if (parameters->outfmt == "sam") {
      ret_aln_lines = GenerateUnmappedSamLine_(read, mapping_data->unmapped_reason, parameters->verbose_sam_output);
    } else if (parameters->outfmt == "afg") {
      // In AFG format there is no need to report 'unmapped' (or non-overlapping) reads).
    } else if (parameters->outfmt == "m5") {

    } else {  // Default to SAM output if the specified format is unknown.
      ret_aln_lines = GenerateUnmappedSamLine_(read, mapping_data->unmapped_reason, parameters->verbose_sam_output);
    }
    return STATE_UNMAPPED;
  }

  ret_aln_lines = ss.str();

  return STATE_MAPPED;
}

int GraphMap::CollectFinalMappingsAndMapQ_(bool generate_final_mapping_ptrs, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters) {
  auto *first_entry = (mapping_data->intermediate_mappings.at(0));

  if (generate_final_mapping_ptrs == true) {
    mapping_data->final_mapping_ptrs.clear();
    if (mapping_data->intermediate_mappings.at(0)->get_mapping_data().is_mapped == true)
      mapping_data->final_mapping_ptrs.push_back((mapping_data->intermediate_mappings.at(0)));
  }

  int64_t num_cov_bases_query = first_entry->get_mapping_data().cov_bases_query;
  int64_t min_num_cov_bases_query = (int64_t) std::floor(((float) num_cov_bases_query) * (1.0f - parameters->margin_for_ambiguity));

  int64_t ambiguity = 1;
  int64_t num_equal_mappings = 1;
  for (int64_t i = 1; i < mapping_data->intermediate_mappings.size(); i++) {
    if (mapping_data->intermediate_mappings.at(i)->get_mapping_data().is_mapped == true) {
      // Find all the paths that are longer than 10% of the longest one.
      // if (mapping_data->intermediate_mappings.at(i).lcs_length >= min_score && mapping_data->intermediate_mappings.at(i).lcs_length <= max_score) { // >= allowed_score_difference_for_ambiguity) {
// This was used until 30.03.2016.      if ((parameters->margin_for_ambiguity == 0.0f && mapping_data->intermediate_mappings.at(i)->get_mapping_data().num_covering_kmers == num_kmers) || (parameters->margin_for_ambiguity != 0.0f && mapping_data->intermediate_mappings.at(i)->get_mapping_data().num_covering_kmers >= min_num_kmers)) {
      if ((parameters->margin_for_ambiguity == 0.0f && mapping_data->intermediate_mappings.at(i)->get_mapping_data().cov_bases_query == num_cov_bases_query) || (parameters->margin_for_ambiguity != 0.0f && mapping_data->intermediate_mappings.at(i)->get_mapping_data().cov_bases_query >= min_num_cov_bases_query)) {
        // Check if the path is different than the best one - sometimes they can be near duplicates.
        if ((mapping_data->intermediate_mappings.at(i)->get_mapping_data().ref_coords.start != first_entry->get_mapping_data().ref_coords.start || mapping_data->intermediate_mappings.at(i)->get_mapping_data().ref_coords.end != first_entry->get_mapping_data().ref_coords.end)) {
          // Count them for ambiguity (used in mapping quality).
          if (generate_final_mapping_ptrs == true) {
            if (parameters->output_multiple_alignments) {
              mapping_data->final_mapping_ptrs.push_back((mapping_data->intermediate_mappings.at(i)));
            }
          }

          ambiguity += 1;

          // Count only the paths that are not the same length of the longest path. Might be useful for metagenomic,
          // as some genomes can have identical regions.
          if (mapping_data->intermediate_mappings.at(i)->get_mapping_data().cov_bases_query == first_entry->get_mapping_data().cov_bases_query) {
            num_equal_mappings += 1;
          }
        }
      }
    }
  }

  // Calculate the mapping quality for all of the similar alignments.
  float mapping_quality_float = (ambiguity > 1) ? (-10.0f * log10((1.0f - (1.0f / ((float) ambiguity))))) : (ambiguity == 1) ? (-10.0f * log10(0.0001f)) : 0.0f;
  uint8_t mapping_quality = (uint8_t) round(std::max(mapping_quality_float, 0.0f));

  int64_t metagen_ambiguity = ambiguity - num_equal_mappings + 1;
  float metagen_alignment_score_float = (metagen_ambiguity > 1) ? (-10.0f * log10((1.0f - (1.0f / ((float) metagen_ambiguity))))) : (metagen_ambiguity == 1) ? (-10.0f * log10(0.0001f)) : 0.0f;
  int64_t metagen_alignment_score = (int64_t) round(std::max(metagen_alignment_score_float, 0.0f));

  mapping_data->num_similar_mappings = ambiguity;
  mapping_data->num_same_mappings = num_equal_mappings;
  mapping_data->mapping_quality = mapping_quality;
  mapping_data->metagen_alignment_score = metagen_alignment_score;

  return 0;
}

std::string GraphMap::GenerateUnmappedSamLine_(const SingleSequence *read, const std::string& unmapped_reason, int64_t verbose_sam_output) const {
  std::stringstream ss;

  std::string qname = ((std::string) (read->get_header()));
  std::string qname_for_out = (verbose_sam_output < 4) ? (TrimToFirstSpace(qname)) : (qname);

  uint32_t flag = SAM_THIS_SEG_UNMAPPED;

  ss << qname_for_out << "\t";
  ss << flag << "\t";
  ss << SAM_DEFAULT_RNAME << "\t";
  ss << SAM_DEFAULT_POS << "\t";
  ss << SAM_DEFAULT_MAPQ << "\t";      // To avoid confusion with the definition of mapping quality, we will use the value of 255, and report the actual quality as AS optional parameter.
  ss << SAM_DEFAULT_CIGAR << "\t";
  ss << SAM_DEFAULT_RNEXT << "\t";
  ss << SAM_DEFAULT_PNEXT << "\t";
  ss << SAM_DEFAULT_TLEN << "\t";

  if (read->get_data_length() > 0) {
    ss << read->get_data() << "\t";
  } else {
    ss << "*\t";
  }

  if (verbose_sam_output < 5 && read->get_quality_length() > 0 &&
      read->get_quality() != NULL) {
    ss << read->get_quality() << "\t";
  } else {
    ss << "*" << "\t";
  }

  std::stringstream ss_optional1;
  ss_optional1 << "NM:i:" << -1 << "\t";  // Specified by SAM format.
  ss_optional1 << "AS:i:" << -((int64_t) read->get_sequence_length()) << "\t";
  ss_optional1 << "H0:i:" << 0 << "\t";  // Specified by SAM format.
  ss_optional1 << "ZE:f:" << std::numeric_limits<float>::infinity() << "\t";
  ss_optional1 << "ZF:f:" << 0.0f << "\t";
  ss_optional1 << "ZQ:i:" << read->get_sequence_length() << "\t";
  ss_optional1 << "ZR:i:" << 0;

  ss << ss_optional1.str();

  if (verbose_sam_output >= 3) {
    ss << "\tX3:Z:" << unmapped_reason;
  }

  return ss.str();
}
