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
//#include "sam/sam_entry.h"
#include "index/index_hash.h"

#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include "alignment/alignment.h"

int GraphMap::ProcessRead(MappingData *mapping_data, const std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params) {
  if (indexes.size() == 0 || indexes[0] == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "No reference indexes are specified."));
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Entered function. [time: %.2f sec, RSS: %ld MB, peakRSS: %ld MB]\n", (((float) (clock())) / CLOCKS_PER_SEC), getCurrentRSS() / (1024 * 1024), getPeakRSS() / (1024 * 1024)), "ProcessRead");

  // If the read length is too short, call it unmapped.
  if (read->get_sequence_length() < 80) {
    std::stringstream ss;
    ss << "Unmapped_5__readlength_too_short" << "__readlength=" << read->get_sequence_length() << "__limit=" << 80 << "__num_region_iterations=" << mapping_data->num_region_iterations;
    mapping_data->unmapped_reason += ss.str();
    return 0;
  }

  ////////////////////////////////////
  ///// Perform Region Selection /////
  ////////////////////////////////////
  clock_t begin_clock = clock();
  int64_t bin_size = (parameters->alignment_approach == "overlapper") ? -1 : read->get_sequence_length() / 3;

//  RegionSelection_(bin_size, mapping_data, indexes, read, parameters);
//  RegionSelectionNoBins_(bin_size, mapping_data, indexes, read, parameters);
  RegionSelectionNoCopy_(bin_size, mapping_data, indexes, read, parameters);

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
  Index *index_read = NULL;
  if (parameters->k_graph < 10) {
    index_read = new IndexHash();
    ((IndexHash *) index_read)->set_k(parameters->k_graph);
  } else {
    index_read = new IndexSA();
  }
  index_read->GenerateFromSingleSequenceOnlyForward(*read);

  //////////////////////////////
  ///// Initialize stuff.  /////
  /////////////////////////////

  // Initialize the vertices of the graph.
  mapping_data->vertices.Resize(read->get_sequence_length());

  // Initialize the iteration counter and the value after which the counter should be reset.
  mapping_data->iteration = 0;

  float threshold_step = 0.10f;
  float bin_value_threshold = mapping_data->bins.front().bin_value;
  bin_value_threshold = (mapping_data->bins.front().bin_value * (1.0f - threshold_step));
  bin_value_threshold = std::max(bin_value_threshold, 2.0f);
  int64_t num_regions_within_threshold = CountBinsWithinThreshold_(mapping_data, bin_value_threshold);
  float min_allowed_bin_value = 0.0f;
  if (parameters->parsimonious_mode) {
    min_allowed_bin_value = (0.50f * mapping_data->bins.front().bin_value);

  } else {
    min_allowed_bin_value = (0.75f * mapping_data->bins.front().bin_value);
//    min_allowed_bin_value = 0.0f;
  }

  int64_t max_num_regions = parameters->max_num_regions;
  if (parameters->alignment_approach == "overlapper") {
    max_num_regions = mapping_data->bins.size();
    min_allowed_bin_value = std::min(0.10f * read->get_sequence_length(), 100.0f);
  }

  mapping_data->num_region_iterations = 0;

  // Just verbose.
  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Top 10 scoring bins:\n"), "ProcessRead");
    for (int64_t i = 0; i < mapping_data->bins.size() && i < 10; i++) {
      Region region = CalcRegionFromBin_(i, mapping_data, read, parameters);
      ScoreRegistry local_score(region, i);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read,
                                   FormatString("[i = %ld] ref_id = %ld, location_start = %ld (%lX), location_end = %ld, is_reverse = %d, vote = %ld, region_index = %ld\n",
                                                i, (mapping_data->bins[i].reference_id % indexes[0]->get_num_sequences_forward()), region.start, region.start, region.end, (int) (region.start >= indexes[0]->get_data_length_forward()), region.region_votes, region.region_index), "ProcessRead");
    }
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n\n", "ProcessRead");

  // Process regions one by one.
  for (int64_t i = 0; i < mapping_data->bins.size() && i < max_num_regions; i++) {
    // If the ret_check value is zero, then just continue as normal.
    int ret_check = 0;

    if (parameters->alignment_approach != "overlapper") {
      ret_check = CheckRegionSearchFinished_(i, min_allowed_bin_value, threshold_step, &bin_value_threshold, mapping_data, read, parameters);
    } else {
      if (mapping_data->bins[i].bin_value < 1 || mapping_data->bins[i].bin_value < min_allowed_bin_value)
        break;
    }

    // Region search needs to stop.
    if (ret_check < 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("CheckRegionSearchFinished returned with value to break! ret_check = %d\n", ret_check), "ProcessRead");
      break;
      // Another iteration needs to be performed.
    } else if (ret_check > 0) {
      mapping_data->num_region_iterations += 1;
    }

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

    if (parameters->alignment_algorithm == "myers" || parameters->alignment_algorithm == "gotoh") {
      int ret_value_lcs = PostProcessRegionWithLCS_(&local_score, mapping_data, indexes, read, parameters);
    } else {
      int ret_value_lcs = AnchoredPostProcessRegionWithLCS_(&local_score, mapping_data, indexes, read, parameters);
    }

    local_score.Clear();

    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("-----\n\n"), "ProcessRead");
    }

//    CheckMinimumMappingConditions(mapping_data, parameters);
  }

  if (index_read)
    delete index_read;
  index_read = NULL;
  mapping_data->vertices.Clear();

  end_clock = clock();
  elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_mapping = elapsed_secs;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n+++++++++++++++++ Read mapping elapsed time: %f sec.\n\n", elapsed_secs), "ProcessRead");

  begin_clock = clock();

//  GenerateAlignments_(mapping_data, indexes[0], read, parameters, evalue_params);
  GenerateAlignments2_(mapping_data, indexes[0], read, parameters, evalue_params);

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
//      split_end = std::min((reference_start + extension_size - 1), (reference_start + reference_length));
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
    // Check if there are "ok" mappings, should we continue iterating.
    int ret_evaluate_mappings = EvaluateMappings_(false, mapping_data, read, parameters);

    if (ret_evaluate_mappings == 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("\nret_evaluate_mappings = %d\n", ret_evaluate_mappings), "RunAlignment");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Stopping.\n\n\n"), "RunAlignment");

      return -1;

    } else {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("\nret_evaluate_mappings = %d\n", ret_evaluate_mappings), "RunAlignment");

      if (mapping_data->bins[current_region].bin_value < 2) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Stopping because bins_[i].bin_value < 1.\n\n\n"), "RunAlignment");
        return -2;
      }

      // We need to go into another iteration.
      *bin_value_threshold -= (mapping_data->bins.front().bin_value * (threshold_step));
      *bin_value_threshold = std::max(*bin_value_threshold, 2.0f);
      int64_t num_regions_within_threshold = CountBinsWithinThreshold_(mapping_data, *bin_value_threshold);

      if (*bin_value_threshold < min_allowed_bin_value) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Stopping because bin_value_threshold < (0.25 * bins_.front().bin_value).\n\n\n"), "RunAlignment");
//          std::cout << "\tbins_.front().bin_value = " << bins_.front().bin_value << "\n";
//          std::cout << "\tbin_value_threshold = " << bin_value_threshold << "\n";
//          std::cout << "\t0.25 * bins_.front().bin_value = " << 0.25 * bins_.front().bin_value << "\n";
//          std::cout << "\tthreshold_step = " << threshold_step << "\n";
//          std::cout << "\tnum_regions_within_threshold = " << num_regions_within_threshold << "\n";
//          std::cout << "\n\n\n";
//          fflush(stdout);
        return -3;
      }

      if (num_regions_within_threshold > parameters->max_num_regions) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, parameters->num_threads == 1 || read->get_sequence_id() == parameters->debug_read, FormatString("Stopping because num_regions_within_threshold > parameters->max_num_regions. (num_regions_within_threshold = %ld, parameters->max_num_regions = %ld)\n\n\n", num_regions_within_threshold, parameters->max_num_regions), "RunAlignment");
        return -4;
      }
    }

    return ret_evaluate_mappings;
  }

  return 0;
}

int GraphMap::EvaluateMappings_(bool evaluate_edit_distance, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters) {
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

//  best_entry->GenerateFPFilter(index_, read, *parameters);

  if (best_entry->get_fpfilter_cov_bases() < 0.50f ||      // 0.20f
      best_entry->get_fpfilter_query_len() < 0.75f ||      // 0.20f
      best_entry->get_fpfilter_std() < 0.20f)              // 0.20f
    return 3;

//  if (best_entry->alignment_score_cov_bases < 0.20f ||      // 0.20f
//      best_entry->alignment_score_query_len < 0.20f ||      // 0.20f
//      best_entry->alignment_score_std < 0.20f)              // 0.20f
//    return 3;

  if (best_entry->CalcDistanceRatioSuppress() > parameters->error_rate)
    return 4;

  if (evaluate_edit_distance == true) {
    /// TODO: Reimplement or remove this chunk.
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Deprecated functionality."));

//    if (best_entry->get_alignment_primary().edit_distance <= 0) {
//      int64_t alignment_position = 0, edit_distance = 0;
//      if (best_entry->get_region_data().is_split == false || parameters->is_reference_circular == false)
//        CalcEditDistanceLinear(MyersEditDistanceWrapper, read, indexes_[0], *parameters, best_entry, &alignment_position, &edit_distance);
//      else
//        CalcEditDistanceCircular(MyersEditDistanceWrapper, read, indexes_[0], *parameters, best_entry, &alignment_position, &edit_distance);
//
//      best_entry->get_alignment_primary().edit_distance = edit_distance;
//
//      // Calculate the alignment score for each alignment separately.
//      float NM_ratio = ((float) edit_distance) / ((float) read->get_sequence_length());
//      float alignment_score_float = (1.0f - sigmoid(NM_ratio, parameters->error_rate, 0.05));  // , (readlength * parameters.error_rate), 0.03f)) * 100.0f;
//      float covered_bases_ratio = (((float) best_entry->get_mapping_data().cov_bases_max) / ((float) read->get_sequence_length()));
//      float factor_covered_bases = (sigmoid(covered_bases_ratio, 0.05, 0.05));
//      float factor_query_length = (((float) (best_entry->get_mapping_data().query_coords.end - best_entry->get_mapping_data().query_coords.start)) / (((float) read->get_sequence_length()) * (1.0f - parameters->error_rate)));
//      if (factor_query_length > 1.0f)
//        factor_query_length = 1.0f;
//      alignment_score_float = alignment_score_float * factor_covered_bases * factor_query_length * 100;
//      int64_t alignment_score = (int64_t) round(alignment_score_float);
//
//    }
//
//    if (best_entry->get_alignment_primary().edit_distance > (parameters->error_rate * read->get_sequence_length()) || best_entry->get_alignment_primary().alignment_score <= 1)
//      return 5;
  }

  return 0;
}

int GraphMap::GenerateAlignments2_(MappingData *mapping_data, const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params) {
  if (mapping_data->intermediate_mappings.size() == 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("mapping_data->intermediate_mappings.size() == 0\n"), "GenerateAlignments_");
    return 1;
  }

  EvaluateMappings_(false, mapping_data, read, parameters);
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
    int ret_aln = AlignRegion(read, index, parameters, evalue_params, true, region_data);
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
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Alignment sane, but ret_aln != 0 (ret_aln = %d)!\n", ret_aln), "GenerateAlignments_");
        continue;
      } else {
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

//    /// TODO: This needs to happen in all cases! Currently only for testing purposes for anchored alignment.
//    if (parameters->alignment_algorithm == "anchor") {
//      /// Check if something went wrong and the read is unmapped.
//      if (evalue_left_part > ((double) 1e+300) || evalue_right_part > ((double) 1e+300) ||
//          (evalue_left_part > ((double) 1e0) && mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().cov_bases_query < 0.10f * read->get_sequence_length())) {
//
//        mapping_data->unmapped_reason += FormatString("__Evalue_is_too_large__evalue_left_part=%e__evalue_right_part=%e", evalue_left_part, evalue_right_part);
//        mapping_data->final_mapping_ptrs.at(i)->SetAligned(false);
//        mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().unmapped_reason = mapping_data->unmapped_reason;
//
//        /// Keep the output if alignment is insane for debug purposes.
//        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Alignment is insane!\n"), "GenerateAlignments_");
//      }
//    }

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

  return 0;
}

int GraphMap::GenerateAlignments_(MappingData *mapping_data, const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params) {
  if (mapping_data->intermediate_mappings.size() == 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("mapping_data->intermediate_mappings.size() == 0\n"), "GenerateAlignments_");
    return 1;
  }

  EvaluateMappings_(false, mapping_data, read, parameters);
  CollectFinalMappingsAndMapQ_(true, mapping_data, read, parameters);

  for (int64_t i = 0; i < mapping_data->final_mapping_ptrs.size(); i++) {
    if (mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().is_mapped == false) {
      mapping_data->final_mapping_ptrs.at(i)->SetAligned(false);
      mapping_data->unmapped_reason += FormatString("__Unaligned_because_get_mapping_data().is_mapped==false");
      mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().unmapped_reason = mapping_data->unmapped_reason;

      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().is_mapped == false\n"), "GenerateAlignments_");
      continue;
    }

    int64_t relative_position_left_part = 0;
    int64_t relative_position_right_part = 0;
    SeqOrientation orientation = kForward;
    int64_t reference_id = -1;
    int64_t position_ambiguity = 0;
    std::string cigar_left_part = "*";
    std::string cigar_right_part = "";

    int64_t num_eq_ops = 0, num_x_ops = 0, num_i_ops = 0, num_d_ops = 0;

    int64_t AS_left_part = 0, AS_right_part = 0;
    double evalue_left_part = 0, evalue_right_part = 0;
    int64_t nonclipped_length_left_part = 0, nonclipped_length_right_part = 0;
    // TODO: Promijeniti interface ove funkcije da ne vraca edit distance kao povratnu vrijednost, nego preko parametra, a da povratna vrijednost bude samo status o uspjehu.

    clock_t begin_clock = clock();
    int edit_distance = HybridRealignment(read, index, *parameters, mapping_data->final_mapping_ptrs.at(i), &relative_position_left_part, &cigar_left_part, &AS_left_part, &nonclipped_length_left_part, &relative_position_right_part, &cigar_right_part, &AS_right_part, &nonclipped_length_right_part, &orientation, &reference_id, &position_ambiguity, &num_eq_ops, &num_x_ops, &num_i_ops, &num_d_ops);
    clock_t end_clock = clock();
    double elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n+++++++++++++++++ Alignment elapsed time: %f sec.\n\n", elapsed_secs), "GenerateAlignments_");

    CalculateEValueDNA(AS_left_part, nonclipped_length_left_part, index->get_data_length_forward(), evalue_params, &evalue_left_part);
    CalculateEValueDNA(AS_right_part, nonclipped_length_right_part, index->get_data_length_forward(), evalue_params, &evalue_right_part);

    int64_t sum_of_error_ops = num_x_ops + num_i_ops + num_d_ops;
    float mismatch_rate = (((float) (num_x_ops + num_i_ops + num_d_ops)) / ((float) (num_eq_ops + num_x_ops + num_d_ops + num_i_ops)));
    float match_rate = (((float) num_eq_ops) / ((float) read->get_sequence_length()));
    float min_match_rate = (1.0f - parameters->error_rate);

    float covered_bases_ratio = (((float) mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().cov_bases_max) / ((float) read->get_sequence_length()));
    float factor_covered_bases = (sigmoid(covered_bases_ratio, 0.05, 0.05));
    float factor_query_length = (((float) (mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().query_coords.end - mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().query_coords.start)) / (((float) read->get_sequence_length()) * (1.0f - parameters->error_rate)));
    if (factor_query_length > 1.0f)
      factor_query_length = 1.0f;
    float factor_readlength = std::min(((float) read->get_sequence_length()) / 2000.0f, 1.0f);

    /// Check if something went wrong and the read is unmapped.
    if (edit_distance < 0) {
      mapping_data->unmapped_reason += FormatString("__HybridRealignment_returned_with_error__ret_value=%d", edit_distance);
      mapping_data->final_mapping_ptrs.at(i)->SetAligned(false);
      mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().unmapped_reason = mapping_data->unmapped_reason;

      /// Keep the output if alignment is insane for debug purposes.
      if (edit_distance != ALIGNMENT_NOT_SANE) {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Alignment sane, but edit_distance < 0!\n"), "GenerateAlignments_");
        continue;
      } else {
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Alignment is insane!\n"), "GenerateAlignments_");
      }
    }

    /// TODO: This needs to happen in all cases! Currently only for testing purposes for anchored alignment.
    if (parameters->alignment_algorithm == "anchor") {
      /// Check if something went wrong and the read is unmapped.
      if (evalue_left_part > ((double) 1e+300) || evalue_right_part > ((double) 1e+300) ||
          (evalue_left_part > ((double) 1e0) && mapping_data->final_mapping_ptrs.at(i)->get_mapping_data().cov_bases_query < 0.10f * read->get_sequence_length())) {

        mapping_data->unmapped_reason += FormatString("__Evalue_is_too_large__evalue_left_part=%e__evalue_right_part=%e", evalue_left_part, evalue_right_part);
        mapping_data->final_mapping_ptrs.at(i)->SetAligned(false);
        mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().unmapped_reason = mapping_data->unmapped_reason;

        /// Keep the output if alignment is insane for debug purposes.
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Alignment is insane!\n"), "GenerateAlignments_");
      }
    }

    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("cigar_right_part.size() = %d\n", cigar_right_part.size()), "GenerateAlignments_");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("evalue_left_part = %f\n", evalue_left_part), "GenerateAlignments_");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("mismatch_rate = %f\n", mismatch_rate), "GenerateAlignments_");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("match_rate = %f\n", match_rate), "GenerateAlignments_");
    }

    AlignmentResults primary_alignment;

    primary_alignment.is_aligned = true;
    primary_alignment.is_reverse = (orientation == kForward) ? false : true;
    primary_alignment.ref_start = relative_position_left_part; // + 1;
    primary_alignment.ref_end = relative_position_left_part; // + 1;
    primary_alignment.cigar = (orientation == kForward) ? (cigar_left_part) : (ReverseCigarString(cigar_left_part));

    primary_alignment.mapping_quality = mapping_data->mapping_quality;
    primary_alignment.edit_distance = edit_distance;
    primary_alignment.alignment_score = AS_left_part;
    primary_alignment.evalue = evalue_left_part;
    primary_alignment.num_secondary_alns = mapping_data->num_same_mappings;

    primary_alignment.num_eq_ops = num_eq_ops;
    primary_alignment.num_x_ops = num_x_ops;
    primary_alignment.num_i_ops = num_i_ops;
    primary_alignment.num_d_ops = num_d_ops;
    primary_alignment.nonclipped_length = nonclipped_length_left_part;

    mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().time_region_selection = mapping_data->time_region_selection;
    mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().time_mapping = mapping_data->time_mapping;
    mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().time_alignment = elapsed_secs;
    mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().time_region_seed_lookup = mapping_data->time_region_seed_lookup;
    mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().time_region_hitsort = mapping_data->time_region_hitsort;
    mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().time_region_conversion = mapping_data->time_region_conversion;

    if (parameters->evalue_threshold >= ((double) 0.0) && primary_alignment.evalue > parameters->evalue_threshold) {
      primary_alignment.is_aligned = false;
      mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().unmapped_reason += FormatString("_evalue=%f>%f", primary_alignment.evalue, parameters->evalue_threshold);
    }
    if (parameters->mapq_threshold >= 0 && primary_alignment.mapping_quality < parameters->mapq_threshold) {
      primary_alignment.is_aligned = false;
      mapping_data->final_mapping_ptrs.at(i)->get_mapping_metadata().unmapped_reason += FormatString("_mapq=%ld<%ld", primary_alignment.mapping_quality, parameters->mapq_threshold);
    }

    mapping_data->final_mapping_ptrs.at(i)->AddAlignmentData(primary_alignment);

    if (cigar_right_part.size() > 0) {
      AlignmentResults secondary_alignment;
      secondary_alignment = primary_alignment;
      secondary_alignment.is_aligned = true;
      secondary_alignment.cigar = (orientation == kForward) ? (cigar_right_part) : (ReverseCigarString(cigar_right_part));
      secondary_alignment.ref_start = relative_position_right_part; // + 1;
      secondary_alignment.ref_end = relative_position_right_part; // + 1;
      secondary_alignment.alignment_score = AS_right_part;
      secondary_alignment.evalue = evalue_right_part;
      secondary_alignment.nonclipped_length = nonclipped_length_right_part;

      if (primary_alignment.is_aligned == false)
        secondary_alignment.is_aligned = false;

      if (parameters->evalue_threshold >= ((double) 0.0) && secondary_alignment.evalue > parameters->evalue_threshold) {
        secondary_alignment.is_aligned = false;
//        secondary_alignment.get_alignment_primary().unmapped_reason += FormatString("_evalue=%f>%f", secondary_alignment.evalue, parameters->evalue_threshold);
      }
      if (parameters->mapq_threshold >= 0 && secondary_alignment.mapping_quality < parameters->mapq_threshold) {
        secondary_alignment.is_aligned = false;
//        primary_alignment.unmapped_reason += FormatString("_mapq=%ld<%ld", secondary_alignment.mapping_quality, parameters->mapq_threshold);
      }

      if (secondary_alignment.is_aligned == true)
        mapping_data->final_mapping_ptrs.at(i)->AddAlignmentData(secondary_alignment);
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
    if (ss.str().size() > 0)
      ss << "\n";

    if (parameters->outfmt == "sam") {
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateSAM((num_mapped_alignments == 0), parameters->verbose_sam_output);  // TODO: Don't make the first alignment primary by default, but the best one.
    } else if (parameters->outfmt == "afg") {
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateAFG();
    } else if (parameters->outfmt == "m5") {
      ss << mapping_data->final_mapping_ptrs.at(i)->GenerateM5((num_mapped_alignments == 0), parameters->verbose_sam_output);
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
      ret_aln_lines = GenerateUnmappedSamLine_(mapping_data, parameters->verbose_sam_output, read);
    } else if (parameters->outfmt == "afg") {
      // In AFG format there is no need to report 'unmapped' (or non-overlapping) reads).
    } else if (parameters->outfmt == "m5") {

    } else {  // Default to SAM output if the specified format is unknown.
      ret_aln_lines = GenerateUnmappedSamLine_(mapping_data, parameters->verbose_sam_output, read);
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

  int64_t num_kmers = first_entry->get_mapping_data().num_covering_kmers;
  int64_t min_num_kmers = (int64_t) std::ceil(((float) num_kmers) * (1.0f - parameters->margin_for_ambiguity));

  int64_t ambiguity = 1;
  int64_t num_equal_mappings = 1;
  for (int64_t i = 1; i < mapping_data->intermediate_mappings.size(); i++) {
    if (mapping_data->intermediate_mappings.at(i)->get_mapping_data().is_mapped == true) {
      // Find all the paths that are longer than 10% of the longest one.
      // if (mapping_data->intermediate_mappings.at(i).lcs_length >= min_score && mapping_data->intermediate_mappings.at(i).lcs_length <= max_score) { // >= allowed_score_difference_for_ambiguity) {
      if ((parameters->margin_for_ambiguity == 0.0f && mapping_data->intermediate_mappings.at(i)->get_mapping_data().num_covering_kmers == num_kmers) || (parameters->margin_for_ambiguity != 0.0f && mapping_data->intermediate_mappings.at(i)->get_mapping_data().num_covering_kmers >= min_num_kmers)) {
        // Check if the path is different than the best one - sometimes they can be near duplicates.
        if ((mapping_data->intermediate_mappings.at(i)->get_mapping_data().ref_coords.start != first_entry->get_mapping_data().ref_coords.start || mapping_data->intermediate_mappings.at(i)->get_mapping_data().ref_coords.end != first_entry->get_mapping_data().ref_coords.end)) {
          // Count them for ambiguity (used in mapping quality).
          if (generate_final_mapping_ptrs == true) {
            if (parameters->output_multiple_alignments) {
              mapping_data->final_mapping_ptrs.push_back((mapping_data->intermediate_mappings.at(i)));
              ambiguity += 1;
            }
          } else {
            ambiguity += 1;
          }

          // Count only the paths that are not the same length of the longest path. Might be useful for metagenomic,
          // as some genomes can have identical regions.
          if (mapping_data->intermediate_mappings.at(i)->get_mapping_data().num_covering_kmers == first_entry->get_mapping_data().num_covering_kmers) {
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

std::string GraphMap::GenerateUnmappedSamLine_(MappingData *mapping_data, int64_t verbose_sam_output, const SingleSequence *read) const {
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
  ss << read->get_data() << "\t";

  if (verbose_sam_output < 5 && read->get_quality_length() > 0 && read->get_quality() != NULL) {
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
    ss << "\tX3:Z:" << mapping_data->unmapped_reason;
  }

  return ss.str();
}
