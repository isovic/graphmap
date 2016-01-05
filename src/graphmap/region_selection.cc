/*
 * region_selection.cc
 *
 *  Created on: Mar 19, 2015
 *      Author: isovic
 */

#include "graphmap/graphmap.h"

int GraphMap::RegionSelection_(int64_t bin_size, MappingData* mapping_data, const std::vector<Index *> indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  if (indexes.size() == 0 || indexes[0] == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "No reference indexes are specified."));
  }

  int64_t readlength = read->get_sequence_length();
  int64_t num_fwd_seqs = indexes[0]->get_num_sequences_forward();
  bool self_overlap = (parameters->alignment_approach == "overlapper" && parameters->reference_path == parameters->reads_path);

  mapping_data->bin_size = bin_size;

  float bin_size_inverse = (bin_size > 0) ? (1.0f / ((float) bin_size)) : (0.0f);

  ////////////////////////////////////////////////////
  ///// This part prepares the bins. /////
  ////////////////////////////////////////////////////
  // Create bins for each chromosome (or reference sequence) separately.
  std::vector<std::vector<float> > bins_chromosome;
  std::vector<std::vector<int64_t> > last_update_chromosome;

  // Resize for the forward and reverse too.
  bins_chromosome.resize(indexes[0]->get_num_sequences_forward() * 2);
  last_update_chromosome.resize(indexes[0]->get_num_sequences_forward() * 2);
  // Resizing containers for each chromosome.
  for (int64_t i = 0; i < indexes[0]->get_num_sequences_forward(); i++) {
    int64_t current_reference_length = indexes[0]->get_reference_lengths()[i];
    int64_t current_num_bins = ceil(((float) current_reference_length) * bin_size_inverse) + 1;
    // Forward strand.
    bins_chromosome[i].resize(current_num_bins, 0);
    last_update_chromosome[i].resize(current_num_bins, 0);
    // Reverse strand.
    bins_chromosome[i + indexes[0]->get_num_sequences_forward()].resize(current_num_bins, 0);
    last_update_chromosome[i + indexes[0]->get_num_sequences_forward()].resize(current_num_bins, 0);
  }

  mapping_data->num_seeds_with_no_hits = 0;
  mapping_data->num_seeds_over_limit = 0;
  mapping_data->num_seeds_errors = 0;

  int64_t k = (int64_t) ((IndexSpacedHash *) indexes[0])->get_shape_index_length();

  ////////////////////////////////////////////////////
  ///// This part counts the occurances in bins. /////
  ////////////////////////////////////////////////////
  // Filling the bins with values, so we get an occurrence map.
  for (int64_t i = 0; i < (readlength - k + 1); i += parameters->kmer_step) {  // i++) {
    int8_t *seed = (int8_t *) &(read->get_data()[i]);
    uint64_t hits_start = 0, num_hits = 0;
    int64_t *hits = NULL;

    for (int64_t index_id = 0; index_id < indexes.size(); index_id++) {
      IndexSpacedHash *index = (IndexSpacedHash *) indexes[index_id];

      if (index != NULL) {
        int ret_search = index->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);

        // Check if there is too many hits (or too few).
        if (ret_search == 1) {
          mapping_data->num_seeds_with_no_hits += 1;
        } else if (ret_search == 2) {
          mapping_data->num_seeds_over_limit += 1;
        } else if (ret_search > 2) {
          mapping_data->num_seeds_errors += 1;
        }

        // Counting kmers in regions of bin_size on the genome
        for (int64_t j = hits_start; j < (hits_start + num_hits); j++) {
          int64_t position = hits[j];
          int64_t local_position = (int64_t) (((uint64_t) position) & MASK_32_BIT);
          int64_t reference_index = (int64_t) (((uint64_t) position) >> 32);  // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);

          if (self_overlap == true &&
              (reference_index % num_fwd_seqs) == read->get_sequence_id()) {
            continue;
          }

          if (reference_index < 0) {
            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, local_position, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
            continue;
          }

          // Convert the absolute coordinates to local coordinates on the hit reference.
          int64_t x = i;          // Coordinate on the read.
          int64_t y_local = local_position;
          int64_t l_local = y_local - x;

          // Compensate for sequence overhangs.
          if (l_local < 0 && parameters->is_reference_circular == false) {
            l_local = 0;
          }
          if (l_local < 0 && parameters->is_reference_circular == true) {
            l_local = index->get_reference_lengths()[reference_index] - 1;
          }

          // Calculate the index of the bin the position belongs to.
          int64_t position_bin = floor(((float) l_local) * bin_size_inverse);

          // We mark the last update with (i + 1) and not only i to avoid the default value of zero that has been set with vector initialization.
          if (parameters->skip_multiple_kmers_per_bin == true &&
              last_update_chromosome[reference_index][position_bin] == (i + 1)) {
            continue;
          }
          if (reference_index >= bins_chromosome.size() ||
              position_bin >= bins_chromosome[reference_index].size()) {
            continue;
          }

          bins_chromosome[reference_index][position_bin] += 1.0f;
          last_update_chromosome[reference_index][position_bin] = (i + 1);
        }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)

        if (index->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
          if (hits)
            free(hits);
          hits = NULL;
        }
        hits = NULL;
      }
    }
  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");

  int64_t num_bins_above_zero = 0;
  for (int64_t i = 0; i < (indexes[0]->get_num_sequences_forward() * 2); i++) {
    for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
      if (bins_chromosome[i][j] > 0.0f) {
        num_bins_above_zero += 1;
      }
    }
  }

  // Convert the bins to a more compact form, which will be easier to sort.
  // The tuple will contain: reference_id, bin_index, bin_count.
  mapping_data->bins.clear();
  mapping_data->bins.resize(num_bins_above_zero);
  for (int64_t i = 0; i < (indexes[0]->get_num_sequences_forward() * 2); i++) {
    for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
      if (bins_chromosome[i][j] > 0.0f) {
        ChromosomeBin new_bin;
        new_bin.reference_id = i;
        new_bin.bin_id = j;
        new_bin.bin_value = bins_chromosome[i][j];
        mapping_data->bins.push_back(new_bin);
      }
    }
  }

  // Sort the bins in the descending order of bins_[i].bin_value;
  std::sort(mapping_data->bins.begin(), mapping_data->bins.end(), bins_greater_than_key());

  // Verbose all bin counts along each chromomsome.
  // This debug part below is commented out for speed during debugging/profiling.
//  if (parameters->verbose_level > 8 && read->get_sequence_id() == parameters->debug_read) {
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Regions kmer count along reference:\n"), "OccuranceStatistics");
//    // The tuple will contain: reference_id, bin_index, bin_count.
//    for (int64_t i = 0; i < (indexes[0]->get_num_sequences_forward()); i++) {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (forward) %s, bins_chromosome.size() = %ld, ", i, indexes[0]->get_headers()[i].c_str(), bins_chromosome[i].size()), "[]");
//      float max_bin_value_forward = 0;
//      for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
//        max_bin_value_forward = std::max(max_bin_value_forward, bins_chromosome[i][j]);
//      }
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_forward), "[]");
//      std::stringstream ss_bins_forward;
//      for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
//        if (bins_chromosome[i][j] > 0) {
//          ss_bins_forward << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i][j]);
//          if (ss_bins_forward.str().size() > 120) {
//            ss_bins_forward << "\n";
//            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_forward.str(), "[]");
//            ss_bins_forward.str("");
//          }
//        }
//      }
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_forward.str(), "[]");
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n\n", "[]");
//
//      std::stringstream ss_bins_reverse;
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (reverse) %s, bins_chromosome.size() = %ld, ", i + indexes[0]->get_num_sequences_forward(), indexes[0]->get_headers()[i].c_str(), bins_chromosome[i + indexes[0]->get_num_sequences_forward()].size()), "[]");
//      float max_bin_value_reverse = 0;
//      for (int64_t j = 0; j < bins_chromosome[i + indexes[0]->get_num_sequences_forward()].size(); j++) {
//        max_bin_value_reverse = std::max(max_bin_value_reverse, bins_chromosome[i + indexes[0]->get_num_sequences_forward()][j]);
//      }
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_reverse), "[]");
//      for (int64_t j = 0; j < bins_chromosome[i + indexes[0]->get_num_sequences_forward()].size(); j++) {
//        if (bins_chromosome[i + indexes[0]->get_num_sequences_forward()][j] > 0) {
//          ss_bins_reverse << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i + indexes[0]->get_num_sequences_forward()][j]);
//          if (ss_bins_reverse.str().size() > 120) {
//            ss_bins_reverse << "\n";
//            LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_reverse.str(), "[]");
//            ss_bins_reverse.str("");
//          }
//        }
//      }
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_reverse.str(), "[]");
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n\n", "[]");
//    }
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
//  }

  return 0;
}
