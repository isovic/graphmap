/*
 * region_selection.cc
 *
 *  Created on: Mar 19, 2015
 *      Author: isovic
 */

#include "graphmap/graphmap.h"



int GraphMap::RegionSelection_(int64_t bin_size, MappingData* mapping_data, const Index* index, const Index* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {
  int64_t readlength = read->get_sequence_length();

  mapping_data->bin_size = bin_size;

  // Verbose all chromosomes.
//  if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//    for (int64_t i = 0; i < index->get_reference_starting_pos().size(); i++) {
//      if (i == 0)
//        ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("(forward)\n"), "[]");
//      else if (i == index->get_num_sequences_forward())
//        ErrorReporting::GetInstance().VerboseLog( VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("(reverse)\n"), "[]");
//      ErrorReporting::GetInstance().VerboseLog( VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[%ld] Chromosome: '%s', absolute start: %ld, length: %ld, starting bin: %ld\n", i, index->get_headers()[i % index->get_num_sequences_forward()].c_str(), index->get_reference_starting_pos()[i], index->get_reference_lengths()[i], index->get_reference_starting_pos()[i] / read->get_sequence_length()), "[]");
//    }
//  }

  ////////////////////////////////////////////////////
  ///// This part counts the occurances in bins. /////
  ////////////////////////////////////////////////////
  // Create bins for each chromosome (or reference sequence) separately.
  std::vector<std::vector<float> > bins_chromosome;
  std::vector<std::vector<int64_t> > last_update_chromosome;

  // Resize for the forward and reverse too.
  bins_chromosome.resize(index->get_num_sequences_forward() * 2);
  last_update_chromosome.resize(index->get_num_sequences_forward() * 2);
  // Resizing containers for each chromosome.
  for (int64_t i = 0; i < index->get_num_sequences_forward(); i++) {
    int64_t current_reference_length = index->get_reference_lengths()[i];
    int64_t current_num_bins = ceil(((float) current_reference_length) / ((float) bin_size));

    // Forward strand.
    bins_chromosome[i].resize(current_num_bins, 0);
    last_update_chromosome[i].resize(current_num_bins, 0);
    // Reverse strand.
    bins_chromosome[i + index->get_num_sequences_forward()].resize(current_num_bins, 0);
    last_update_chromosome[i + index->get_num_sequences_forward()].resize(current_num_bins, 0);
  }

  mapping_data->num_seeds_with_no_hits = 0;
  mapping_data->num_seeds_over_limit = 0;
  mapping_data->num_seeds_errors = 0;

  int64_t k = (int64_t) ((IndexSpacedHash *) index_)->get_shape_index_length();

  // Filling the bins with values, so we get an occurrence map.
  for (int64_t i = 0; i < (readlength - k + 1); i += parameters->kmer_step) {  // i++) {
    int8_t *seed = (int8_t *) &(read->get_data()[i]);
    uint64_t hits_start = 0, num_hits = 0;
    int64_t *hits = NULL;

    if (index_ != NULL) {
      int ret_search = index->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);

      // Check if there is too many hits (or too few).
      if (ret_search == 1) {
        mapping_data->num_seeds_with_no_hits += 1;

      } else if (ret_search == 2) {
        mapping_data->num_seeds_over_limit += 1;

//        if (parameters->max_num_hits > 0) {
//          if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//            ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
//            ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Skipping kmer, too many hits: num_hits = %ld, parameters->max_num_hits = %ld\n", num_hits, parameters->max_num_hits), "SelectRegionsWithHoughAndCircular");
//          }
//
//          if (index->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
//            if (hits) {
//              free(hits);
//            }
//            hits = NULL;
//          }
//          continue;
//        }
      } else if (ret_search > 2) {
        mapping_data->num_seeds_errors += 1;
      }

      // Counting kmers in regions of bin_size on the genome
      for (int64_t j = hits_start; j < (hits_start + num_hits); j++) {
        int64_t position = hits[j];

        int64_t x = i;          // Coordinate on the read.
        int64_t y = position;   // Coordinate on the reference.
        int64_t l = y - x;      // Hough projection on the reference.

        // Find the index of the reference that was hit. This also includes the reverse sequences.
        // Reverse sequences are considered the same as any other reference sequence.
        int64_t reference_index = index->RawPositionToReferenceIndexWithReverse(y);
        int64_t reference_starting_pos = index->get_reference_starting_pos()[reference_index];
        if (reference_index < 0) {
          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, y, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
          continue;
        }

        // Convert the absolute coordinates to local coordinates on the hit reference.
        int64_t y_local = y - reference_starting_pos;
        int64_t l_local = l - reference_starting_pos;

        // Compensate for sequence overhangs.
        if (l_local < 0 && parameters->is_reference_circular == false) {
          l_local = 0;
        }
        if (l_local < 0 && parameters->is_reference_circular == true) {
          l_local = index->get_reference_lengths()[reference_index] - 1;
        }

        // Calculate the index of the bin the position belongs to.
        int64_t position_bin = floor(((float) l_local) / ((float) bin_size));

        // We mark the last update with (i + 1) and not only i to avoid the default value of zero that has been set with vector initialization.
        if (parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)) {
//          ErrorReporting::GetInstance().VerboseLog(
//          VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)\n\n"), "SelectRegionsWithHoughAndCircular");
          continue;
        }

        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
        int64_t reference_end = index->get_reference_starting_pos()[reference_index] + index->get_reference_lengths()[reference_index];

        if (reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()) {
//          ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()\n\n"), "SelectRegionsWithHoughAndCircular");
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

    if (index_secondary != NULL) {
      int ret_search = index_secondary->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);

      // Check if there is too many hits (or too few).
      if (ret_search == 1) {
        mapping_data->num_seeds_with_no_hits += 1;

      } else if (ret_search == 2) {
        mapping_data->num_seeds_over_limit += 1;

//        if (parameters->max_num_hits > 0) {
//          if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//            ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n"), "[]");
//            ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Skipping kmer, too many hits: num_hits = %ld, parameters->max_num_hits = %ld\n", num_hits, parameters->max_num_hits), "SelectRegionsWithHoughAndCircular");
//          }
//
//          if (index_secondary->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
//            if (hits) {
//              free(hits);
//            }
//            hits = NULL;
//          }
//          continue;
//        }
      } else if (ret_search > 2) {
        mapping_data->num_seeds_errors += 1;
      }

      // Counting kmers in regions of bin_size on the genome
      for (int64_t j = hits_start; j < (hits_start + num_hits); j++) {
        int64_t position = hits[j];

        int64_t x = i;          // Coordinate on the read.
        int64_t y = position;   // Coordinate on the reference.
        int64_t l = y - x;      // Hough projection on the reference.

        // Find the index of the reference that was hit. This also includes the reverse sequences.
        // Reverse sequences are considered the same as any other reference sequence.
        int64_t reference_index = index_secondary->RawPositionToReferenceIndexWithReverse(y);
        int64_t reference_starting_pos = index_secondary->get_reference_starting_pos()[reference_index];
        if (reference_index < 0) {
          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, y, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
          continue;
        }

        // Convert the absolute coordinates to local coordinates on the hit reference.
        int64_t y_local = y - reference_starting_pos;
        int64_t l_local = l - reference_starting_pos;

        // Compensate for sequence overhangs.
        if (l_local < 0 && parameters->is_reference_circular == false) {
          l_local = 0;
        }
        if (l_local < 0 && parameters->is_reference_circular == true) {
          l_local = index_secondary->get_reference_lengths()[reference_index] - 1;
        }

        // Calculate the index of the bin the position belongs to.
        int64_t position_bin = floor(((float) l_local) / ((float) bin_size));

        // We mark the last update with (i + 1) and not only i to avoid the default value of zero that has been set with vector initialization.
        if (parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)) {
//          ErrorReporting::GetInstance().VerboseLog(
//          VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)\n\n"), "SelectRegionsWithHoughAndCircular");
          continue;
        }

        int64_t reference_start = index_secondary->get_reference_starting_pos()[reference_index];
        int64_t reference_end = index_secondary->get_reference_starting_pos()[reference_index] + index_secondary->get_reference_lengths()[reference_index];

        if (reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()) {
//          ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()\n\n"), "SelectRegionsWithHoughAndCircular");
          continue;
        }

        bins_chromosome[reference_index][position_bin] += 1.0f;
        last_update_chromosome[reference_index][position_bin] = (i + 1);
      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)

      if (index_secondary->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
        if (hits)
          free(hits);
        hits = NULL;
      }
      hits = NULL;
    }
  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");

  int64_t num_bins_above_zero = 0;
  for (int64_t i = 0; i < (index->get_num_sequences_forward() * 2); i++) {
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
  for (int64_t i = 0; i < (index->get_num_sequences_forward() * 2); i++) {
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
// Ovaj debug sam maknuo za brzinu kod profiliranja!
  if (parameters->verbose_level > 8 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Regions kmer count along reference:\n"), "OccuranceStatistics");
    // The tuple will contain: reference_id, bin_index, bin_count.
    for (int64_t i = 0; i < (index->get_num_sequences_forward()); i++) {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (forward) %s, bins_chromosome.size() = %ld, ", i, index->get_headers()[i].c_str(), bins_chromosome[i].size()), "[]");
      float max_bin_value_forward = 0;
      for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
        max_bin_value_forward = std::max(max_bin_value_forward, bins_chromosome[i][j]);
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_forward), "[]");
      std::stringstream ss_bins_forward;
      for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
        if (bins_chromosome[i][j] > 0) {
          ss_bins_forward << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i][j]);
          if (ss_bins_forward.str().size() > 120) {
            ss_bins_forward << "\n";
            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_forward.str(), "[]");
            ss_bins_forward.str("");
          }
        }
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_forward.str(), "[]");
      LogSystem::GetInstance().VerboseLog(
      VERBOSE_LEVEL_ALL_DEBUG,
                                               read->get_sequence_id() == parameters->debug_read, "\n\n", "[]");

      std::stringstream ss_bins_reverse;
      LogSystem::GetInstance().VerboseLog(
      VERBOSE_LEVEL_ALL_DEBUG,
                                               read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (reverse) %s, bins_chromosome.size() = %ld, ", i + index->get_num_sequences_forward(), index->get_headers()[i].c_str(), bins_chromosome[i + index->get_num_sequences_forward()].size()), "[]");
      float max_bin_value_reverse = 0;
      for (int64_t j = 0; j < bins_chromosome[i + index->get_num_sequences_forward()].size(); j++) {
        max_bin_value_reverse = std::max(max_bin_value_reverse, bins_chromosome[i + index->get_num_sequences_forward()][j]);
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_reverse), "[]");
      for (int64_t j = 0; j < bins_chromosome[i + index->get_num_sequences_forward()].size(); j++) {
        if (bins_chromosome[i + index->get_num_sequences_forward()][j] > 0) {
          ss_bins_reverse << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i + index->get_num_sequences_forward()][j]);
          if (ss_bins_reverse.str().size() > 120) {
            ss_bins_reverse << "\n";
            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_reverse.str(), "[]");
            ss_bins_reverse.str("");
          }
        }
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_reverse.str(), "[]");
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n\n", "[]");
    }
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
  }

  return 0;
}
