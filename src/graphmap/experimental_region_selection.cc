/*
 * experimental_region_selection.cc
 *
 *  Created on: Aug 18, 2015
 *      Author: isovic
 */

#include "graphmap/graphmap.h"

int GraphMap::RegionSelectionSpacedHashv2_(int64_t bin_size, MappingData* mapping_data, const IndexSpacedHash* index_primary, const IndexSpacedHash* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {
  int64_t readlength = read->get_sequence_length();

  mapping_data->bin_size = bin_size;

  std::vector<const IndexSpacedHash *> indexes = { index_primary, index_secondary };

  ////////////////////////////////////////////////////
  ///// This part counts the occurances in bins. /////
  ////////////////////////////////////////////////////
  // Create bins for each chromosome (or reference sequence) separately.
  std::vector<std::vector<float> > bins_chromosome;
  std::vector<std::vector<int64_t> > last_update_chromosome;

  // Resize for the forward and reverse too.
  bins_chromosome.resize(index_primary->get_num_sequences_forward() * 2);
  last_update_chromosome.resize(index_primary->get_num_sequences_forward() * 2);
  // Resizing containers for each chromosome.
  for (int64_t i = 0; i < index_primary->get_num_sequences_forward(); i++) {
    int64_t current_reference_length = index_primary->get_reference_lengths()[i];
    int64_t current_num_bins = ceil(((float) current_reference_length) / ((float) bin_size));

    // Forward strand.
    bins_chromosome[i].resize(current_num_bins, 0);
    last_update_chromosome[i].resize(current_num_bins, 0);
    // Reverse strand.
    bins_chromosome[i + index_primary->get_num_sequences_forward()].resize(current_num_bins, 0);
    last_update_chromosome[i + index_primary->get_num_sequences_forward()].resize(current_num_bins, 0);
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
//          printf ("[i = %ld] l_local = %ld\n", i, l_local);
//          fflush(stdout);

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
    }
  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");

  int64_t num_bins_above_zero = 0;
  for (int64_t i = 0; i < (index_primary->get_num_sequences_forward() * 2); i++) {
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
  for (int64_t i = 0; i < (index_primary->get_num_sequences_forward() * 2); i++) {
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
    for (int64_t i = 0; i < (index_primary->get_num_sequences_forward()); i++) {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (forward) %s, bins_chromosome.size() = %ld, ", i, index_primary->get_headers()[i].c_str(), bins_chromosome[i].size()), "[]");
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
                                          read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (reverse) %s, bins_chromosome.size() = %ld, ", i + index_primary->get_num_sequences_forward(), index_primary->get_headers()[i].c_str(), bins_chromosome[i + index_primary->get_num_sequences_forward()].size()), "[]");
      float max_bin_value_reverse = 0;
      for (int64_t j = 0; j < bins_chromosome[i + index_primary->get_num_sequences_forward()].size(); j++) {
        max_bin_value_reverse = std::max(max_bin_value_reverse, bins_chromosome[i + index_primary->get_num_sequences_forward()][j]);
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_reverse), "[]");
      for (int64_t j = 0; j < bins_chromosome[i + index_primary->get_num_sequences_forward()].size(); j++) {
        if (bins_chromosome[i + index_primary->get_num_sequences_forward()][j] > 0) {
          ss_bins_reverse << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i + index_primary->get_num_sequences_forward()][j]);
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

///// Handles a special case when bin_size = -1 as well, in which case the references will not be split into smaller regions.
int GraphMap::RegionSelectionSpacedHashFast_(int64_t bin_size, MappingData* mapping_data, const Index* index_primary, const Index* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {
  int64_t readlength = read->get_sequence_length();
  int64_t num_fwd_seqs = index_primary->get_num_sequences_forward();

  mapping_data->bin_size = bin_size;

  std::vector<const Index *> indexes = { index_primary, index_secondary };

  ////////////////////////////////////////////////////
  ///// This part counts the occurances in bins. /////
  ////////////////////////////////////////////////////
  // Create bins for each chromosome (or reference sequence) separately.
  std::vector<std::vector<float> > bins_chromosome;
  std::vector<std::vector<int64_t> > last_update_chromosome;

  // Resize for the forward and reverse too.
  bins_chromosome.resize(num_fwd_seqs * 2);
  last_update_chromosome.resize(num_fwd_seqs * 2);
  // Resizing containers for each chromosome.
  for (int64_t i = 0; i < num_fwd_seqs; i++) {
    int64_t current_reference_length = index_primary->get_reference_lengths()[i];
    int64_t current_num_bins = (bin_size > 0) ? ceil(((float) current_reference_length) / ((float) bin_size)) : 1;

    // Forward strand.
    bins_chromosome[i].resize(current_num_bins, 0);
    last_update_chromosome[i].resize(current_num_bins, 0);
    // Reverse strand.
    bins_chromosome[i + num_fwd_seqs].resize(current_num_bins, 0);
    last_update_chromosome[i + num_fwd_seqs].resize(current_num_bins, 0);
  }

  mapping_data->num_seeds_with_no_hits = 0;
  mapping_data->num_seeds_over_limit = 0;
  mapping_data->num_seeds_errors = 0;

  int64_t k = (int64_t) ((IndexSpacedHashFast *) index_)->get_shape_index_length();

  // Filling the bins with values, so we get an occurrence map.
  for (int64_t i = 0; i < (readlength - k + 1); i += parameters->kmer_step) {  // i++) {
    int8_t *seed = (int8_t *) &(read->get_data()[i]);
    uint64_t hits_start = 0, num_hits = 0;
    int64_t *hits = NULL;

    for (int64_t index_id = 0; index_id < indexes.size(); index_id++) {
      if (indexes[index_id] == NULL)
        continue;

      IndexSpacedHashFast *index = (IndexSpacedHashFast *) indexes[index_id];
      int ret_search = ((IndexSpacedHashFast *) index)->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);
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

//        int64_t reference_index = (int64_t) (position & MASK_REF_ID);
//        int64_t local_position = position >> 32;  // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
        int64_t local_position = (int64_t) (position & MASK_32_BIT);
        int64_t reference_index = (int64_t) ((uint64_t) position) >> 32;  // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);

        if (reference_index < 0) {
          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "3Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, local_position, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
          continue;
        }

        int64_t x = i;          // Coordinate on the read.
        int64_t y_local = local_position;
        int64_t l_local = y_local - x;

        if (parameters->alignment_approach == "overlapper") {
          if (index->get_headers()[reference_index % num_fwd_seqs] == ((std::string) read->get_header())) {
            continue;
          }
        }

        // Compensate for sequence overhangs.
        if (l_local < 0 && parameters->is_reference_circular == false) {
          l_local = 0;
        }
        if (l_local < 0 && parameters->is_reference_circular == true) {
          l_local = index->get_reference_lengths()[reference_index] - 1;
        }

        // Calculate the index of the bin the position belongs to.
        int64_t position_bin = (bin_size > 0) ? floor(((float) l_local) / ((float) bin_size)) : 0;
        if (reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()) {
          continue;
        }
        // We mark the last update with (i + 1) and not only i to avoid the default value of zero that has been set with vector initialization.
        if (parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)) {
          continue;
        }

        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
        int64_t reference_end = index->get_reference_starting_pos()[reference_index] + index->get_reference_lengths()[reference_index];

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
  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");

  int64_t num_bins_above_zero = 0;
  for (int64_t i = 0; i < (num_fwd_seqs * 2); i++) {
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
  for (int64_t i = 0; i < (num_fwd_seqs * 2); i++) {
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
  if (parameters->verbose_level > 8 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Regions kmer count along reference:\n"), "OccuranceStatistics");
    // The tuple will contain: reference_id, bin_index, bin_count.
    for (int64_t i = 0; i < (num_fwd_seqs); i++) {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (forward) %s, bins_chromosome.size() = %ld, ", i, index_primary->get_headers()[i].c_str(), bins_chromosome[i].size()), "[]");
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
                                          read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (reverse) %s, bins_chromosome.size() = %ld, ", i + num_fwd_seqs, index_primary->get_headers()[i].c_str(), bins_chromosome[i + num_fwd_seqs].size()), "[]");
      float max_bin_value_reverse = 0;
      for (int64_t j = 0; j < bins_chromosome[i + num_fwd_seqs].size(); j++) {
        max_bin_value_reverse = std::max(max_bin_value_reverse, bins_chromosome[i + num_fwd_seqs][j]);
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_reverse), "[]");
      for (int64_t j = 0; j < bins_chromosome[i + num_fwd_seqs].size(); j++) {
        if (bins_chromosome[i + num_fwd_seqs][j] > 0) {
          ss_bins_reverse << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i + num_fwd_seqs][j]);
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

/////// Handles a special case when bin_size = -1 as well, in which case the references will not be split into smaller regions.
//int GraphMap::ExperimentalRegionSelectionInverseDivision_(int64_t bin_size, MappingData* mapping_data, const Index* index_primary, const Index* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {
//  int64_t readlength = read->get_sequence_length();
//  int64_t num_fwd_seqs = index_primary->get_num_sequences_forward();
//
//  mapping_data->bin_size = bin_size;
//
//  float bin_size_inverse = (bin_size != 0) ? ((float) 1.0) / ((float) bin_size) : 0;
//
//  std::vector<const Index *> indexes = { index_primary, index_secondary };
//
//  ////////////////////////////////////////////////////
//  ///// This part counts the occurances in bins. /////
//  ////////////////////////////////////////////////////
//  // Create bins for each chromosome (or reference sequence) separately.
//  std::vector<std::vector<float> > bins_chromosome;
//  std::vector<std::vector<int64_t> > last_update_chromosome;
//
//  // Resize for the forward and reverse too.
//  bins_chromosome.resize(num_fwd_seqs * 2);
//  last_update_chromosome.resize(num_fwd_seqs * 2);
//  // Resizing containers for each chromosome.
//  for (int64_t i = 0; i < num_fwd_seqs; i++) {
//    int64_t current_reference_length = index_primary->get_reference_lengths()[i];
//    int64_t current_num_bins = (bin_size > 0) ? ceil(((float) current_reference_length) * bin_size_inverse) : 1;
//
//    // Forward strand.
//    bins_chromosome[i].resize(current_num_bins, 0);
//    last_update_chromosome[i].resize(current_num_bins, 0);
//    // Reverse strand.
//    bins_chromosome[i + num_fwd_seqs].resize(current_num_bins, 0);
//    last_update_chromosome[i + num_fwd_seqs].resize(current_num_bins, 0);
//  }
//
//  mapping_data->num_seeds_with_no_hits = 0;
//  mapping_data->num_seeds_over_limit = 0;
//  mapping_data->num_seeds_errors = 0;
//
//  int64_t k = (int64_t) ((IndexSpacedHashFast *) index_)->get_shape_index_length();
//
//  // Filling the bins with values, so we get an occurrence map.
//  for (int64_t i = 0; i < (readlength - k + 1); i += parameters->kmer_step) {  // i++) {
//    int8_t *seed = (int8_t *) &(read->get_data()[i]);
//    uint64_t hits_start = 0, num_hits = 0;
//    int64_t *hits = NULL;
//
//    for (int64_t index_id = 0; index_id < indexes.size(); index_id++) {
//      if (indexes[index_id] == NULL)
//        continue;
//
//      IndexSpacedHashFast *index = (IndexSpacedHashFast *) indexes[index_id];
//
////      std::vector<unsigned __int128> coded_hash_keys;
////      std::vector<int64_t> hash_keys, key_counts;
////      index->CalcAllKeysFromSequence(read, parameters->kmer_step, hash_keys, key_counts);
//
//      int ret_search = ((IndexSpacedHashFast *) index)->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);
//      // Check if there is too many hits (or too few).
//      if (ret_search == 1) {
//        mapping_data->num_seeds_with_no_hits += 1;
//      } else if (ret_search == 2) {
//        mapping_data->num_seeds_over_limit += 1;
//      } else if (ret_search > 2) {
//        mapping_data->num_seeds_errors += 1;
//      }
//
//      // Counting kmers in regions of bin_size on the genome
//      for (int64_t j = hits_start; j < (hits_start + num_hits); j++) {
//        int64_t position = hits[j];
//
//        int64_t local_position = (int64_t) (position & MASK_64_BIT);
//        int64_t reference_index = position >> 32;  // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
//
//        if (reference_index < 0) {
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "3Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, local_position, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
//          continue;
//        }
//
//        int64_t x = i;          // Coordinate on the read.
//        int64_t y_local = local_position;
//        int64_t l_local = y_local - x;
//
//        if (parameters->alignment_approach == "overlapper") {
//          if (index->get_headers()[reference_index % num_fwd_seqs] == ((std::string) read->get_header())) {
//            continue;
//          }
//        }
//
//        // Compensate for sequence overhangs.
//        if (l_local < 0 && parameters->is_reference_circular == false) {
//          l_local = 0;
//        }
//        if (l_local < 0 && parameters->is_reference_circular == true) {
//          l_local = index->get_reference_lengths()[reference_index] - 1;
//        }
//
//        // Calculate the index of the bin the position belongs to.
//        int64_t position_bin = (bin_size > 0) ? floor(((float) l_local) * bin_size_inverse) : 0;
//        if (reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()) {
//          continue;
//        }
//        // We mark the last update with (i + 1) and not only i to avoid the default value of zero that has been set with vector initialization.
//        if (parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)) {
//          continue;
//        }
//
//        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
//        int64_t reference_end = index->get_reference_starting_pos()[reference_index] + index->get_reference_lengths()[reference_index];
//
//        bins_chromosome[reference_index][position_bin] += 1.0f;
//        last_update_chromosome[reference_index][position_bin] = (i + 1);
//      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
//
//      if (index->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
//        if (hits)
//          free(hits);
//        hits = NULL;
//      }
//      hits = NULL;
//    }
//  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)
//
//  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");
//
//  int64_t num_bins_above_zero = 0;
//  for (int64_t i = 0; i < (num_fwd_seqs * 2); i++) {
//    for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
//      if (bins_chromosome[i][j] > 0.0f) {
//        num_bins_above_zero += 1;
//      }
//    }
//  }
//
//  // Convert the bins to a more compact form, which will be easier to sort.
//  // The tuple will contain: reference_id, bin_index, bin_count.
//  mapping_data->bins.clear();
//  mapping_data->bins.resize(num_bins_above_zero);
//  for (int64_t i = 0; i < (num_fwd_seqs * 2); i++) {
//    for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
//      if (bins_chromosome[i][j] > 0.0f) {
//        ChromosomeBin new_bin;
//        new_bin.reference_id = i;
//        new_bin.bin_id = j;
//        new_bin.bin_value = bins_chromosome[i][j];
//        mapping_data->bins.push_back(new_bin);
//      }
//    }
//  }
//
//  // Sort the bins in the descending order of bins_[i].bin_value;
//  std::sort(mapping_data->bins.begin(), mapping_data->bins.end(), bins_greater_than_key());
//
//  // Verbose all bin counts along each chromomsome.
//  if (parameters->verbose_level > 8 && read->get_sequence_id() == parameters->debug_read) {
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Regions kmer count along reference:\n"), "OccuranceStatistics");
//    // The tuple will contain: reference_id, bin_index, bin_count.
//    for (int64_t i = 0; i < (num_fwd_seqs); i++) {
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (forward) %s, bins_chromosome.size() = %ld, ", i, index_primary->get_headers()[i].c_str(), bins_chromosome[i].size()), "[]");
//      float max_bin_value_forward = 0;
//      for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
//        max_bin_value_forward = std::max(max_bin_value_forward, bins_chromosome[i][j]);
//      }
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_forward), "[]");
//      std::stringstream ss_bins_forward;
//      for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
//        if (bins_chromosome[i][j] > 0) {
//          ss_bins_forward << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i][j]);
//          if (ss_bins_forward.str().size() > 120) {
//            ss_bins_forward << "\n";
//            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_forward.str(), "[]");
//            ss_bins_forward.str("");
//          }
//        }
//      }
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_forward.str(), "[]");
//      LogSystem::GetInstance().VerboseLog(
//      VERBOSE_LEVEL_ALL_DEBUG,
//                                          read->get_sequence_id() == parameters->debug_read, "\n\n", "[]");
//
//      std::stringstream ss_bins_reverse;
//      LogSystem::GetInstance().VerboseLog(
//      VERBOSE_LEVEL_ALL_DEBUG,
//                                          read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (reverse) %s, bins_chromosome.size() = %ld, ", i + num_fwd_seqs, index_primary->get_headers()[i].c_str(), bins_chromosome[i + num_fwd_seqs].size()), "[]");
//      float max_bin_value_reverse = 0;
//      for (int64_t j = 0; j < bins_chromosome[i + num_fwd_seqs].size(); j++) {
//        max_bin_value_reverse = std::max(max_bin_value_reverse, bins_chromosome[i + num_fwd_seqs][j]);
//      }
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_reverse), "[]");
//      for (int64_t j = 0; j < bins_chromosome[i + num_fwd_seqs].size(); j++) {
//        if (bins_chromosome[i + num_fwd_seqs][j] > 0) {
//          ss_bins_reverse << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i + num_fwd_seqs][j]);
//          if (ss_bins_reverse.str().size() > 120) {
//            ss_bins_reverse << "\n";
//            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_reverse.str(), "[]");
//            ss_bins_reverse.str("");
//          }
//        }
//      }
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ss_bins_reverse.str(), "[]");
//      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n\n", "[]");
//    }
//    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, "\n", "[]");
//  }
//
//  return 0;
//}

///// Handles a special case when bin_size = -1 as well, in which case the references will not be split into smaller regions.
int GraphMap::ExperimentalRegionSelection_(int64_t bin_size, MappingData* mapping_data, const Index* index_primary, const Index* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {
  int64_t readlength = read->get_sequence_length();
  int64_t num_fwd_seqs = index_primary->get_num_sequences_forward();

  mapping_data->bin_size = bin_size;

  float bin_size_inverse = (bin_size != 0) ? ((float) 1.0) / ((float) bin_size) : 0;

  std::vector<const Index *> indexes = { index_primary, index_secondary };

  ////////////////////////////////////////////////////
  ///// This part counts the occurances in bins. /////
  ////////////////////////////////////////////////////
  // Create bins for each chromosome (or reference sequence) separately.
  std::vector<std::vector<float> > bins_chromosome;
  std::vector<std::vector<int64_t> > last_update_chromosome;

  // Resize for the forward and reverse too.
  bins_chromosome.resize(num_fwd_seqs * 2);
  last_update_chromosome.resize(num_fwd_seqs * 2);
  // Resizing containers for each chromosome.
  for (int64_t i = 0; i < num_fwd_seqs; i++) {
    int64_t current_reference_length = index_primary->get_reference_lengths()[i];
    int64_t current_num_bins = (bin_size > 0) ? ceil(((float) current_reference_length) * bin_size_inverse) : 1;

    // Forward strand.
    bins_chromosome[i].resize(current_num_bins, 0);
    last_update_chromosome[i].resize(current_num_bins, 0);
    // Reverse strand.
    bins_chromosome[i + num_fwd_seqs].resize(current_num_bins, 0);
    last_update_chromosome[i + num_fwd_seqs].resize(current_num_bins, 0);
  }

  mapping_data->num_seeds_with_no_hits = 0;
  mapping_data->num_seeds_over_limit = 0;
  mapping_data->num_seeds_errors = 0;

  std::vector<SeedHit3> seed_hits;

  for (int64_t index_id = 0; index_id < indexes.size(); index_id++) {
    if (indexes[index_id] == NULL)
      continue;

    IndexSpacedHashFast *index = (IndexSpacedHashFast *) indexes[index_id];

    clock_t begin_clock = clock();
    std::vector<unsigned __int128> coded_hash_keys;
    std::vector<int64_t> hash_keys, key_counts, seed_hits_y, seed_hits_x;
    index->CalcAllKeysFromSequence(read, parameters->kmer_step, hash_keys, key_counts);
    index->LookUpHashKeys(bin_size, read, hash_keys, key_counts, seed_hits);
    clock_t end_clock = clock();
    double elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n+++++++++++++++++ [Current index_id = %ld] New region selection time: %f sec.\n\n", index_id, elapsed_secs), "ProcessRead");
  }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)

  std::sort(seed_hits.begin(), seed_hits.end(), seed_hit3_compare());

  // Counting kmers in regions of bin_size on the genome
  for (int64_t j = 0; j < seed_hits.size(); j++) {
    int64_t reference_index = seed_hits[j].ref_id;
    int64_t x = seed_hits[j].x;
    int64_t y_local = seed_hits[j].y;
    int64_t l_local = y_local - x;

    if (reference_index < 0) {
      //        LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "3Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, y_local, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
      continue;
    }

    if (parameters->alignment_approach == "overlapper") {
      if (index_primary->get_headers()[reference_index % num_fwd_seqs] == ((std::string) read->get_header())) {
        continue;
      }
    }

    // Compensate for sequence overhangs.
    if (l_local < 0 && parameters->is_reference_circular == false) {
      l_local = 0;
    }
    if (l_local < 0 && parameters->is_reference_circular == true) {
      l_local = index_primary->get_reference_lengths()[reference_index] - 1;
    }

    // Calculate the index of the bin the position belongs to.
    int64_t position_bin = (bin_size > 0) ? floor(((float) l_local) * bin_size_inverse) : 0;
    if (reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()) {
      continue;
    }
    // We mark the last update with (i + 1) and not only i to avoid the default value of zero that has been set with vector initialization.
    if (parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (x + 1)) {
      continue;
    }

    int64_t reference_start = index_primary->get_reference_starting_pos()[reference_index];
    int64_t reference_end = index_primary->get_reference_starting_pos()[reference_index] + index_primary->get_reference_lengths()[reference_index];

    bins_chromosome[reference_index][position_bin] += 1.0f;
    last_update_chromosome[reference_index][position_bin] = (x + 1);

    //    printf ("Tu sam 1!\n");
    //    fflush(stdout);
    //    exit(1);
    //    break;
  }

  int64_t k = (int64_t) ((IndexSpacedHashFast *) index_)->get_shape_index_length();

//  // Filling the bins with values, so we get an occurrence map.
//  for (int64_t i = 0; i < (readlength - k + 1); i += parameters->kmer_step) {  // i++) {
//    int8_t *seed = (int8_t *) &(read->get_data()[i]);
//    uint64_t hits_start = 0, num_hits = 0;
//    int64_t *hits = NULL;
//
//    for (int64_t index_id = 0; index_id < indexes.size(); index_id++) {
//      if (indexes[index_id] == NULL)
//        continue;
//
//      IndexSpacedHashFast *index = (IndexSpacedHashFast *) indexes[index_id];
//
//      int ret_search = ((IndexSpacedHashFast *) index)->FindAllRawPositionsOfSeed(seed, k, parameters->max_num_hits, &hits, &hits_start, &num_hits);
//      // Check if there is too many hits (or too few).
//      if (ret_search == 1) {
//        mapping_data->num_seeds_with_no_hits += 1;
//      } else if (ret_search == 2) {
//        mapping_data->num_seeds_over_limit += 1;
//      } else if (ret_search > 2) {
//        mapping_data->num_seeds_errors += 1;
//      }
//
//      // Counting kmers in regions of bin_size on the genome
//      for (int64_t j = hits_start; j < (hits_start + num_hits); j++) {
//        int64_t position = hits[j];
//
//        int64_t local_position = (int64_t) (position & MASK_32_BIT);
//        int64_t reference_index = position >> 32;  // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
//
//        if (reference_index < 0) {
//          LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "3Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, local_position, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
//          continue;
//        }
//
//        int64_t x = i;          // Coordinate on the read.
//        int64_t y_local = local_position;
//        int64_t l_local = y_local - x;
//
//        if (parameters->alignment_approach == "overlapper") {
//          if (index->get_headers()[reference_index % num_fwd_seqs] == ((std::string) read->get_header())) {
//            continue;
//          }
//        }
//
//        // Compensate for sequence overhangs.
//        if (l_local < 0 && parameters->is_reference_circular == false) {
//          l_local = 0;
//        }
//        if (l_local < 0 && parameters->is_reference_circular == true) {
//          l_local = index->get_reference_lengths()[reference_index] - 1;
//        }
//
//        // Calculate the index of the bin the position belongs to.
//        int64_t position_bin = (bin_size > 0) ? floor(((float) l_local) * bin_size_inverse) : 0;
//        if (reference_index >= bins_chromosome.size() || position_bin >= bins_chromosome[reference_index].size()) {
//          continue;
//        }
//        // We mark the last update with (i + 1) and not only i to avoid the default value of zero that has been set with vector initialization.
//        if (parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)) {
//          continue;
//        }
//
//        int64_t reference_start = index->get_reference_starting_pos()[reference_index];
//        int64_t reference_end = index->get_reference_starting_pos()[reference_index] + index->get_reference_lengths()[reference_index];
//
//        bins_chromosome[reference_index][position_bin] += 1.0f;
//        last_update_chromosome[reference_index][position_bin] = (i + 1);
//      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
//
//      if (index->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
//        if (hits)
//          free(hits);
//        hits = NULL;
//      }
//      hits = NULL;
//    }
//  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");

  int64_t num_bins_above_zero = 0;
  for (int64_t i = 0; i < (num_fwd_seqs * 2); i++) {
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
  for (int64_t i = 0; i < (num_fwd_seqs * 2); i++) {
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
  if (parameters->verbose_level > 8 && read->get_sequence_id() == parameters->debug_read) {
    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Regions kmer count along reference:\n"), "OccuranceStatistics");
    // The tuple will contain: reference_id, bin_index, bin_count.
    for (int64_t i = 0; i < (num_fwd_seqs); i++) {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (forward) %s, bins_chromosome.size() = %ld, ", i, index_primary->get_headers()[i].c_str(), bins_chromosome[i].size()), "[]");
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
                                          read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (reverse) %s, bins_chromosome.size() = %ld, ", i + num_fwd_seqs, index_primary->get_headers()[i].c_str(), bins_chromosome[i + num_fwd_seqs].size()), "[]");
      float max_bin_value_reverse = 0;
      for (int64_t j = 0; j < bins_chromosome[i + num_fwd_seqs].size(); j++) {
        max_bin_value_reverse = std::max(max_bin_value_reverse, bins_chromosome[i + num_fwd_seqs][j]);
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_reverse), "[]");
      for (int64_t j = 0; j < bins_chromosome[i + num_fwd_seqs].size(); j++) {
        if (bins_chromosome[i + num_fwd_seqs][j] > 0) {
          ss_bins_reverse << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i + num_fwd_seqs][j]);
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















int GraphMap::RegionSelectionSpacedHashFastv2_(int64_t bin_size, MappingData* mapping_data, const IndexSpacedHashFast* index_primary, const IndexSpacedHashFast* index_secondary, const SingleSequence* read, const ProgramParameters* parameters) {
  int64_t readlength = read->get_sequence_length();

  mapping_data->bin_size = bin_size;

  float bin_size_inverse = 1.0f / ((float) bin_size);

  std::vector<const IndexSpacedHashFast *> indexes = { index_primary, index_secondary };

  ////////////////////////////////////////////////////
  ///// This part counts the occurances in bins. /////
  ////////////////////////////////////////////////////
  // Create bins for each chromosome (or reference sequence) separately.
  std::vector<std::vector<float> > bins_chromosome;
  std::vector<std::vector<int64_t> > last_update_chromosome;

  // Resize for the forward and reverse too.
  bins_chromosome.resize(index_primary->get_num_sequences_forward() * 2);
  last_update_chromosome.resize(index_primary->get_num_sequences_forward() * 2);
  // Resizing containers for each chromosome.
  for (int64_t i = 0; i < index_primary->get_num_sequences_forward(); i++) {
    int64_t current_reference_length = index_primary->get_reference_lengths()[i];
    int64_t current_num_bins = ceil(((float) current_reference_length) * bin_size_inverse);

    // Forward strand.
    bins_chromosome[i].resize(current_num_bins, 0);
    last_update_chromosome[i].resize(current_num_bins, 0);
    // Reverse strand.
    bins_chromosome[i + index_primary->get_num_sequences_forward()].resize(current_num_bins, 0);
    last_update_chromosome[i + index_primary->get_num_sequences_forward()].resize(current_num_bins, 0);
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

          int64_t x = i;          // Coordinate on the read.
//          int64_t y = position;   // Coordinate on the reference.
//          int64_t l = y - x;      // Hough projection on the reference.

          // Find the index of the reference that was hit. This also includes the reverse sequences.
          // Reverse sequences are considered the same as any other reference sequence.
//          int64_t reference_index = index->RawPositionToReferenceIndexWithReverse(y);
//          int64_t reference_starting_pos = index->get_reference_starting_pos()[reference_index];
          if (reference_index < 0) {
            LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, local_position, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
            continue;
          }

          // Convert the absolute coordinates to local coordinates on the hit reference.
          int64_t y_local = local_position;
          int64_t l_local = y_local - x;
//          printf ("[i = %ld] l_local = %ld\n", i, l_local);
//          fflush(stdout);

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
          if (parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)) {
            //          ErrorReporting::GetInstance().VerboseLog(
            //          VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "parameters->skip_multiple_kmers_per_bin == true && last_update_chromosome[reference_index][position_bin] == (i + 1)\n\n"), "SelectRegionsWithHoughAndCircular");
            continue;
          }

//          int64_t reference_start = index->get_reference_starting_pos()[reference_index];
//          int64_t reference_end = index->get_reference_starting_pos()[reference_index] + index->get_reference_lengths()[reference_index];

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
    }
  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");

  int64_t num_bins_above_zero = 0;
  for (int64_t i = 0; i < (index_primary->get_num_sequences_forward() * 2); i++) {
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
  for (int64_t i = 0; i < (index_primary->get_num_sequences_forward() * 2); i++) {
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
    for (int64_t i = 0; i < (index_primary->get_num_sequences_forward()); i++) {
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (forward) %s, bins_chromosome.size() = %ld, ", i, index_primary->get_headers()[i].c_str(), bins_chromosome[i].size()), "[]");
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
                                          read->get_sequence_id() == parameters->debug_read, FormatString("[ref %ld] (reverse) %s, bins_chromosome.size() = %ld, ", i + index_primary->get_num_sequences_forward(), index_primary->get_headers()[i].c_str(), bins_chromosome[i + index_primary->get_num_sequences_forward()].size()), "[]");
      float max_bin_value_reverse = 0;
      for (int64_t j = 0; j < bins_chromosome[i + index_primary->get_num_sequences_forward()].size(); j++) {
        max_bin_value_reverse = std::max(max_bin_value_reverse, bins_chromosome[i + index_primary->get_num_sequences_forward()][j]);
      }
      LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("max_bin_value = %.2f\n", max_bin_value_reverse), "[]");
      for (int64_t j = 0; j < bins_chromosome[i + index_primary->get_num_sequences_forward()].size(); j++) {
        if (bins_chromosome[i + index_primary->get_num_sequences_forward()][j] > 0) {
          ss_bins_reverse << FormatString("[%4ld] = %.2f\t", j, bins_chromosome[i + index_primary->get_num_sequences_forward()][j]);
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
