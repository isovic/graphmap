/*
 * experimental.cc
 *
 *  Created on: Jan 20, 2016
 *      Author: isovic
 */

#include "graphmap/graphmap.h"

typedef unsigned __int128 uint128_t;

int GraphMap::RegionSelectionNoBins_(int64_t bin_size, MappingData* mapping_data, const std::vector<Index *> indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  if (indexes.size() == 0 || indexes[0] == NULL) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "No reference indexes are specified."));
  }

  clock_t begin_clock = clock();

  int64_t readlength = read->get_sequence_length();
  int64_t num_fwd_seqs = indexes[0]->get_num_sequences_forward();
  bool self_overlap = (parameters->alignment_approach == "overlapper" && parameters->reference_path == parameters->reads_path);

  float bin_size_inverse = (bin_size > 0) ? (1.0f / ((float) bin_size)) : (0.0f);
  int64_t k = (int64_t) ((IndexSpacedHash *) indexes[0])->get_shape_index_length();

  mapping_data->bin_size = bin_size;
  mapping_data->num_seeds_with_no_hits = 0;
  mapping_data->num_seeds_over_limit = 0;
  mapping_data->num_seeds_errors = 0;

  std::vector<uint128_t> hit_coords;
  hit_coords.reserve(100000);

  for (int64_t index_id = 0; index_id < indexes.size(); index_id++) {
    IndexSpacedHash *index = (IndexSpacedHash *) indexes[index_id];
    if (index == NULL) { continue; }

    for (int64_t i = 0; i < (readlength - k + 1); i += parameters->kmer_step) {
      int8_t *seed = (int8_t *) &(read->get_data()[i]);
      uint64_t hits_start = 0, num_hits = 0;
      int64_t *hits = NULL;

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
        // Convert the absolute coordinates to local coordinates on the hit reference.
        int64_t x = i;          // Coordinate on the read.
        int64_t y = local_position;
        int64_t l = y - x;
        // Compensate for sequence overhangs.
        if (l < 0 && parameters->is_reference_circular == false) { l = 0; }
        else if (l < 0 && parameters->is_reference_circular == true) { l = index->get_reference_lengths()[reference_index] - 1; }
        // Calculate the index of the bin the position belongs to.
        int64_t position_bin = floor(((float) l) * bin_size_inverse);

        if (self_overlap == true && (reference_index % num_fwd_seqs) == read->get_sequence_id()) {
          continue;
        }

        if (reference_index < 0) {
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, local_position, j, hits_start, num_hits), "SelectRegionsWithHoughAndCircular");
          continue;
        }

        uint128_t hit_coord = (((uint128_t) reference_index) & ((uint128_t) 0x00000000FFFFFFFF)) << 64 |
                              (((uint128_t) position_bin) & ((uint128_t) 0x00000000FFFFFFFF)) << 32 |
                              (((uint128_t) i) & ((uint128_t) 0x00000000FFFFFFFF)) << 0;
        hit_coords.push_back(hit_coord);

      }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)


      if (index->IsManualCleanupRequired("FindAllRawPositionsOfSeed") == 0) {
        if (hits)
          free(hits);
        hits = NULL;
      }
      hits = NULL;

    }
  }

  double elapsed_secs_lookup = double(clock() - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_region_seed_lookup = elapsed_secs_lookup;
  begin_clock = clock();

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("Finished fetching seed hits. Sorting %ld hits.\n", hit_coords.size()), std::string(__FUNCTION__));

  std::sort(hit_coords.begin(), hit_coords.end());

  double elapsed_secs_hitsort = double(clock() - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_region_hitsort = elapsed_secs_hitsort;
  begin_clock = clock();

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld, hit_coords.size() = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit, hit_coords.size()), std::string(__FUNCTION__));

  // Convert the bins to a more compact form. Count occurrences of all bins, and discard the same hits (same query coordinate and same hit position (rasterized to a bin)).
  // The tuple will contain: reference_id, bin_index, bin_count.
  mapping_data->bins.clear();
  mapping_data->bins.reserve(10000);

  if (hit_coords.size() > 0) {
    ChromosomeBin new_bin;
    new_bin.reference_id = (hit_coords[0] >> 64) & ((uint128_t) 0x00000000FFFFFFFF);
    new_bin.bin_id = (hit_coords[0] >> 32) & ((uint128_t) 0x00000000FFFFFFFF);
    new_bin.bin_value = 1;
    mapping_data->bins.push_back(new_bin);
  }
  int64_t prev_query_pos = 0;
  for (int64_t i = 1; i < hit_coords.size(); i++) {
    int64_t ref_id = (hit_coords[i] >> 64) & ((uint128_t) 0x00000000FFFFFFFF);
    int64_t bin_id = (hit_coords[i] >> 32) & ((uint128_t) 0x00000000FFFFFFFF);
    int64_t query_pos = (hit_coords[i] >> 0) & ((uint128_t) 0x00000000FFFFFFFF);
    if (bin_id != mapping_data->bins.back().bin_id || ref_id != mapping_data->bins.back().reference_id) {
      ChromosomeBin new_bin;
      new_bin.reference_id = ref_id;
      new_bin.bin_id = bin_id;
      new_bin.bin_value = 1;
      mapping_data->bins.push_back(new_bin);
    } else if (i == 0 || (i > 0 && query_pos != prev_query_pos)) {
      mapping_data->bins.back().bin_value += 1;
    }
    prev_query_pos = query_pos;
  }

  // Sort the bins in the descending order of bins_[i].bin_value;
  std::sort(mapping_data->bins.begin(), mapping_data->bins.end(), bins_greater_than_key());

  double elapsed_secs_conversion = double(clock() - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_region_conversion = elapsed_secs_conversion;
  begin_clock = clock();

//  for (int64_t i=0; i<mapping_data->bins.size(); i++) {
//    printf ("[%ld] ref_id = %ld\tbin_id = %ld\tbin_value = %.2f\n", i, mapping_data->bins[i].reference_id, mapping_data->bins[i].bin_id, mapping_data->bins[i].bin_value);
//  }
//  printf ("\n");
//  fflush(stdout);

  return 0;
}
