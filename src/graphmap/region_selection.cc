/*
 * region_selection.cc
 *
 *  Created on: Mar 19, 2015
 *      Author: isovic
 */

#include <string.h>
#include <algorithm>
#include "graphmap/graphmap.h"
#include "log_system/log_system.h"
#include "sparsehash/dense_hash_map"

using google::dense_hash_map;      // namespace where class lives by default

struct DenseType2 {
  int32_t timestamp = 0;
  float count = 0.0;
};

typedef unsigned __int128 uint128_t;

int GraphMap::RegionSelectionNoCopy_(int64_t bin_size, MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
  clock_t begin_clock = clock();
  clock_t diff_clock = begin_clock;

  if (indexes.size() == 0 || indexes[0] == nullptr) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "No reference indexes are specified."));
  }

  int64_t readlength = read->get_sequence_length();
  int64_t num_fwd_seqs = indexes[0]->get_num_sequences_forward();
  bool is_overlapper = (parameters->overlapper == true && parameters->reference_path == parameters->reads_path);
  bool no_self_overlap = (parameters->no_self_hits == true);

  mapping_data->bin_size = bin_size;

  float bin_size_inverse = (bin_size > 0) ? (1.0f / ((float) bin_size)) : (0.0f);

  ////////////////////////////////////////////////////
  ///// This part prepares the bins. /////
  ////////////////////////////////////////////////////
  // Create bins for each chromosome (or reference sequence) separately.
  diff_clock = clock();
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

  int64_t k = (int64_t) indexes[0]->get_shape_max_width(); // get_shape_index_length();

  mapping_data->time_region_alloc = ((double) clock() - diff_clock) / CLOCKS_PER_SEC;
  diff_clock = clock();

  ////////////////////////////////////////////////////
  ///// This part counts the occurrences in bins. /////
  ////////////////////////////////////////////////////
  // Filling the bins with values, so we get an occurrence map.
  float max_bin_value = -1.0f;
  mapping_data->time_region_seed_lookup = 0.0;
  int64_t total_num_hits = 0;
  diff_clock = clock();
  for (int64_t i = 0; i < (readlength - k + 1); i += parameters->kmer_step) {  // i++) {
    int8_t *seed = (int8_t *) &(read->get_data()[i]);

    for (int64_t index_id = 0; index_id < indexes.size(); index_id++) {
      auto& index = indexes[index_id];

      if (index != NULL) {
        clock_t diff_find_seeds = clock();
        std::vector<const uint128_t *> hit_vector;
        std::vector<int64_t> hit_counts;
        int ret_search = index->Find(seed, k, parameters->threshold_hits, hit_vector, hit_counts);
        mapping_data->time_region_seed_lookup += ((double) clock() - diff_find_seeds) / CLOCKS_PER_SEC;

        // Check if there is too many hits (or too few).
        if (ret_search == 1) {
          mapping_data->num_seeds_with_no_hits += 1;
        } else if (ret_search == 2) {
          mapping_data->num_seeds_over_limit += 1;
          continue;
        } else if (ret_search > 2) {
          mapping_data->num_seeds_errors += 1;
        }

        // Counting kmers in regions of bin_size on the genome
//        printf ("[%ld[ num_hits = %ld\n", i, num_hits);

        for (int64_t hits_id = 0; hits_id < hit_vector.size(); hits_id++) {
          auto hits = hit_vector[hits_id];
          total_num_hits += hit_counts[hits_id];

          for (int64_t j = 0; j < hit_counts[hits_id]; j++) {
//            int64_t position = hits[j];
            int64_t local_position = is::MinimizerIndex::seed_position(hits[j]);
            int64_t reference_index = is::MinimizerIndex::seed_seq_id(hits[j]);
//            int64_t local_position = (int64_t) (((uint64_t) position) & MASK_32_BIT);
//            int64_t reference_index = (int64_t) (((uint64_t) position) >> 32);  // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);

            if ((is_overlapper == true && (reference_index % num_fwd_seqs) == read->get_sequence_id()) ||
                (no_self_overlap == true && index->get_headers()[reference_index % num_fwd_seqs] == std::string(read->get_header()))) {
              continue;
            }

            if (reference_index < 0) {
              LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. reference_index = %ld, y = %ld, j = %ld / (%ld, %ld)\n", reference_index, local_position, j, 0, hit_counts[hits_id]), "SelectRegionsWithHoughAndCircular");
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
//            printf ("\nreference_index = %ld, position_bin = %ld, hits_id = %ld/%ld, j = %ld/%ld\n", reference_index, position_bin, hits_id, hit_vector.size(), j, hit_counts[hits_id]);
//            fflush(stdout);

            if (last_update_chromosome[reference_index][position_bin] == (i + 1)) {
              continue;
            }
            if (reference_index >= bins_chromosome.size() ||
                position_bin >= bins_chromosome[reference_index].size()) {
              continue;
            }

            bins_chromosome[reference_index][position_bin] += 1.0f;
            if (bins_chromosome[reference_index][position_bin] > max_bin_value) { max_bin_value = bins_chromosome[reference_index][position_bin]; }

            last_update_chromosome[reference_index][position_bin] = (i + 1);
          }  // for (int64_t j=hits_start; j<(hits_start + num_hits); j++)
        }

      }
    }
  }  // for (int64_t i=0; i<(readlength - parameters->k_region + 1); i++)

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL_DEBUG, read->get_sequence_id() == parameters->debug_read, FormatString("\n[BuildOccuranceMap] k_region = %d, num_seeds_with_no_hits = %ld, num_seeds_over_limit = %ld\n", parameters->k_region, mapping_data->num_seeds_with_no_hits, mapping_data->num_seeds_over_limit), "ProcessKmersInBins_");
//  LOG_DEBUG_HIGH("total_num_hits = %ld\n", total_num_hits);

  mapping_data->time_region_counting = ((double) clock() - diff_clock) / CLOCKS_PER_SEC;
  diff_clock = clock();

  float min_allowed_bin_value = std::max(2.0f, (float) std::floor(parameters->min_bin_percent * max_bin_value));
  int64_t num_bins_above_min = 0;
  for (int64_t i = 0; i < (indexes[0]->get_num_sequences_forward() * 2); i++) {
    for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
      if (bins_chromosome[i][j] > min_allowed_bin_value) { num_bins_above_min += 1; }
    }
  }

  // Convert the bins to a more compact form, which will be easier to sort.
  // The tuple will contain: reference_id, bin_index, bin_count.
  mapping_data->bins.clear();
  mapping_data->bins.resize(num_bins_above_min);
  for (int64_t i = 0; i < (indexes[0]->get_num_sequences_forward() * 2); i++) {
    for (int64_t j = 0; j < bins_chromosome[i].size(); j++) {
      if (bins_chromosome[i][j] > min_allowed_bin_value) {
        ChromosomeBin new_bin;
        new_bin.reference_id = i;
        new_bin.bin_id = j;
        new_bin.bin_value = bins_chromosome[i][j];
        if (j > 0) { new_bin.bin_value += bins_chromosome[i][j-1] / 2.0f; }
        if ((j + 1) < bins_chromosome[i].size()) { new_bin.bin_value += bins_chromosome[i][j+1] / 2.0f; }
        mapping_data->bins.push_back(new_bin);
      }
    }
  }

  mapping_data->time_region_conversion = ((double) clock() - diff_clock) / CLOCKS_PER_SEC;
  diff_clock = clock();

  // Sort the bins in the descending order of bins_[i].bin_value;
  std::sort(mapping_data->bins.begin(), mapping_data->bins.end(), bins_greater_than_key());

  mapping_data->time_region_hitsort = ((double) clock() - diff_clock) / CLOCKS_PER_SEC;
  diff_clock = clock();

  clock_t end_clock = clock();
  double elapsed_secs = double(end_clock - begin_clock) / CLOCKS_PER_SEC;
  mapping_data->time_region_selection = elapsed_secs;
  LOG_DEBUG_SPEC("Region selection timings:\n");
  LOG_DEBUG_SPEC("    time_region_seed_lookup = %f\n", mapping_data->time_region_seed_lookup);
  LOG_DEBUG_SPEC("    time_region_alloc = %f\n", mapping_data->time_region_alloc);
  LOG_DEBUG_SPEC("    time_region_counting = %f\n", mapping_data->time_region_counting);
  LOG_DEBUG_SPEC("    time_region_conversion = %f\n", mapping_data->time_region_conversion);
  LOG_DEBUG_SPEC("    time_region_sort = %f\n", mapping_data->time_region_hitsort);
  LOG_DEBUG_SPEC("\n");
  LOG_DEBUG_SPEC("    total_num_hits = %ld\n", total_num_hits);
  LOG_DEBUG_SPEC("    read_len = %ld\n", read->get_sequence_length());
//  exit(1);

  return 0;
}
