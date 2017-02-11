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

void GraphMap::AppendSeedHits_(const uint128_t& seed, std::shared_ptr<is::MinimizerIndex> index, bool threshold_hits, double count_cutoff, std::vector<uint128_t> &all_hits) {
  int64_t key = is::MinimizerIndex::seed_key(seed);
  int32_t pos_read = is::MinimizerIndex::seed_position(seed);

  const uint128_t *found_seeds = NULL;
  int64_t num_found_seeds = 0;

  int lookup_ret = index->KeyLookup(key, &found_seeds, &num_found_seeds);

//  printf ("key = %s = %ld, pos_read = %d, num_found_seeds = %ld\n", is::SeedToString(key, index->get_shape_max_width()).c_str(), key, pos_read, num_found_seeds);

  // Something went wrong.
  if (lookup_ret) {
    return;
  }

  // Filter seeds with excessive hits if necessary.
  if (threshold_hits &&
      num_found_seeds >= count_cutoff) {
    return;
  }

  all_hits.insert(all_hits.end(), num_found_seeds, 0);

  // Reformat the new hits. The key part is no longer needed,
  // but source position on read is.
  auto *phits = &all_hits[0] + (all_hits.size() - num_found_seeds);
  for (int64_t k=0; k<num_found_seeds; k++) {
    if (found_seeds[k] == kInvalidSeed) {
      continue;
//      printf ("Invalid seed found!!\n");
//      fflush(stdout);
//      exit(1);
    }

    uint128_t pos_ref = is::MinimizerIndex::seed_position(found_seeds[k]);
    uint128_t diag = pos_ref - pos_read;
    uint128_t seq_id = is::MinimizerIndex::seed_seq_id(found_seeds[k]);
    phits[k] = (seq_id << 64) | (diag << 32) | (pos_ref);
  }
}

inline int32_t hit_pos_ref(const uint128_t& hit) {
  return (int32_t) (hit & kSeedMask32_1);
}

inline int32_t hit_diag(const uint128_t& hit) {
  return (int32_t) ((hit & kSeedMask32_2) >> 32);
}

inline int32_t hit_seq_id(const uint128_t& hit) {
  return (int32_t) ((hit & kSeedMask32_3) >> 64);
}

std::string PrintSeed(uint128_t seed) {
  std::stringstream ss;

  ss << (uint32_t) ((seed & kSeedMask32_4) >> 96) << " " << (uint32_t) ((seed & kSeedMask32_3) >> 64) << " " << (uint32_t) ((seed & kSeedMask32_2) >> 32) << " " << (uint32_t) ((seed & kSeedMask32_1));

  return ss.str();
}

void FindRegionBounds(const std::vector<uint128_t> &all_hits, int64_t begin_hit, int64_t end_hit, int64_t &start, int64_t &end) {
  if ((end_hit - begin_hit + 1) <= 0) {
    return;
  }

  start = hit_pos_ref(all_hits[begin_hit]);
  end = start;
  for (int64_t i=begin_hit; i<end_hit; i++) {
    int64_t curr_pos = hit_pos_ref(all_hits[i]);
    start = std::min(start, curr_pos);
    end = std::max(end, curr_pos);
  }
}

int GraphMap::RegionSelectionWithSort_(int64_t bin_size, MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes,
                                       const SingleSequence* read, const ProgramParameters* parameters, std::vector<Region>& regions) {
  clock_t begin_clock = clock();
  clock_t diff_clock = begin_clock;

  if (indexes.size() == 0 || indexes[0] == nullptr) {
    FATAL_REPORT(ERR_UNEXPECTED_VALUE, "No reference indexes are specified.");
  }

  int64_t max_shape_width = indexes[0]->get_shape_max_width();
  for (int32_t i=0; i<indexes.size(); i++) {
    max_shape_width = std::max(max_shape_width, indexes[i]->get_shape_max_width());
  }

  int64_t readlength = read->get_sequence_length();
  int64_t num_fwd_seqs = indexes[0]->get_num_sequences_forward();
  bool is_overlapper = (parameters->overlapper == true && parameters->reference_path == parameters->reads_path);
  bool no_self_overlap = (parameters->no_self_hits == true);
  float bin_size_inverse = (bin_size > 0) ? (1.0f / ((float) bin_size)) : (0.0f);

  mapping_data->bin_size = bin_size;


//  exit(1);

  std::vector<uint128_t> all_hits;
  all_hits.reserve(indexes[0]->avg_seed_occurrence() * readlength);

  for (int32_t i=0; i<indexes.size(); i++) {
    std::vector<uint128_t> seeds;       // All seeds for a given index and a given sequence.
    std::vector<int8_t> seed_key_lens;  // Lengths of seed keys. Since multiple lookup keys can be used, each can be of different length.
    std::vector<int8_t> seed_lens;    // Since there are several
    auto index = indexes[i];
    index->CollectLookupSeeds(read->get_data(), read->get_quality(), read->get_sequence_length(), 0.0f, false, false, 1, seeds);
//    index->CollectLookupSeeds(read->get_data(), read->get_quality(), read->get_sequence_length(), 0.0f, false, parameters->use_minimizers, parameters->minimizer_window, seeds);

//    printf ("seeds.size() = %ld\n", seeds.size());
//    fflush(stdout);
    for (int64_t j=0; j<seeds.size(); j++) {
//      printf ("%s\n", PrintSeed(seeds[j]).c_str());
//      fflush(stdout);
      AppendSeedHits_(seeds[j], index, parameters->threshold_hits, index->count_cutoff(), all_hits);
//      printf ("--------------\n");
//      fflush(stdout);
    }
  }

  std::sort(all_hits.begin(), all_hits.end());

  // Cluster minimizer hits.
  int32_t diag_epsilon = 500;     // TODO: Parametrize this.
  int64_t begin_hit = 0;
  int64_t total_num_hits = all_hits.size();
//  printf ("all_hits.size() = %ld\n", all_hits.size());
//  fflush(stdout);

  mapping_data->bins.clear();
  for (int64_t end_hit=0; end_hit<total_num_hits; end_hit++) {
    if ((end_hit + 1) == total_num_hits || hit_seq_id(all_hits[end_hit+1]) != hit_seq_id(all_hits[end_hit]) || (hit_diag(all_hits[end_hit+1]) - hit_diag(all_hits[end_hit])) >= diag_epsilon) {
      // Add a cluster.
      int64_t ref_id = hit_seq_id(all_hits[end_hit]);
      int64_t ref_start = indexes[0]->get_reference_starting_pos()[ref_id];
      int64_t ref_len = indexes[0]->get_reference_lengths()[ref_id];
      int64_t read_len = read->get_sequence_length();

      int64_t start = 0, end = 0;
      FindRegionBounds(all_hits, begin_hit, end_hit, start, end);

      Region region;
      region.start = std::max(ref_start + start - read_len, ref_start);
      region.end = ref_start + std::min(end + max_shape_width + read_len, ref_len);; // TODO: instead of max_seed_len, the exact seed length (for a sepecific lookup key) should be used here.
      region.reference_id = ref_id;
      region.region_index = 0;
      region.region_votes = (end_hit - begin_hit + 1);
      region.is_split = false;
      region.split_start = 0;
      region.split_end = 0;

      begin_hit = end_hit + 1;

      if (region.reference_id == 0xFFFFFFFFFFFFFFFF ||
          region.region_votes < 5) {
        continue;
      }
      regions.emplace_back(region);

//      ChromosomeBin new_bin;
//      new_bin.reference_id = hit_seq_id(all_hits[end_hit]);
//      new_bin.bin_id = j;
//      new_bin.bin_value = bins_chromosome[i][j];
//      if (j > 0) { new_bin.bin_value += bins_chromosome[i][j-1] / 2.0f; }
//      if ((j + 1) < bins_chromosome[i].size()) { new_bin.bin_value += bins_chromosome[i][j+1] / 2.0f; }
//      mapping_data->bins.push_back(new_bin);
    }
  }

  std::sort(regions.begin(), regions.end(), [](const Region& r1, const Region& r2) { return r1.region_votes > r2.region_votes; });

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
//  LOG_DEBUG_SPEC("    total_num_hits = %ld\n", total_num_hits);
//  LOG_DEBUG_SPEC("    read_len = %ld\n", read->get_sequence_length());
//  exit(1);

  return 0;
}

