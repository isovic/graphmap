/*
 * index_gapped_minimizer.cc
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#include "index_gapped_minimizer.h"
#include "log_system/log_system.h"
#include "omp_sort.hpp"
#include <omp.h>
#include <algorithm>
#include <tuple>
#include <deque>

IndexGappedMinimizer::IndexGappedMinimizer() {
  hash_.set_empty_key(empty_hash_key);
  Clear();
}

IndexGappedMinimizer::~IndexGappedMinimizer() {
  Clear();
}

void IndexGappedMinimizer::Clear() {
  count_cutoff_ = 0.0;
  num_sequences_ = 0;
  num_sequences_forward_ = 0;
  seeds_.clear();
  hash_.clear();
  data_.clear();
  reference_lengths_.clear();
  headers_.clear();
}

int IndexGappedMinimizer::CreateFromSequenceFile(const SequenceFile& seqs, const std::vector<CompiledShape>& compiled_shapes, float min_avg_seed_qv, bool index_reverse_strand, bool use_minimizers, int32_t minimizer_window_len, int32_t num_threads) {
  clock_t absolute_time = clock();
  clock_t diff_time = clock();

  Clear();

  /// Determine the maximum length of all shapes, to limit the last kmer being checked.
  int64_t max_seed_len = 0;
  int32_t max_incl_bits = 0;
  for (int32_t i=0; i<compiled_shapes.size(); i++) {
    max_seed_len = std::max(max_seed_len, (int64_t) compiled_shapes[i].shape.size());
    max_incl_bits = std::max(max_incl_bits, compiled_shapes[i].num_incl_bits);
  }

  AssignData_(seqs, index_reverse_strand);

  int64_t total_num_seeds = 0;
  std::vector<int64_t> seed_starts_for_seq;
  AllocateSpaceForSeeds_(seqs, index_reverse_strand, compiled_shapes.size(), max_seed_len, num_sequences_forward_, seed_starts_for_seq, &total_num_seeds);

  LOG_ALL("Allocated memory for a list of %ld seeds (128 bits each) (%.5f sec, diff: %.5f sec).\n", seeds_.size(), (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  diff_time = clock();
  LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
  LOG_ALL("Collecting seeds with %ld threads.\n", num_threads);

  // Process all reads in parallel.
  #pragma omp parallel for num_threads(num_threads) shared(seqs) schedule(dynamic, 1)
  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    uint32_t thread_id = omp_get_thread_num();
    if (thread_id == 0) {
      LOG_ALL("\rRead %ld/%ld", (i + 1), seqs.get_sequences().size());
    }

    int8_t *seqdata = (int8_t *) seqs.get_sequences()[i]->get_data();
    int8_t *seqqual = (int8_t *) seqs.get_sequences()[i]->get_quality();
    int64_t seqlen = seqs.get_sequences()[i]->get_data_length();
//    printf ("\nseqlen = %ld, %s\n", seqlen, GetSubstring((char *) seqdata, 10).c_str());
//    fflush(stdout);

    uint64_t seq_id = seqs.get_sequences()[i]->get_sequence_absolute_id() + 0; // The '+ 0' is actually short for this: + ((orientation == kReverse) ? num_sequences_forward_ : 0);
    int64_t num_seeds_processed = CollectSeedsForSeq_(seqdata, seqqual, seqlen, min_avg_seed_qv, seq_id, compiled_shapes, &(seeds_[seed_starts_for_seq[i]]));
    if (use_minimizers) {
      MakeMinimizers_(&(seeds_[seed_starts_for_seq[i]]), num_seeds_processed, 2, minimizer_window_len);     // 2 only refers to the fwd and rev complement (2 seeds per base).
    }
  }

//  DumpSeeds("temp/seeds.sparse.minimizers.csv", max_incl_bits/2);

  if (use_minimizers) {
    // Remove all excess seeds so that sorting will be faster.
    LOG_ALL("Removing excess seeds.\n");
    int64_t num_dense_seeds = MakeSeedListDense_(&(seeds_[0]), seeds_.size());
    seeds_.resize(num_dense_seeds);
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("\n"), std::string("[]"));

//  DumpSeeds("temp/seeds.dense.minimizers.csv", max_incl_bits/2);

  LOG_ALL("Sorting the seeds (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  diff_time = clock();
  pquickSort(&(seeds_[0]), seeds_.size(), num_threads);
  //    FlagDuplicates_(&(seeds_[seed_starts_for_seq[i]]), num_seeds_processed);
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Making unique (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
//  diff_time = clock();

  LOG_ALL("Generating the hash table (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  diff_time = clock();

//  key_occ_.clear();
//  key_occ_.resize(pow(2, max_incl_bits), 0.0f);
//  for (int64_t i=0; i<seed_list_.size(); i++) {
//    key_occ_[seed_list_[i] >> 32] += 1.0f;
//  }

  /// This part generates the lookup table for each seed key.
//  hash_.resize(pow(2, max_incl_bits), NULL);
//  hash_counts_.resize(pow(2, max_incl_bits), 0);
//  for (int64_t i=0; i<hash_.size(); i++) { hash_[i] = NULL; hash_counts_[i] = 0; };
  uint64_t prev_key = 0;
  bool init = false;
  SeedHashValue new_hash_val;
  for (int64_t i=0; i<seeds_.size(); i++) {
    if (seeds_[i] == 0) { continue; }
    uint64_t key = (uint64_t) (seeds_[i] >> 64);
    if (init == false) {
      new_hash_val.start = i;
      new_hash_val.num = 1;
      init = true;
    } else if (key != prev_key && new_hash_val.num > 0) {
      hash_[prev_key] = new_hash_val;
      new_hash_val.start = i;
      new_hash_val.num = 0;
    } else {
      new_hash_val.num += 1;
    }
    prev_key = key;
  }
  // Store the last hash key in the list.
  if (new_hash_val.num > 0) {
    hash_[prev_key] = new_hash_val;
    new_hash_val.start = 0;
    new_hash_val.num = 0;
  }

  double avg = 0.0f, stddev = 0.0f, perc9999 = 0.0f;
  LOG_ALL("Calculating the distribution statistics for key counts (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  OccurrenceStatistics_(0.9999, num_threads, &avg, &stddev, &perc9999);
  count_cutoff_ = perc9999;
  LOG_ALL("Index statistics: average key count = %f, std dev = %f, percentil(99.99%) = %f\n", avg, stddev, perc9999);

//  DumpSeeds("temp/seeds.csv", max_incl_bits/2);
//  DumpHash("temp/hash.csv", max_incl_bits/2);
//  DumpSortedHash("temp/hash.minimizers.sorted.csv", max_incl_bits/2);

  return 0;
}

void IndexGappedMinimizer::AssignData_(const SequenceFile& seqs, bool index_reverse_strand) {
  /// The array seq_seeds_starts will hold starting positions for seeds comming from
  /// individual reads. I.e. each read will get a part of seq_seeds_ array which it will
  /// modify, and this array marks the part which is designated for this particular read.
  num_sequences_ = (index_reverse_strand == false) ? seqs.get_sequences().size() : seqs.get_sequences().size()*2;
  num_sequences_forward_ = seqs.get_sequences().size();
  reference_lengths_.resize(num_sequences_);
  headers_.resize(num_sequences_);

  data_.resize(num_sequences_);
  for (int64_t i=0; i<num_sequences_forward_; i++) {
    if (seqs.get_sequences()[i]->get_data_format() != kDataFormat2BitSparse) {
      ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Data format is not 2-bit sparse!\n");
    }
    data_[i].assign(seqs.get_sequences()[i]->get_data(), seqs.get_sequences()[i]->get_data() + seqs.get_sequences()[i]->get_data_length());
    headers_[i] = seqs.get_sequences()[i]->get_header();
    reference_lengths_[i] = seqs.get_sequences()[i]->get_sequence_length();
  }

  if (index_reverse_strand == true) {
    for (int64_t i=0; i<num_sequences_forward_; i++) {
      int8_t *revcmp = seqs.get_sequences()[i]->GetReverseComplement();
      if (revcmp == NULL) {
        printf ("revcmp == NULL!\n");
        exit(1);
      }
      data_[num_sequences_forward_ + i].assign(revcmp, revcmp + seqs.get_sequences()[i]->get_data_length());
      if (revcmp) {
        delete[] revcmp;
        revcmp = NULL;
      }

      reference_lengths_[i + num_sequences_forward_] = reference_lengths_[i];
      headers_[i + num_sequences_forward_] = headers_[i];
    }
  }
}

void IndexGappedMinimizer::AllocateSpaceForSeeds_(const SequenceFile& seqs, bool index_reverse_strand, int64_t num_shapes, int64_t max_seed_len, int64_t num_fwd_seqs, std::vector<int64_t>& seed_starts_for_seq, int64_t* total_num_seeds) {
  seeds_.clear();
  seed_starts_for_seq.clear();
  seed_starts_for_seq.resize(seqs.get_sequences().size() * 2);
//  seed_ends_for_seq.clear();
//  seed_ends_for_seq.resize(seqs.get_sequences().size() * 2);
  *total_num_seeds = 0;
  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    seed_starts_for_seq[i] = *total_num_seeds;
    *total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * num_shapes;
    if (index_reverse_strand == true) {
//      seed_starts_for_seq[i + num_sequences_forward_] = *total_num_seeds + (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * num_shapes;
      *total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * num_shapes;
    }
  }

  /// The seed_list_will contain all seeds that were obtained from the input sequences.
  seeds_.resize(*total_num_seeds);
}

int IndexGappedMinimizer::CollectSeedsForSeq_(int8_t* seqdata, int8_t* seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id, const std::vector<CompiledShape>& compiled_shapes, uint128_t* seed_list) {
  int64_t max_seed_len = 0;
  for (int32_t i=0; i<compiled_shapes.size(); i++) {
    max_seed_len = std::max(max_seed_len, (int64_t) compiled_shapes[i].shape.size());
  }

  /// The seqdata will be split in parts separated by N bases (similar to DALIGNER).
  /// Parts shorter than seed_len will be skipped.
  std::vector<int64_t> split_start;
  std::vector<int64_t> split_len;
  int64_t start = 0;
  for (int64_t i=0; i<seqlen; i++) {
    if (seqdata[i] > 3 || (i+1) >= seqlen) {
      int64_t current_len = ((i+1) >= seqlen) ? ((i - start) + 1) : (i - start);
      if (current_len >= max_seed_len) {
        split_start.push_back(start);
        split_len.push_back(current_len);
      }
      start = i + 1;
    }
  }

  int64_t seed_id = 0;

  double qv_sum = 0.0f;
  if (seqqual != NULL) {
    for (int32_t j=0; j<max_seed_len && j<seqlen; j++) {
      qv_sum += (seqqual[j] - 33.0);
    }
  }

  /// Iterate through all split regions.
  for (int64_t i=0; i<split_start.size(); i++) {
    /// We will use a buffer of 2-bit encoded bases, which will keep up to 32 bases at once.
    /// Fill the buffer for gapped spaced seed making. Buffer holds the next 8 bytes of data (max. 32 bases).
    /// Only seed_len bases are used, so the rest of the buffer is just to have the data ready.
    uint64_t buffer = 0;
    int64_t num_bases_ahead = sizeof(buffer) * 8 / 2;

    for (uint64_t pos=0; pos<(split_len[i] - max_seed_len); pos++) {
      /// Prepare the base-buffer. This would be an equivalent of a full-seed, including the
      /// don't care bases. This buffer will be used to extract the gapped spaced seed.
      /// Initialize the buffer.
      if (pos == 0) {
        for (int32_t j=0; j<num_bases_ahead && j<split_len[i]; j++) {
          int8_t seqbase = seqdata[j + split_start[i]];
          buffer |= (((uint64_t) seqbase) << (sizeof(buffer)*8 - j*2 - 2));
        }
        if (seqqual != NULL) {
          qv_sum = 0.0;
          for (int32_t j=0; j<max_seed_len && j<split_len[i]; j++) {
            qv_sum += (seqqual[j] - 33.0);
          }
        }
      } else {  /// Else, shift the buffer and re-fill.
        buffer = buffer << 2;
        if ((pos + num_bases_ahead) < split_len[i]) {
          int8_t seqbase = seqdata[pos + split_start[i] + num_bases_ahead - 1];
          buffer |= (((uint64_t) seqbase) << (0));
        }
        if (seqqual != NULL) {
          qv_sum = qv_sum - seqqual[pos-1] + seqqual[pos+max_seed_len-1];
        }
      }

      if (min_avg_seed_qv > 0 && (qv_sum / max_seed_len) < min_avg_seed_qv) {
        continue;
      }

      /// Extract gapped spaced seeds from the buffer.
      for (int64_t shape_id=0; shape_id<compiled_shapes.size(); shape_id++) {
        uint64_t seed = compiled_shapes[shape_id].CreateSeedFromShape(buffer);
        uint64_t key = SeedHashFunction_(seed);
        uint64_t position = pos + split_start[i] + 1;     // Make the position 1-based.
        seed_list[seed_id] = (((uint128_t) key) << 64) | (((uint128_t) seq_id) << 32) | (((uint128_t) position) << 0);
        seed_id += 1;

        uint64_t rev_seed = ReverseComplementSeed_(seed, compiled_shapes[shape_id].num_incl_bits/2);
        uint64_t rev_key = SeedHashFunction_(rev_seed);
        uint64_t rev_position = position | kIndexIdReverse;
        seed_list[seed_id] = (((uint128_t) rev_key) << 64) | (((uint128_t) seq_id) << 32) | (((uint128_t) rev_position) << 0);
        seed_id += 1;
      }

    }
  }

  return seed_id;
}

inline uint64_t IndexGappedMinimizer::SeedHashFunction_(uint64_t seed) {
  return seed;
}

inline uint64_t IndexGappedMinimizer::ReverseComplementSeed_(uint64_t seed, int32_t num_bases) {
  uint64_t rev_seed = 0;
  uint64_t complement_seed = ((1 << (num_bases * 2)) - 1) - seed;      // Create a complement of the seed.
  for (int32_t i=0; i<num_bases; i++) {
    rev_seed |= (complement_seed & 0x03) << ((num_bases - i - 1) * 2);  // Reverse the bases.
    complement_seed >>= 2;
  }
  return rev_seed;
}

// Parameter window_len specifies the length of the window in the number of bases. For each base, there may be more than one seed
// (e.g. rev. complement, multiple indexes, etc.), and for this reason the parameter num_seeds_per_base is given.
int IndexGappedMinimizer::MakeMinimizers_(uint128_t* seed_list, int64_t num_seeds, int64_t num_seeds_per_base, int32_t window_len) {
  if (window_len > num_seeds) { return 1; }

  std::deque<int>  q(window_len);

  // Setup the initial deque.
  for (int64_t i=0; i<window_len*num_seeds_per_base; i++) {
    // Remove smaller elements if any.
    while ( (!q.empty()) && GET_KEY_FROM_CODED_SEED(seed_list[i]) <= GET_KEY_FROM_CODED_SEED(seed_list[q.back()])) {
      q.pop_back();
    }
    q.push_back(i);
  }

  std::vector<int64_t> minimizer_indices;
  minimizer_indices.reserve(num_seeds);

  // Every other seed is the reverse complement of the previous one.
  // Thus the sliding window will skip 2 instead of 1 seed.
  for (int64_t i=window_len; i<num_seeds; i+=num_seeds_per_base) {
    // Store the largest element of the previous window.
    minimizer_indices.push_back(q.front());
    // Remove the elements which are out of this window-
    while ((!q.empty()) && q.front() <= (i - window_len)) {
      q.pop_front();
    }
    // Remove smaller elements if any.
    while ( (!q.empty()) && GET_KEY_FROM_CODED_SEED(seed_list[i]) <= GET_KEY_FROM_CODED_SEED(seed_list[q.back()])) {
      q.pop_back();
    }

    q.push_back(i);
  }

  minimizer_indices.push_back(q.front());

  // Remove excess seeds.
  std::vector<bool> keep;
  keep.resize(num_seeds, false);
  for (int64_t i=0; i<minimizer_indices.size(); i++) {
    keep[minimizer_indices[i]] = true;
  }
  for (int64_t i=0; i<num_seeds; i++) {
    if (keep[i] == false) {
      seed_list[i] = 0;
    }
  }

  return 0;
}

int IndexGappedMinimizer::FlagDuplicates_(uint128_t* seed_list, int64_t num_seeds) {
  for (int64_t i=1; i<num_seeds; i++) {
    if (seed_list[i-1] == seed_list[i]) {
      seed_list[i-1] = 0;
    }
  }
  return 0;
}

void IndexGappedMinimizer::DumpHash(std::string out_path, int32_t num_bases) {
  FILE *fp = fopen(out_path.c_str(), "w");
  for (auto it=hash_.begin(); it!=hash_.end(); it++) {
    uint64_t key = it->first;
    SeedHashValue shv = it->second;
//    fprintf (fp, "%6X\t%ld\t%ld\n", key, shv.start, shv.num);
    std::string key_string = SeedToString_(key, num_bases);
    fprintf (fp, "%s\t%ld\t%ld\n", key_string.c_str(), shv.start, shv.num);
  }
  fclose(fp);
}

void IndexGappedMinimizer::DumpSortedHash(std::string out_path, int32_t num_bases) {
  std::vector<std::tuple<std::string, int64_t, int64_t> > sorted_hash;
  sorted_hash.reserve(hash_.size());
  for (auto it=hash_.begin(); it!=hash_.end(); it++) {
    uint64_t key = it->first;
    SeedHashValue shv = it->second;
    std::string key_string = SeedToString_(key, num_bases);
    sorted_hash.push_back(std::make_tuple(key_string, shv.start, shv.num));
  }
  std::sort(sorted_hash.begin(), sorted_hash.end(), [](const std::tuple<std::string, int64_t, int64_t>& a, const std::tuple<std::string, int64_t, int64_t>& b) { return ((std::get<2>(a)) > (std::get<2>(b))); });

  FILE *fp = fopen(out_path.c_str(), "w");
  for (int64_t i=0; i<sorted_hash.size(); i++) {
    fprintf (fp, "%s\t%ld\t%ld\n", std::get<0>(sorted_hash[i]).c_str(), std::get<1>(sorted_hash[i]), std::get<2>(sorted_hash[i]));
  }
  fclose(fp);
}

void IndexGappedMinimizer::DumpSeeds(std::string out_path, int32_t num_bases) {
  FILE *fp = fopen(out_path.c_str(), "w");
  for (int64_t i=0; i<seeds_.size(); i++) {
    uint64_t key = (uint64_t) (seeds_[i] >> 64);
    uint64_t coded_pos = seeds_[i] & 0x0000000000000000FFFFFFFFFFFFFFFF;
    IndexPos ipos(coded_pos);
//    fprintf (fp, "%6X %ld\n", key, ipos.get_pos());
    std::string key_string = SeedToString_(key, num_bases);
    fprintf (fp, "%s %ld\n", key_string.c_str(), ipos.get_pos());
  }
  fclose(fp);
}

int IndexGappedMinimizer::OccurrenceStatistics_(double percentil, int32_t num_threads, double* ret_avg, double* ret_stddev, double *ret_percentil_val) {
  if (percentil < 0.0 || percentil > 1.0) { return 1; }

  std::vector<int32_t> key_counts;
  key_counts.resize(hash_.size(), 0);
  int64_t currkey = 0;
  double avg = 0.0, stddev = 0.0, sum = 0.0;
  // Initialize the array of counts. Needed for percentil calculation.
  // Also, calculate the avg on the fly.
  for (auto it = hash_.begin(); it != hash_.end(); it++) {
    key_counts[currkey++] = it->second.num;
    avg += it->second.num;
    sum += it->second.num;
  }
  if (key_counts.size() > 0) { avg /= key_counts.size(); }

  // Calculate the standard deviation.
  for (int64_t i=0; i<key_counts.size(); i++) {
    stddev += (key_counts[i] - avg) * (key_counts[i] - avg);
  }
  if (key_counts.size() > 1) { stddev /= (key_counts.size() - 1); } // Unbiased estimator.
  stddev = sqrt(stddev);

  // Calculate the percentil.
  pquickSort(&(key_counts[0]), key_counts.size(), num_threads);
  double perc_val = key_counts[percentil * (key_counts.size() - 1)];
//  LOG_ALL("key_counts.size() = %ld\n", key_counts.size());
//  LOG_ALL("Simple percentil (percentage of the total array): %f\n", perc_val);

//  double perc_sum = 0.0;
//  for (int64_t i=0; i<key_counts.size(); i++) {
//    if ((perc_sum + key_counts[i]) >= (sum * percentil)) {
//      perc_val = key_counts[i];
//      break;
//    }
//    perc_sum += key_counts[i];
//  }
//  LOG_ALL("Weighted percentil (percentage of the total array): %f\n", perc_val);

  *ret_avg = avg;
  *ret_stddev = stddev;
  *ret_percentil_val = perc_val;

  return 0;
}

int IndexGappedMinimizer::MakeSeedListDense_(uint128_t* seed_list, int64_t num_seeds) {
  int64_t offset = 0;
  int64_t i = 0;
  for (i=0; (i+offset)<num_seeds; i++) {
    while (GET_POS_FROM_CODED_SEED(seed_list[i+offset]) == 0) {
      offset += 1;
    }
    seed_list[i] = seed_list[i+offset];
  }
  return i;
}

void IndexGappedMinimizer::CollectMinimizers(int8_t* seqdata, int8_t* seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id, const std::vector<CompiledShape>& compiled_shapes, int64_t minimizer_window_len, std::vector<uint128_t> &seed_list) {
  int64_t num_seeds_processed = CollectSeedsForSeq_(seqdata, seqqual, seqlen, min_avg_seed_qv, seq_id, compiled_shapes, &(seed_list[0]));
  MakeMinimizers_(&(seed_list[0]), num_seeds_processed, 2*compiled_shapes.size(), minimizer_window_len);     // 2 only refers to the fwd and rev complement (2 seeds per base).
  int64_t num_dense_seeds = MakeSeedListDense_(&(seed_list[0]), seed_list.size());
  seed_list.resize(num_dense_seeds);
}

inline std::string IndexGappedMinimizer::SeedToString_(uint64_t seed,
                                                       int32_t num_bases) {
  std::stringstream ss;
  for (int32_t i=0; i<num_bases; i++) {
    ss << "ACGT"[(seed & 0x03)];
    seed >>= 2;
  }
  std::string seed_string = ss.str();
  std::reverse(seed_string.begin(), seed_string.end());
  return seed_string;
}
