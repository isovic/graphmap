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

IndexGappedMinimizer::IndexGappedMinimizer() {
  // TODO Auto-generated constructor stub

}

IndexGappedMinimizer::~IndexGappedMinimizer() {
  // TODO Auto-generated destructor stub
}

void IndexGappedMinimizer::Clear() {
  seeds_.clear();
  hash_.clear();
}

int IndexGappedMinimizer::CreateFromSequenceFile(const SequenceFile& seqs, const std::vector<CompiledShape>& compiled_shapes, float min_avg_seed_qv, bool index_reverse_strand, bool use_minimizers, int32_t num_threads) {
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

  /// The seed_list_will contain all seeds that were obtained from the input sequences.
  seeds_.resize(total_num_seeds);

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
    int64_t num_seeds_processed = AddSeedsForSeq_(seqdata, seqqual, seqlen, min_avg_seed_qv, seq_id, compiled_shapes, &(seeds_[seed_starts_for_seq[i]]));
    if (use_minimizers) {
      MakeMinimizers_(&(seeds_[seed_starts_for_seq[i]]), num_seeds_processed);
    }
    FlagDuplicates_(&(seeds_[seed_starts_for_seq[i]]), num_seeds_processed);
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("\n"), std::string("[]"));

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Sorting the seeds (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();
  pquickSort(&(seeds_[0]), seeds_.size(), num_threads);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Making unique (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Generating occurence counts (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();

//  key_occ_.clear();
//  key_occ_.resize(pow(2, max_incl_bits), 0.0f);
//  for (int64_t i=0; i<seed_list_.size(); i++) {
//    key_occ_[seed_list_[i] >> 32] += 1.0f;
//  }

//  //   For debugging purposes, output the occurrence table to a file.
//    auto temp = key_occ_;
//    std::sort(temp.begin(), temp.end(), [](float a, float b) { return b < a; });
//    FILE *fp = fopen("dist.csv", "w");
//    for (int64_t i=0; i<temp.size(); i++) {
//      fprintf (fp, "%ld\n", (int64_t) temp[i]);
//    }
//    fclose(fp);

//  for (int64_t i=0;i <seed_list_.size(); i++) {
//    int32_t key = (int32_t) (seed_list_[i] >> ((uint64_t) 32));
//    int32_t ref = seed_list_[i] & 0x00000000FFFFFFFF;
//    printf ("seed[%ld] = %6X %d\n", i, key, ref);
//  }

  /// This part generates the lookup table for each seed key.
//  hash_.resize(pow(2, max_incl_bits), NULL);
//  hash_counts_.resize(pow(2, max_incl_bits), 0);
//  for (int64_t i=0; i<hash_.size(); i++) { hash_[i] = NULL; hash_counts_[i] = 0; };
  hash_.set_empty_key(0);
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
      new_hash_val.start = 0;
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

//  // Sort the counts so that we can get an occurance histogram.
//  // We will select a 99.99 percentile of the data as the cutoff value.
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Sorting the hash_couts_ to obtain the occurrence histogram (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
//  diff_time = clock();
//  std::vector<size_t> hash_counts_indices;
//  ordered_sort_vector(hash_counts_, hash_counts_indices);
//  int64_t num_nonzero = 0;
//  for (int64_t i=0; i<hash_counts_.size(); i++) {
//    if (hash_counts_[hash_counts_indices[i]] != 0) {
//      num_nonzero = hash_counts_.size() - i;
//      break;
//    }
//  }
//
//  // Determine the index of the 99.99 percentile.
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Calculating the 99.99 percentile (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
//  diff_time = clock();
//  int64_t id99 = hash_counts_.size() - (int64_t) round(((double) 1.0 - 0.9999) * ((double) num_nonzero));
//  int64_t val99 = hash_counts_[hash_counts_indices[id99]];
//
//  count_cutoff_ = val99;
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("num_nonzero = %ld, id99 = %ld (%.5f sec, diff: %.5f sec).\n", num_nonzero, id99, (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
//
//  for (int64_t i=0; i<hash_counts_.size(); i++) {
//    if (hash_counts_[i] >= count_cutoff_) {
//      hash_counts_[i] = 0;
//    }
//  }
//
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("The 99.99 percentile is %ld (%.5f sec, diff: %.5f sec).\n", val99, (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));

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
}

int IndexGappedMinimizer::AddSeedsForSeq_(int8_t* seqdata, int8_t* seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id, const std::vector<CompiledShape>& compiled_shapes, uint128_t* seed_list) {
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

uint64_t IndexGappedMinimizer::SeedHashFunction_(uint64_t seed) {
  return seed;
}

uint64_t IndexGappedMinimizer::ReverseComplementSeed_(uint64_t seed, int32_t num_bases) {
  uint64_t rev_seed = 0;
  uint64_t complement_seed = ((1 << (num_bases * 2)) - 1) - seed;      // Create a complement of the seed.
  for (int32_t i=0; i<num_bases; i++) {
    rev_seed |= (complement_seed & 0x03) << ((num_bases - i - 1) * 2);  // Reverse the bases.
  }
  return rev_seed;
}

int IndexGappedMinimizer::MakeMinimizers_(uint128_t* seed_list, int64_t num_seeds) {
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
