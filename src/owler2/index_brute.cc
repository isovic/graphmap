/*
 * index_brute.cc
 *
 *  Created on: Mar 3, 2016
 *      Author: isovic
 */

#include "index_brute.h"
#include "log_system/log_system.h"
#include "qsort.h"
#include "omp_sort.hpp"
#include <omp.h>

IndexBrute::IndexBrute() {
  Clear();
}

IndexBrute::IndexBrute(std::string seq_path, std::vector<std::string> &shapes, float min_avg_seed_qv, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads) {
  Create(seq_path, shapes, min_avg_seed_qv, index_reverse_strand, index_desig, batch_size, num_threads);
}

IndexBrute::~IndexBrute() {
  Clear();
}

void IndexBrute::Clear() {
  seed_list_.clear();
  hash_.clear();
  hash_counts_.clear();
  num_sequences_ = 0;
  num_sequences_forward_ = 0;
  reference_lengths_.clear();
  headers_.clear();
  count_cutoff_ = 0;
  data_.clear();
}

int IndexBrute::Create(std::string seq_path, std::vector<std::string> &shapes, float min_avg_seed_qv, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads) {
  SequenceFile seqs;
  seqs.OpenFileForBatchLoading(seq_path);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Started creating the index.\n"), std::string(__FUNCTION__));

  std::vector<CompiledShape> compiled_shapes = CompileShapes(shapes);
  int64_t max_seed_len = 0;
  for (int32_t i=0; i<compiled_shapes.size(); i++) {
    max_seed_len = std::max(max_seed_len, (int64_t) compiled_shapes[i].shape.size());
  }

  clock_t absolute_time = clock();
  clock_t last_batch_loading_time = clock();

  /// Load sequences in batch (if requested), or all at once.
  while ((batch_size <= 0 && !seqs.LoadAllAsBatch(SEQ_FORMAT_AUTO, true, false)) || (batch_size > 0 && !seqs.LoadNextBatchInMegabytes(SEQ_FORMAT_AUTO, batch_size, true, false))) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Batch of %ld reads (%ld MiB) loaded in %.2f sec. (%ld bases)\n", seqs.get_sequences().size(), seqs.CalculateTotalSize(MEMORY_UNIT_MEGABYTE), (((float) (clock() - last_batch_loading_time))/CLOCKS_PER_SEC), seqs.GetNumberOfBases()), "ProcessReads");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), std::string(__FUNCTION__));

    clock_t time_before_processing = clock();

    CreateFromSequenceFile(seqs, min_avg_seed_qv, compiled_shapes, index_reverse_strand, index_desig, num_threads);

    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Processed a batch of %ld reads (%ld MiB) in %.2f sec.\n", seqs.get_sequences().size(), seqs.CalculateTotalSize(MEMORY_UNIT_MEGABYTE), (((float) (clock() - time_before_processing))/CLOCKS_PER_SEC), seqs.GetNumberOfBases()), std::string(__FUNCTION__));
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), std::string(__FUNCTION__));
    if (batch_size > 0) {  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("\n"), "[]");  }
    last_batch_loading_time = clock();
  }

  seqs.CloseFileAfterBatchLoading();

  return 0;
}

const std::vector<float>& IndexBrute::get_key_occ() const {
  return key_occ_;
}

void IndexBrute::set_key_occ(const std::vector<float>& keyOcc) {
  key_occ_ = keyOcc;
}

const std::vector<int64_t>& IndexBrute::get_seq_seed_starts() const {
  return seq_seed_starts_;
}

void IndexBrute::set_seq_seed_starts(const std::vector<int64_t>& seqSeedStarts) {
  seq_seed_starts_ = seqSeedStarts;
}

int IndexBrute::ProcessKmerSpectrum_(int8_t *seqdata, int8_t *seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id, const std::vector<CompiledShape> &compiled_shapes, uint64_t *seed_list) {
  int64_t max_seed_len = 0;
  for (int32_t i=0; i<compiled_shapes.size(); i++) {
    max_seed_len = std::max(max_seed_len, (int64_t) compiled_shapes[i].shape.size());
  }

//  if (seq->get_data_format() != kDataFormat2BitSparse) return 1;

//  int8_t *seqdata = (int8_t *) seq->get_data();
//  int64_t seqlen = seq->get_data_length();
//  uint64_t seq_id = seq->get_sequence_absolute_id() + ((orientation == kReverse) ? num_forward_seqs : 0);

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
        int64_t seed = CreateSeedFromShape(compiled_shapes[shape_id], buffer);
        int32_t key = (int32_t) seed;
        int32_t position = pos + split_start[i];
        seed_list[seed_id] = ((uint64_t) key << 32) | (seq_id);

//        auto it = seeds_to_reads_.find(key);
//        if (it == seeds_to_reads_.end()) {
//          std::vector<int32_t> new_vec = {seq_id};
//          seeds_to_reads_[key] = new_vec;
//        } else {
//          it->second.push_back(seq_id);
//        }

//        seed_list[seed_id].seq_id = seq->get_sequence_absolute_id();
//        seed_list[seed_id].orientation = orientation;
        seed_id += 1;
      }

    }
  }

  return seed_id;
}

const std::vector<int64_t>& IndexBrute::get_seq_seed_counts() const {
  return seq_seed_counts_;
}

void IndexBrute::set_seq_seed_counts(const std::vector<int64_t>& seqSeedCounts) {
  seq_seed_counts_ = seqSeedCounts;
}

int64_t IndexBrute::get_count_cutoff() const {
  return count_cutoff_;
}

void IndexBrute::set_count_cutoff(int64_t countCutoff) {
  count_cutoff_ = countCutoff;
}

const std::vector<std::vector<int8_t> >& IndexBrute::get_data() const {
  return data_;
}

void IndexBrute::set_data(const std::vector<std::vector<int8_t> >& data) {
  data_ = data;
}

int IndexBrute::ProcessKmerSpectrumSeparateReads_(int8_t *seqdata, int8_t *seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id, const std::vector<CompiledShape> &compiled_shapes, uint64_t *seed_list) {
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
        int64_t seed = CreateSeedFromShape(compiled_shapes[shape_id], buffer);
        int32_t key = (int32_t) seed;
        int32_t position = pos + split_start[i];
        seed_list[seed_id] = ((uint64_t) key << 32) | (position);
        seed_id += 1;
      }
    }
  }

  std::sort(seed_list, seed_list + seed_id);

  return seed_id;
}

std::vector<CompiledShape> IndexBrute::CompileShapes(const std::vector<std::string> &shapes) {
  std::vector<CompiledShape> compiled_shapes;
  for (int32_t i=0; i<shapes.size(); i++) {
    CompiledShape compiled_shape(shapes[i]);
    compiled_shapes.push_back(compiled_shape);
  }
  return compiled_shapes;
}

const std::vector<uint64_t *>& IndexBrute::get_hash() const {
  return hash_;
}

void IndexBrute::set_hash(const std::vector<uint64_t *>& hash) {
  hash_ = hash;
}

const std::vector<int64_t>& IndexBrute::get_hash_counts() const {
  return hash_counts_;
}

void IndexBrute::set_hash_counts(const std::vector<int64_t>& hashCounts) {
  hash_counts_ = hashCounts;
}

const std::vector<uint64_t>& IndexBrute::get_seed_list() const {
  return seed_list_;
}

void IndexBrute::set_seed_list(const std::vector<uint64_t>& seedList) {
  seed_list_ = seedList;
}

//int IndexBrute::CreateFromSequence(SingleSequence& seq, std::vector<CompiledShape> &compiled_shapes, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads) {
//  Clear();
//
//  /// Determine the maximum length of all shapes, to limit the last kmer being checked.
//  int64_t max_seed_len = 0;
//  int32_t max_incl_bits = 0;
//  for (int32_t i=0; i<compiled_shapes.size(); i++) {
//    max_seed_len = std::max(max_seed_len, (int64_t) compiled_shapes[i].shape.size());
//    max_incl_bits = std::max(max_incl_bits, compiled_shapes[i].num_incl_bits);
//  }
//
//  /// The array seq_seeds_starts will hold starting positions for seeds comming from
//  /// individual reads. I.e. each read will get a part of seq_seeds_ array which it will
//  /// modify, and this array marks the part which is designated for this particular read.
//  int64_t total_num_seeds = 0;
//  if (index_reverse_strand == true) {
//    total_num_seeds += (seq.get_sequence_length() - max_seed_len) * 2 * compiled_shapes.size();
//  } else {
//    total_num_seeds += (seq.get_sequence_length() - max_seed_len) * compiled_shapes.size();
//  }
//
//  /// The seed_list_will contain all seeds that were obtained from the input sequences.
//  seed_list_.resize(total_num_seeds);
//
////  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("\n"), std::string(__FUNCTION__));
////  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Allocated memory for a list of %ld seeds (128 bits each).\n", seed_list_.size()), std::string(__FUNCTION__));
////  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), std::string(__FUNCTION__));
//
//  // Process all reads in parallel.
//    SingleSequence *seq2bit = seq.ConvertDataFormatAndReturn(kDataFormat2BitSparse);
////    if (seq2bit == NULL) continue;
//    int64_t num_seeds_processed = ProcessKmerSpectrum_(seq2bit, compiled_shapes, kForward, &(seed_list_[0]));
//    if (seq2bit) delete seq2bit;
//    seq2bit = NULL;
//
//    if (index_reverse_strand == true) {
//      SingleSequence seq_reverse;
//      seq_reverse.CopyFrom(seq);
//      seq_reverse.ReverseComplement();
//      SingleSequence *seq2bit_rev = seq_reverse.ConvertDataFormatAndReturn(kDataFormat2BitSparse);
////      if (seq2bit_rev == NULL) continue;
//      int64_t num_seeds_processed_rev = ProcessKmerSpectrum_(seq2bit_rev, compiled_shapes, kReverse, &(seed_list_[0 + num_seeds_processed]));
//      if (seq2bit_rev) delete seq2bit_rev;
//      seq2bit_rev = NULL;
//    }
//
////  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Sorting seeds.\n"), std::string(__FUNCTION__));
//
//  /// After sorting, there will be seeds which are equal to 0 (position, seq_id and key). These are caused by non ACTG bases
//  /// because the space is preallocated to the total number of bases in the sequence. If any of those bases were 'N', some
//  /// seed positions would not be filled, and therefore remained empty.
//  /// 1) Standard C implementation, 3-4x faster than std::sort.
////   qsort(&seed_list_[0], seed_list_.size(), sizeof(HashLocation), compare_seed_hits);
//  /// 2) Slow, it took about ~430 CPU sec to process ~238 million 128-bit elements.
//   std::sort(seed_list_.begin(), seed_list_.end());
//  /// 3) This is a very fast serial implementation (inline). Faster than C's qsort (comparison on ~238 elements: 89CPUsec vs. 72CPUsec).
////   QSORT(HashLocation, &(seed_list_[0]), seed_list_.size(), islt);
//  /// 4) OpenMP version of the quicksort algorithm. Could be faster (need to find implementations which do not use STL templates).
////   pquickSort(&(seed_list_[0]), seed_list_.size());
//
//  /// This part generates the lookup table for each seed key.
////  hash_.resize(pow(2, max_incl_bits));
////  hash_counts_.resize(pow(2, max_incl_bits));
////  for (int64_t i=0; i<hash_.size(); i++) { hash_[i] = -1; hash_counts_[i] = 0; };
//  int64_t num_keys = 0;
//  int32_t prev_key = 0;
//  for (int64_t i=0; i<seed_list_.size(); i++) {
//    int32_t key = seed_list_[i].key;
//    if (i == 0 || (i > 0 && key != prev_key)) {
//      num_keys += 1;
//    }
//    prev_key = key;
//  }
//
//  hash_.resize(num_keys);
//  hash_counts_.resize(num_keys);
////  printf ("num_keys = %ld\n", num_keys);
////  fflush(stdout);
//
//  int64_t count = 0;
//  num_keys = 0;
//  prev_key = 0;
//  for (int64_t i=0; i<seed_list_.size(); i++) {
//    int32_t key = seed_list_[i].key;
//    if (i == 0 || (i > 0 && key != prev_key)) {
//      hash_[num_keys] = i;
//      if (num_keys > 0) {
//        hash_counts_[num_keys-1] = i - hash_[num_keys-1];
//      }
//      num_keys += 1;
//
//      count = 0;
//    }
//    prev_key = key;
//  }
//  if (seed_list_.size() > 0) {
//    hash_counts_.back() = seed_list_.size() - hash_[num_keys-1];
//  }
//
////  printf ("Tu sam 1!\n");
////  fflush(stdout);
////
////  for (int64_t i=0; i<hash_.size(); i++) {
////    printf ("i = %ld\thash_[i] = %d\thash_counts_[i] = %d\n", i, hash_[i], hash_counts_[i]);
////    fflush(stdout);
////  }
////  exit(1);
//
//
//
////  /// Fix the edge condition (the last seed key).
////  if (seed_list_.size() > 0) {
////    int32_t prev_key = seed_list_.back().key;
////    hash_counts_[prev_key] = seed_list_.size() - hash_[prev_key];
////  }
//
//  return 0;
//}

uint64_t IndexBrute::CreateSeedFromShape(const CompiledShape &compiled_shape, uint64_t bases2bit) {
  uint64_t seed = 0;

  for (int32_t i=0; i<compiled_shape.masks.size(); i++) {
    uint64_t buffer_part = bases2bit & compiled_shape.masks[i].bits;
    buffer_part = buffer_part << compiled_shape.masks[i].shift;
    seed |= buffer_part;
  }

  seed = seed >> (64 - compiled_shape.num_incl_bits);

  return seed;
}


int IndexBrute::CreateFromSequenceFile(const SequenceFile &seqs, float min_avg_seed_qv, const std::vector<CompiledShape> &compiled_shapes, bool index_reverse_strand, std::string index_desig, int32_t num_threads) {
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

  /// The array seq_seeds_starts will hold starting positions for seeds comming from
  /// individual reads. I.e. each read will get a part of seq_seeds_ array which it will
  /// modify, and this array marks the part which is designated for this particular read.
  num_sequences_ = (index_reverse_strand == false) ? seqs.get_sequences().size() : seqs.get_sequences().size()*2;
  num_sequences_forward_ = seqs.get_sequences().size();
  reference_lengths_.resize(num_sequences_);
  headers_.resize(num_sequences_);

  data_.resize(num_sequences_);
  for (int64_t i=0; i<num_sequences_forward_; i++) {
    data_[i].assign(seqs.get_sequences()[i]->get_data(), seqs.get_sequences()[i]->get_data() + seqs.get_sequences()[i]->get_data_length());
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
    }
  }

  int64_t total_num_seeds = 0;
//  std::vector<int64_t> seq_seed_starts;
  seq_seed_starts_.clear();
  seq_seed_starts_.resize(seqs.get_sequences().size() * 2);
  seq_seed_counts_.clear();
  seq_seed_counts_.resize(seqs.get_sequences().size() * 2);
  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    if (seqs.get_sequences()[i]->get_data_format() != kDataFormat2BitSparse) {
      ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Data format is not 2-bit sparse!\n");
    }

    reference_lengths_[i] = seqs.get_sequences()[i]->get_sequence_length();
    headers_[i] = seqs.get_sequences()[i]->get_header();
    if (index_reverse_strand == true) {
      reference_lengths_[i + num_sequences_forward_] = reference_lengths_[i];
      headers_[i + num_sequences_forward_] = headers_[i];
    }

    seq_seed_starts_[i] = total_num_seeds;
    if (index_reverse_strand == true) {
      seq_seed_starts_[i + num_sequences_forward_] = total_num_seeds + (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * compiled_shapes.size();
      total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * 2 * compiled_shapes.size();
    } else {
      total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * compiled_shapes.size();
    }
  }

  /// The seed_list_will contain all seeds that were obtained from the input sequences.
  seed_list_.resize(total_num_seeds);

  LOG_ALL("Allocated memory for a list of %ld seeds (128 bits each) (%.5f sec, diff: %.5f sec).\n", seed_list_.size(), (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
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
    int64_t num_seeds_processed = ProcessKmerSpectrum_(seqdata, seqqual, seqlen, min_avg_seed_qv, seq_id, compiled_shapes, &(seed_list_[seq_seed_starts_[i]]));

    if (index_reverse_strand == true) {
      int8_t *seqrevdata = (int8_t *) seqs.get_sequences()[i]->GetReverseComplement();
      int8_t *seqrevqual = reverse_data(seqqual, seqlen);
      int64_t seqlen = seqs.get_sequences()[i]->get_data_length();
      uint64_t seq_id = seqs.get_sequences()[i]->get_sequence_absolute_id() + num_sequences_forward_; // The '+ num_sequences_forward_' is actually short for this: + ((orientation == kReverse) ? num_sequences_forward_ : 0);

      if (seqrevdata == NULL) continue;
      int64_t num_seeds_processed_rev = ProcessKmerSpectrum_(seqrevdata, seqrevqual, seqlen, min_avg_seed_qv, seq_id, compiled_shapes, &(seed_list_[seq_seed_starts_[i + num_sequences_forward_]]));
      if (seqrevdata) delete seqrevdata;
      seqrevdata = NULL;
    }
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("\n"), std::string("[]"));

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Sorting the seeds (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();
//  std::sort(seed_list_.begin(), seed_list_.end() );
  pquickSort(&(seed_list_[0]), seed_list_.size(), num_threads);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Making unique (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();
  seed_list_.erase(std::unique(seed_list_.begin(), seed_list_.end()), seed_list_.end());

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true,
                               FormatString("Generating occurence counts (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();

  key_occ_.clear();
  key_occ_.resize(pow(2, max_incl_bits), 0.0f);
  for (int64_t i=0; i<seed_list_.size(); i++) {
    key_occ_[seed_list_[i] >> 32] += 1.0f;
  }

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
  hash_.resize(pow(2, max_incl_bits), NULL);
  hash_counts_.resize(pow(2, max_incl_bits), 0);
//  for (int64_t i=0; i<hash_.size(); i++) { hash_[i] = NULL; hash_counts_[i] = 0; };
  int32_t prev_key = 0;
  for (int64_t i=0; i<seed_list_.size(); i++) {
    int32_t key = (int32_t) (seed_list_[i] >> ((uint64_t) 32));
    if (hash_[key] == NULL) { hash_[key]= &(seed_list_[i]); hash_counts_[key] = 1; }
    else { hash_counts_[key] += 1; }
  }

  // Sort the counts so that we can get an occurance histogram.
  // We will select a 99.99 percentile of the data as the cutoff value.
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Sorting the hash_couts_ to obtain the occurrence histogram (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();
  std::vector<size_t> hash_counts_indices;
  ordered_sort_vector(hash_counts_, hash_counts_indices);
  int64_t num_nonzero = 0;
  for (int64_t i=0; i<hash_counts_.size(); i++) {
    if (hash_counts_[hash_counts_indices[i]] != 0) {
      num_nonzero = hash_counts_.size() - i;
      break;
    }
  }

  // Determine the index of the 99.99 percentile.
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Calculating the 99.99 percentile (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();
  int64_t id99 = hash_counts_.size() - (int64_t) round(((double) 1.0 - 0.9999) * ((double) num_nonzero));
  int64_t val99 = hash_counts_[hash_counts_indices[id99]];

  count_cutoff_ = val99;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("num_nonzero = %ld, id99 = %ld (%.5f sec, diff: %.5f sec).\n", num_nonzero, id99, (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));

  for (int64_t i=0; i<hash_counts_.size(); i++) {
    if (hash_counts_[i] >= count_cutoff_) {
      hash_counts_[i] = 0;
    }
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("The 99.99 percentile is %ld (%.5f sec, diff: %.5f sec).\n", val99, (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));

  return 0;
}

int IndexBrute::CreateFromSequenceFileSeparateReads(const SequenceFile &seqs, float min_avg_seed_qv, const std::vector<CompiledShape> &compiled_shapes, bool index_reverse_strand, std::string index_desig, int64_t batch_size, int32_t num_threads) {
  Clear();

  clock_t absolute_time = clock();
  clock_t diff_time = clock();

  /// Determine the maximum length of all shapes, to limit the last kmer being checked.
  int64_t max_seed_len = 0;
  int32_t max_incl_bits = 0;
  for (int32_t i=0; i<compiled_shapes.size(); i++) {
    max_seed_len = std::max(max_seed_len, (int64_t) compiled_shapes[i].shape.size());
    max_incl_bits = std::max(max_incl_bits, compiled_shapes[i].num_incl_bits);
  }

  /// The array seq_seeds_starts will hold starting positions for seeds comming from
  /// individual reads. I.e. each read will get a part of seq_seeds_ array which it will
  /// modify, and this array marks the part which is designated for this particular read.
  num_sequences_ = (index_reverse_strand == false) ? seqs.get_sequences().size() : seqs.get_sequences().size()*2;
  num_sequences_forward_ = seqs.get_sequences().size();
  reference_lengths_.resize(num_sequences_);
  headers_.resize(num_sequences_);

  data_.resize(num_sequences_);
  for (int64_t i=0; i<num_sequences_forward_; i++) {
    data_[i].assign(seqs.get_sequences()[i]->get_data(), seqs.get_sequences()[i]->get_data() + seqs.get_sequences()[i]->get_data_length());
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
    }
  }

  int64_t total_num_seeds = 0;
  seq_seed_starts_.clear();
  seq_seed_starts_.resize(seqs.get_sequences().size() * 2);
  seq_seed_counts_.clear();
  seq_seed_counts_.resize(seqs.get_sequences().size() * 2);

  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    reference_lengths_[i] = seqs.get_sequences()[i]->get_sequence_length();
    headers_[i] = seqs.get_sequences()[i]->get_header();
    if (index_reverse_strand == true) {
      reference_lengths_[i + num_sequences_forward_] = reference_lengths_[i];
      headers_[i + num_sequences_forward_] = headers_[i];
    }

    seq_seed_starts_[i] = total_num_seeds;
    if (index_reverse_strand == true) {
      seq_seed_starts_[i + num_sequences_forward_] = total_num_seeds + (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * compiled_shapes.size();
      total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * 2 * compiled_shapes.size();
    } else {
      total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * compiled_shapes.size();
    }
  }

  /// The seed_list_will contain all seeds that were obtained from the input sequences.
  seed_list_.resize(total_num_seeds);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Allocated memory for a list of %ld seeds (128 bits each) (%.5f sec, diff: %.5f sec).\n",
                                                                                          seed_list_.size(), (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), std::string(__FUNCTION__));

  // Process all reads in parallel.
  #pragma omp parallel for num_threads(num_threads) shared(seqs) schedule(dynamic, 1)
  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("\rRead %ld/%ld", (i + 1), seqs.get_sequences().size()), std::string("[]"));
    int8_t *seqdata = (int8_t *) seqs.get_sequences()[i]->get_data();
    int8_t *seqqual = (int8_t *) seqs.get_sequences()[i]->get_quality();
    int64_t seqlen = seqs.get_sequences()[i]->get_data_length();
    uint64_t seq_id = seqs.get_sequences()[i]->get_sequence_absolute_id() + 0; // The '+ 0' is actually short for this: + ((orientation == kReverse) ? num_sequences_forward_ : 0);
    int64_t num_seeds_processed = ProcessKmerSpectrumSeparateReads_(seqdata, seqqual, seqlen, min_avg_seed_qv, seq_id, compiled_shapes, &(seed_list_[seq_seed_starts_[i]]));
    seq_seed_counts_[i] = num_seeds_processed;

    if (index_reverse_strand == true) {
      int8_t *seqrevdata = (int8_t *) seqs.get_sequences()[i]->GetReverseComplement();
      int8_t *seqrevqual = reverse_data(seqqual, seqlen);
      int64_t seqlen = seqs.get_sequences()[i]->get_data_length();
      uint64_t seq_id = seqs.get_sequences()[i]->get_sequence_absolute_id() + num_sequences_forward_; // The '+ num_sequences_forward_' is actually short for this: + ((orientation == kReverse) ? num_sequences_forward_ : 0);

      if (seqrevdata == NULL) continue;
      int64_t num_seeds_processed_rev = ProcessKmerSpectrumSeparateReads_(seqrevdata, seqrevqual, seqlen, min_avg_seed_qv, seq_id, compiled_shapes, &(seed_list_[seq_seed_starts_[i + num_sequences_forward_]]));
      seq_seed_counts_[i+num_sequences_forward_] = num_seeds_processed;
      if (seqrevdata) delete seqrevdata;
      seqrevdata = NULL;
    }
  }

//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Counting occ (%.5f sec, diff: %.5f sec).\n",
//                                                                                          (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
//  diff_time = clock();
//
//  key_occ_.clear();
//  key_occ_.resize(pow(2, max_incl_bits), 0.0f);
//  for (int64_t i=0; i<seed_list_.size(); i++) {
//    key_occ_[seed_list_[i] >> 32] += 1.0f;
//  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Finished constructing the index (%.5f sec, diff: %.5f sec).\n",
                                                                                          (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
  diff_time = clock();

//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("\n"), std::string("[]"));

  return 0;
}

int IndexBrute::CreateFromSequenceFile(const SequenceFile& seqs, float min_avg_seed_qv, const std::vector<std::string> &shapes, bool index_reverse_strand, std::string index_desig, int32_t num_threads) {
  std::vector<CompiledShape> compiled_shapes = CompileShapes(shapes);
  return CreateFromSequenceFile(seqs, min_avg_seed_qv, compiled_shapes, index_reverse_strand, index_desig, num_threads);
}

const std::vector<std::string>& IndexBrute::get_headers() const {
  return headers_;
}

void IndexBrute::set_headers(const std::vector<std::string>& headers) {
  headers_ = headers;
}

int64_t IndexBrute::get_num_sequences() const {
  return num_sequences_;
}

void IndexBrute::set_num_sequences(int64_t numSequences) {
  num_sequences_ = numSequences;
}

int64_t IndexBrute::get_num_sequences_forward() const {
  return num_sequences_forward_;
}

void IndexBrute::set_num_sequences_forward(int64_t numSequencesForward) {
  num_sequences_forward_ = numSequencesForward;
}

const std::vector<int64_t>& IndexBrute::get_reference_lengths() const {
  return reference_lengths_;
}

void IndexBrute::set_reference_lengths(const std::vector<int64_t>& referenceLengths) {
  reference_lengths_ = referenceLengths;
}
