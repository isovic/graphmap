/*
 * index_spaced_hash.cc
 *
 *  Created on: July 11, 2015
 *      Author: ivan
 */

#include <algorithm>
#include <memory>

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "index/index_spaced_hash_fast.h"

IndexSpacedHashFast::IndexSpacedHashFast() {
  data_ = NULL;
  kmer_hash_array_ = NULL;
  kmer_counts_ = NULL;
  all_kmers_ = NULL;
  shape_index_ = NULL;

  Clear();

//  InitShapesPredefined(SHAPE_TYPE_444);
  InitShapesPredefined(SHAPE_TYPE_66);

  //  std::string shapes[] = {
  //                          "111110001111100",
  //                          "111110011111000",
  //                          "111111111100000",
  //                          "111110000111110",
  //                          "111110000011111"
  //                         };
  //  int64_t num_shapes = 5;
  //  std::string shapes[] = {
  //                          "11111101111110",
  //                          "11111100111111",
  //                          "11111111111100",
  //                         };
  //  int64_t num_shapes = 3;
//    std::string shapes[] = {
//                            "111110111110",
//                            "111110011111",
//                            "111111111100",
//                           };
//    int64_t num_shapes = 3;

}

IndexSpacedHashFast::IndexSpacedHashFast(uint32_t shape_type) {
  data_ = NULL;
  kmer_hash_array_ = NULL;
  kmer_counts_ = NULL;
  all_kmers_ = NULL;
  shape_index_ = NULL;

  Clear();

  InitShapesPredefined(shape_type);
}

IndexSpacedHashFast::~IndexSpacedHashFast() {
  Clear();

  shapes_lookup_.clear();
  if (shape_index_)
    free(shape_index_);
  shape_index_ = NULL;
  shape_index_length_ = 0;

}

void IndexSpacedHashFast::Clear() {
  if (kmer_hash_array_)
    free(kmer_hash_array_);
  kmer_hash_array_ = NULL;
  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  if (kmer_counts_)
    free(kmer_counts_);
  kmer_counts_ = NULL;

  reference_starting_pos_.clear();
  reference_lengths_.clear();
  data_length_ = 0;
  data_length_forward_ = 0;
  data_ptr_ = 0;
  num_sequences_ = 0;

  if (data_)
    delete[] data_;
  data_ = NULL;

//  k_ = 15;
//  k_ = 12;
////  k_ = 18;
////  k_ = 15;
////  k_ = 13;
////  k_ = 10;
////  k_ = 8;
//  k_ = 14;
////  k_ = 16;
////  k_ = 10;

//  kmer_hash_.clear();

  num_kmers_ = 0;
  all_kmers_size_ = 0;
}

//int64_t IndexSpacedHashFast::GenerateHashKey(int8_t *seed, uint64_t seed_length) {
//  int64_t ret = 0;
//  int64_t current_accepted_base = 0;
//
//  return GenerateHashKeySplit(seed, seed_length, 1, 2);
//}
//
//int64_t IndexSpacedHashFast::GenerateHashKeySplit(int8_t *seed, uint64_t seed_length, int num_split_spaces, int max_num_split_spaces) {
//  int64_t ret = 0;
//  int64_t current_accepted_base = 0;
//
//  int64_t num_consecutive_bases = (seed_length - max_num_split_spaces) / 2;
//
//  for (uint64_t current_base = 0; current_base < num_consecutive_bases; current_base++) {
//    if (!kIsBase[seed[current_base]])
//      return -1;
//    ret |= (kBaseToBwa[seed[current_base]] << (current_accepted_base * 2));
//    current_accepted_base += 1;
//  }
//
//  for (uint64_t current_base = 0; current_base < num_consecutive_bases; current_base++) {
//    if (!kIsBase[seed[num_consecutive_bases + num_split_spaces + current_base]])
//      return -1;
//    ret |= (kBaseToBwa[seed[num_consecutive_bases + num_split_spaces + current_base]] << (current_accepted_base * 2));
//    current_accepted_base += 1;
//  }
//
//  return ret;
//}

int64_t IndexSpacedHashFast::GenerateHashKeyFromShape(int8_t *seed, const char *shape, int64_t shape_length) const {
  uint64_t ret = 0;
  uint64_t current_accepted_base = 0;

//  int64_t num_consecutive_bases = (seed_length - max_num_split_spaces) / 2;

  for (uint64_t current_base = 0; current_base < shape_length; current_base++) {
    if (shape[current_base] != '1')
      continue;

    if (!kIsBase[seed[current_base]])
      return -1;

    ret |= (((uint64_t) kBaseToBwa[seed[current_base]]) << ((current_accepted_base << 1)));
    current_accepted_base += 1;
  }

  return ((int64_t) ret);
}

//void IndexSpacedHashFast::CountKmers(int8_t *sequence_data, int64_t sequence_length, int k, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) {  // std::vector<int64_t> &ret_kmer_counts) {
//  int64_t hash_key = -1;
//
////  ret_kmer_counts.resize(std::pow(2, (2 * k)));
//  int64_t num_kmers = CalcNumHashKeys();
//  int64_t *kmer_counts = (int64_t *) calloc(sizeof(int64_t), num_kmers);
//
//  for (uint64_t i = 0; i < (sequence_length - k + 1); i++) {
//    int8_t *seed_start = &(sequence_data[i]);
//    hash_key = GenerateHashKey(seed_start, k);
//
//    if (hash_key < 0)
//      continue;
//    kmer_counts[hash_key] += 1;
//  }
//
//  *ret_kmer_counts = kmer_counts;
//  *ret_num_kmers = num_kmers;
//}
//
//int64_t IndexSpacedHashFast::CalcNumHashKeys() {
//  // Every third base from the kmer is left out.
////  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * (2.0f * (float) k_) / 3.0f)));
////  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * (k_ - 1.0f))));
//  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * (shape_index_length_ - 2))));
//
//  return num_sparse_kmers;
//}

void IndexSpacedHashFast::CountKmersFromShape(int8_t *sequence_data, int64_t sequence_length, const char *shape, int64_t shape_length, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) const {  // std::vector<int64_t> &ret_kmer_counts) {
  int64_t hash_key = -1;

//  ret_kmer_counts.resize(std::pow(2, (2 * k)));
  int64_t num_kmers = CalcNumHashKeysFromShape(shape, shape_length);
  int64_t *kmer_counts = (int64_t *) calloc(sizeof(int64_t), num_kmers);

  for (uint64_t i = 0; i < (sequence_length - shape_length + 1); i++) {
    int8_t *seed_start = &(sequence_data[i]);
    hash_key = GenerateHashKeyFromShape(seed_start, shape, shape_length);

    if (hash_key < 0)
      continue;
    kmer_counts[hash_key] += 1;
  }

  *ret_kmer_counts = kmer_counts;
  *ret_num_kmers = num_kmers;
}

int64_t IndexSpacedHashFast::CalcNumHashKeysFromShape(const char *shape, int64_t shape_length) const {
  int64_t num_accepted_bases = 0;

  for (int64_t i=0; i<shape_length; i++) {
    if (shape[i] != '1')
      continue;
    num_accepted_bases += 1;
  }

  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * num_accepted_bases)));

  return num_sparse_kmers;
}

//int IndexSpacedHashFast::FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) {  //, std::vector<int64_t> &return_positions) {
////  printf ("Tu sam 1!\n");
////  fflush(stdout);
//
//  seed_length = k_;
//  *ret_hits = NULL;
//  *ret_start_hit = 0;
//  *ret_num_hits = 0;
//  int64_t num_all_hits = 0;
//
////  printf ("Tu sam 1.1!\n");
////  fflush(stdout);
//
//
//  int64_t hash_key_mismatch = GenerateHashKeySplit(seed, k_, 1, 2);
//
////  printf ("Tu sam 1.2!\n");
////  fflush(stdout);
//
//  int64_t hash_key_insertion = GenerateHashKeySplit(seed, k_, 2, 2);
//
////  printf ("Tu sam 1.3!\n");
////  fflush(stdout);
//
//  int64_t hash_key_deletion = GenerateHashKeySplit(seed, k_, 0, 2);
//
////  printf ("Tu sam 1.4!\n");
////  fflush(stdout);
//
//  if (hash_key_mismatch < 0 && hash_key_insertion < 0 && hash_key_deletion < 0) {
//    *ret_hits = NULL;
//    *ret_start_hit = 0;
//    *ret_num_hits = 0;
////    printf ("return 3\n");
////    fflush(stdout);
//    return 3;
//  }
//
////  if (hash_key > kmer_hash_.size()) {
////    printf ("%s = %ld\n", GetSubstring((char *) seed, seed_length).c_str(), hash_key);
////    fflush(stdout);
////  }
//
////  printf ("Tu sam 1.5!\n");
////  fflush(stdout);
////
////  if (hash_key_mismatch > kmer_hash_.size()) {
////    printf ("hash_key_mismatch > kmer_hash_.size(), hash_key_mismatch = %ld, kmer_hash_.size() = %ld\n", hash_key_mismatch, kmer_hash_.size());
////    fflush(stdout);
////  }
////  if (hash_key_insertion > kmer_hash_.size()) {
////    printf ("hash_key_insertion > kmer_hash_.size(), hash_key_insertion = %ld, kmer_hash_.size() = %ld\n", hash_key_insertion, kmer_hash_.size());
////    fflush(stdout);
////  }
////  if (hash_key_deletion > kmer_hash_.size()) {
////    printf ("hash_key_deletion > kmer_hash_.size(), hash_key_deletion = %ld, kmer_hash_.size() = %ld\n", hash_key_deletion, kmer_hash_.size());
////    fflush(stdout);
////  }
//
////  num_all_hits = kmer_hash_[hash_key_mismatch].size() + kmer_hash_[hash_key_insertion].size() + kmer_hash_[hash_key_deletion].size();
//  if (hash_key_mismatch >= 0) {
//    num_all_hits += kmer_counts_[hash_key_mismatch];
////    printf ("kmer_hash_[hash_key_mismatch].size() = %ld\thash_key_mismatch = %ld\n", kmer_hash_[hash_key_mismatch].size(), hash_key_mismatch);
//  }
//  if (hash_key_insertion >= 0) {
//    num_all_hits += kmer_counts_[hash_key_insertion];
////    printf ("kmer_hash_[hash_key_insertion].size() = %ld\thash_key_insertion = %ld\n", kmer_hash_[hash_key_insertion].size(), hash_key_insertion);
//  }
//  if (hash_key_deletion >= 0) {
//    num_all_hits += kmer_counts_[hash_key_deletion];
////    printf ("kmer_hash_[hash_key_deletion].size() = %ld\thash_key_deletion = %ld\n", kmer_hash_[hash_key_deletion].size(), hash_key_deletion);
//  }
//
////  printf ("\n");
//
//  int64_t *all_hits = new int64_t[num_all_hits + 1];
//
////  printf ("num_all_hits = %ld\n", num_all_hits);
////  fflush(stdout);
//
//  int64_t current_data_ptr = 0;
//  if (hash_key_mismatch >= 0 && kmer_counts_[hash_key_mismatch] > 0) {
//    memmove(&(all_hits[current_data_ptr]), &(kmer_hash_array_[hash_key_mismatch][0]), kmer_counts_[hash_key_mismatch] * sizeof(int64_t));
//    current_data_ptr += kmer_counts_[hash_key_mismatch];
//
////    for (int i=0; i<kmer_hash_[hash_key_mismatch].size(); i++) {
////      printf ("%ld == %ld, num_all_hits = %ld, kmer_hash_[hash_key_mismatch].size() = %ld\n", kmer_hash_[hash_key_mismatch][i], all_hits[i], num_all_hits, kmer_hash_[hash_key_mismatch].size());
////    }
//  }
//  if (hash_key_insertion >= 0 && kmer_counts_[hash_key_insertion] > 0) {
//    memmove(&(all_hits[current_data_ptr]), &(kmer_hash_array_[hash_key_insertion][0]), kmer_counts_[hash_key_insertion] * sizeof(int64_t));
//    current_data_ptr += kmer_counts_[hash_key_insertion];
//  }
//  if (hash_key_deletion >= 0 && kmer_counts_[hash_key_deletion] > 0) {
//    memmove(&(all_hits[current_data_ptr]), &(kmer_hash_array_[hash_key_deletion][0]), kmer_counts_[hash_key_deletion] * sizeof(int64_t));
//    current_data_ptr += kmer_counts_[hash_key_deletion];
//  }
//
//  std::sort(all_hits, all_hits + num_all_hits);
////  printf ("Tu sam 2!\n");
////  fflush(stdout);
//
//  *ret_hits = all_hits;
//  *ret_start_hit = 0;
//  *ret_num_hits = current_data_ptr;
//
////  printf ("Tu sam 3!\n");
////  fflush(stdout);
//
//  if (num_all_hits == 0) {
//    return 1;
//  }
//
////  printf ("Tu sam 4!\n");
////  fflush(stdout);
//
//  if ((max_num_of_hits > 0 && num_all_hits > ((int64_t) max_num_of_hits))) {
//    return 2;
//  }
//
////  printf ("Tu sam 5!\n");
////  fflush(stdout);
//
//  return 0;
//}

//int IndexSpacedHashFast::CreateIndex_(int8_t *data, uint64_t data_length) {
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Creating spaced hash index.\n"), "CreateIndex_");
//
//  if (kmer_hash_array_)
//    free(kmer_hash_array_);
//  kmer_hash_array_ = NULL;
//  if (all_kmers_)
//    free(all_kmers_);
//  all_kmers_ = NULL;
//  if (kmer_counts_)
//    free(kmer_counts_);
//  kmer_counts_ = NULL;
//
//  int64_t num_kmers = 0;
//  CountKmers(data_, data_length_, k_, &kmer_counts_, &num_kmers);
//  int64_t *kmer_countdown = (int64_t *) malloc(sizeof(int64_t) * num_kmers);
//  memmove(kmer_countdown, kmer_counts_, sizeof(int64_t) * num_kmers);
//  num_kmers_ = num_kmers;
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Kmer counting finished (kmer_counts.size() = %ld)\n", num_kmers_), "CreateIndex_");
//
//  int64_t total_num_kmers = 0;
//  for (uint64_t i = 0; i < num_kmers; i++) {
////    kmer_hash_[i].reserve(kmer_counts[i]);
////    kmer_hash_array_[i] = (int64_t *) malloc(sizeof(int64_t) * kmer_counts_[i]);
//    kmer_countdown[i] -= 1;
//    total_num_kmers += kmer_counts_[i];
////    printf ("kmer_counts_[%ld] = %ld,\ttotal_num_kmers = %ld\n", i, kmer_counts_[i], total_num_kmers);
////    fflush(stdout);
//  }
//
////  kmer_hash_.clear();
////  kmer_hash_.resize(num_kmers);
//  kmer_hash_array_ = (int64_t **) malloc(sizeof(int64_t *) * num_kmers);
//  all_kmers_ = (int64_t *) malloc(sizeof(int64_t) * total_num_kmers);
//  all_kmers_size_ = total_num_kmers;
//
//  int64_t kmer_hash_ptr = 0;
//  for (uint64_t i = 0; i < num_kmers; i++) {
//    if (kmer_counts_[i] > 0)
//      kmer_hash_array_[i] = (all_kmers_ + kmer_hash_ptr);
//    else
//      kmer_hash_array_[i] = NULL;
//    kmer_hash_ptr += kmer_counts_[i];
//  }
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Index memory allocated.\n"), "CreateIndex_");
//
//  int64_t hash_key = -1;
//
//  for (uint64_t i = 0; i < (data_length_ - k_ + 1); i++) {
//    int8_t *seed_start = &(data_[i]);
//    hash_key = GenerateHashKey(seed_start, k_);
//
////    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("%s = %ld = %X\n", GetSubstring((char *) seed_start, k_).c_str(), hash_key, hash_key), "[]");
//
//    if (hash_key < 0)
//      continue;
//
//    kmer_hash_array_[hash_key][kmer_countdown[hash_key]] = ((int64_t) i);
//    kmer_countdown[hash_key] -= 1;
//  }
//
//  if (kmer_countdown)
//    free(kmer_countdown);
//  kmer_countdown = NULL;
//
//  ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Finished creating spaced hash index.\n"), "CreateIndex_");
//
//  return 0;
//}

int IndexSpacedHashFast::FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const {  //, std::vector<int64_t> &return_positions) {
//  printf ("Tu sam 1!\n");
//  fflush(stdout);

  seed_length = shape_index_length_;
  *ret_hits = NULL;
  *ret_start_hit = 0;
  *ret_num_hits = 0;

//  printf ("Tu sam 1.1!\n");
//  fflush(stdout);



  int64_t current_data_ptr = 0;
  int64_t *all_hits = NULL;

  for (int64_t i=0; i<shapes_lookup_.size(); i++) {
    int64_t hash_key = GenerateHashKeyFromShape(seed, shapes_lookup_[i].c_str(), shapes_lookup_[i].size());

    if (hash_key >= 0 && hash_key < num_kmers_ && kmer_counts_[hash_key] > 0) {
      if (all_hits == NULL)
        all_hits = (int64_t *) malloc(sizeof(int64_t) * (current_data_ptr + kmer_counts_[hash_key]));
      else
        all_hits = (int64_t *) realloc(all_hits, (sizeof(int64_t) * (current_data_ptr + kmer_counts_[hash_key])));
      memmove(&(all_hits[current_data_ptr]), &(kmer_hash_array_[hash_key][0]), kmer_counts_[hash_key] * sizeof(int64_t));
      current_data_ptr += kmer_counts_[hash_key];
    }
  }

//  std::sort(all_hits, all_hits + current_data_ptr);

  *ret_hits = all_hits;
  *ret_start_hit = 0;
  *ret_num_hits = current_data_ptr;

  if (current_data_ptr == 0) {
    *ret_hits = NULL;
    *ret_start_hit = 0;
    *ret_num_hits = 0;
    return 1;
  }

  if ((max_num_of_hits > 0 && current_data_ptr > ((int64_t) max_num_of_hits))) {
    return 2;
  }

  return 0;

}

//int64_t IndexSpacedHashFast::RawPositionToReferenceIndexWithReverse(int64_t raw_position) const {
//  /// int64_t coded_position = (int64_t) (local_pos << 32) | ref_id;
//  return ((int64_t) (raw_position & MASK_REF_ID));
//}

//int64_t IndexSpacedHashFast::RawPositionConverter(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position, int64_t *ret_relative_position, SeqOrientation *ret_orientation, int64_t *ret_reference_index_with_reverse) const {
//  //  if (raw_position < 0 || raw_position >= data_length_)
//  //    return -2;
//
//    int64_t reference_index = (int64_t) (((uint64_t) raw_position) & MASK_REF_ID);
//    if (reference_index < 0) {
//      LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "1Offending variable: reference_index. Values: reference_index = %ld, raw_position = %ld, data_length = %ld.", reference_index, raw_position, data_length_));
//      return reference_index;
//    }
//
//    int64_t relative_pos = raw_position >> 32; // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
//    int64_t abs_pos = relative_pos + reference_starting_pos_[(uint64_t) reference_index];
//    SeqOrientation orientation = kForward;
//
//  if (((uint64_t) reference_index) >= num_sequences_forward_) {
//    // Relative position has to be changed, because, from the outside, it is expected that we have reverse-complemented the seed and not the reference sequence.
////    relative_pos = reference_lengths_[(uint64_t) reference_index] - relative_pos - query_length - 1 - reference_index;      // The '-1' is to compensate for the '!' character added at the end of every sequence in the data array.
//    relative_pos = reference_lengths_[(uint64_t) reference_index] - relative_pos - query_length - 1;
//
//    // Unlike BWA, we haven't reversed the order of sequences when their reverse complements were added to the index. That is why we only need to subtract the number of forward sequences, and not do (2*num_forward_sequences - gene_idx - 1).
//    reference_index = reference_index - ((int64_t) num_sequences_forward_);
//    orientation = kReverse;
//  }
//
//  if (ret_reference_index_with_reverse != NULL)
//    *ret_reference_index_with_reverse = reference_index;
//
//  if (ret_absolute_position != NULL)
//    *ret_absolute_position = abs_pos;
//
//  if (ret_relative_position != NULL)
//    *ret_relative_position = relative_pos;
//
//  if (ret_orientation != NULL)
//    *ret_orientation = orientation;
//
//  return reference_index;
//}

int64_t IndexSpacedHashFast::RawPositionConverter2(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position, int64_t *ret_relative_position, SeqOrientation *ret_orientation, int64_t *ret_reference_index_with_reverse) const {
//  if (raw_position < 0 || raw_position >= data_length_)
//    return -2;

  int64_t reference_index = (int64_t) (raw_position & MASK_REF_ID);
  if (reference_index < 0) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "1Offending variable: reference_index. Values: reference_index = %ld, raw_position = %ld, data_length = %ld.", reference_index, raw_position, data_length_));
    return reference_index;
  }

  int64_t relative_pos = raw_position >> 32; // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
  int64_t abs_pos = relative_pos + reference_starting_pos_[(uint64_t) reference_index];
  SeqOrientation orientation = kForward;

  if (((uint64_t) reference_index) >= num_sequences_forward_) {
    orientation = kReverse;
  }

  if (ret_reference_index_with_reverse != NULL)
    *ret_reference_index_with_reverse = reference_index;

  if (ret_absolute_position != NULL)
    *ret_absolute_position = abs_pos;

  if (ret_relative_position != NULL)
    *ret_relative_position = relative_pos;

  if (ret_orientation != NULL)
    *ret_orientation = orientation;

  return reference_index;
}

int IndexSpacedHashFast::FindAllRawPositionsOfSeedKey(int64_t hash_key, int64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const {  //, std::vector<int64_t> &return_positions) {
  int64_t *all_hits = NULL;
  int64_t num_hits = 0;

  if (hash_key >= 0 && hash_key < num_kmers_ && kmer_counts_[hash_key] > 0) {
    all_hits = kmer_hash_array_[hash_key];
    num_hits =  kmer_counts_[hash_key];
  }

  *ret_hits = all_hits;
  *ret_start_hit = 0;
  *ret_num_hits = num_hits;

  if (ret_hits == NULL) {
    return 1;
  }

  if ((max_num_of_hits > 0 && num_hits > ((int64_t) max_num_of_hits))) {
    return 2;
  }

  return 0;
}

int IndexSpacedHashFast::CreateIndex_(int8_t *data, uint64_t data_length) {
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Creating spaced hash index.\n"), "CreateIndex_");

  if (kmer_hash_array_)
    free(kmer_hash_array_);
  kmer_hash_array_ = NULL;
  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  if (kmer_counts_)
    free(kmer_counts_);
  kmer_counts_ = NULL;

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_), "CreateIndex_");

  int64_t num_kmers = 0;
  CountKmersFromShape(data_, data_length_, shape_index_, shape_index_length_, &kmer_counts_, &num_kmers);
  int64_t *kmer_countdown = (int64_t *) malloc(sizeof(int64_t) * num_kmers);
  memmove(kmer_countdown, kmer_counts_, sizeof(int64_t) * num_kmers);
  num_kmers_ = num_kmers;

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Kmer counting finished (kmer_counts.size() = %ld)\n", num_kmers_), "CreateIndex_");

  int64_t total_num_kmers = 0;
  for (uint64_t i = 0; i < num_kmers; i++) {
//    kmer_hash_[i].reserve(kmer_counts[i]);
//    kmer_hash_array_[i] = (int64_t *) malloc(sizeof(int64_t) * kmer_counts_[i]);
//    kmer_countdown[i] -= 1;
    kmer_countdown[i] = 0;
    total_num_kmers += kmer_counts_[i];
//    printf ("kmer_counts_[%ld] = %ld,\ttotal_num_kmers = %ld\n", i, kmer_counts_[i], total_num_kmers);
//    fflush(stdout);
  }

//  kmer_hash_.clear();
//  kmer_hash_.resize(num_kmers);
  kmer_hash_array_ = (int64_t **) malloc(sizeof(int64_t *) * num_kmers);
  all_kmers_ = (int64_t *) malloc(sizeof(int64_t) * total_num_kmers);
  all_kmers_size_ = total_num_kmers;

  int64_t kmer_hash_ptr = 0;
  for (uint64_t i = 0; i < num_kmers; i++) {
    if (kmer_counts_[i] > 0)
      kmer_hash_array_[i] = (all_kmers_ + kmer_hash_ptr);
    else
      kmer_hash_array_[i] = NULL;
    kmer_hash_ptr += kmer_counts_[i];
  }

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Index memory allocated.\n"), "CreateIndex_");

  int64_t hash_key = -1;

  uint64_t current_ref_id = 0;

  for (uint64_t i = 0; i < (data_length_ - shape_index_length_ + 1); i++) {
    if (i >= (reference_starting_pos_[current_ref_id] + reference_lengths_[current_ref_id]))
      current_ref_id += 1;

    int8_t *seed_start = &(data_[i]);
    hash_key = GenerateHashKeyFromShape(seed_start, shape_index_, shape_index_length_);

//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("%s = %ld = %X\n", GetSubstring((char *) seed_start, k_).c_str(), hash_key, hash_key), "[]");

    if (hash_key < 0)
      continue;

//    int64_t coded_position = (int64_t) (((uint64_t) i) & ((uint64_t) 0xFFFFFFFF)) | (((uint64_t) current_ref_id) << 32);
    uint64_t local_pos = ((uint64_t) (i - reference_starting_pos_[current_ref_id])) & ((uint64_t) 0x00000000FFFFFFFF);
    uint64_t ref_id = ((uint64_t) current_ref_id) & ((uint64_t) 0x00000000FFFFFFFF);
    int64_t coded_position = (int64_t) (local_pos << 32) | ref_id;

    kmer_hash_array_[hash_key][kmer_countdown[hash_key]] = coded_position;
//    kmer_countdown[hash_key] -= 1;
    kmer_countdown[hash_key] += 1;
  }

  if (kmer_countdown)
    free(kmer_countdown);
  kmer_countdown = NULL;

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Finished creating spaced hash index.\n"), "CreateIndex_");

  return 0;
}



int IndexSpacedHashFast::SerializeIndex_(FILE* fp_out) {
  int64_t vector_length = 0;
  fwrite(&shape_index_length_, sizeof(int64_t), 1, fp_out);
  fwrite(shape_index_, sizeof(char), shape_index_length_, fp_out);

  fwrite(&num_kmers_, sizeof(int64_t), 1, fp_out);

  fwrite(kmer_counts_, sizeof(int64_t), num_kmers_, fp_out);
  fwrite(&all_kmers_size_, sizeof(int64_t), 1, fp_out);
  fwrite(all_kmers_, sizeof(int64_t), all_kmers_size_, fp_out);

//
//  for (int64_t i=0; i<num_kmers_; i++) {
//    vector_length = kmer_counts_[i];
//    if (vector_length > 0)
//      fwrite(&(kmer_hash_array_[i][0]), sizeof(int64_t), vector_length, fp_out);
//  }

  return 0;
}

int IndexSpacedHashFast::IsManualCleanupRequired(std::string function_name) const {
  if (function_name == "FindAllRawPositionsOfSeed")
    return 0;

  return 1;
}

char* IndexSpacedHashFast::get_shape_index() const {
  return shape_index_;
}

void IndexSpacedHashFast::set_shape_index(char* shapeIndex) {
  shape_index_ = shapeIndex;
}

int64_t IndexSpacedHashFast::get_shape_index_length() const {
  return shape_index_length_;
}

void IndexSpacedHashFast::set_shape_index_length(int64_t shapeIndexLength) {
  shape_index_length_ = shapeIndexLength;
}

//int IndexSpacedHashFast::get_k() const {
//  return k_;
//}
//
//void IndexSpacedHashFast::set_k(int k) {
//  k_ = k;
//}

int IndexSpacedHashFast::DeserializeIndex_(FILE* fp_in) {
  if (shape_index_)
    free(shape_index_);
  shape_index_ = NULL;
  shape_index_length_ = 0;
  if (kmer_hash_array_)
    free(kmer_hash_array_);
  kmer_hash_array_ = NULL;
  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  if (kmer_counts_)
    free(kmer_counts_);
  kmer_counts_ = NULL;

  int64_t vector_length = 0;

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- k_...\n"), "DeserializeIndex_");
  if (fread(&shape_index_length_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable shape_index_length_."));
    return 1;
  }

  shape_index_ = (char *) malloc (sizeof(char) * (shape_index_length_ + 1));
  if (fread(shape_index_, sizeof(char), shape_index_length_, fp_in) != shape_index_length_) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable shape_index_."));
    return 1;
  }
  shape_index_[shape_index_length_] = '\0';

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_), "DeserializeIndex_");

//  printf ("shapes_lookup_.size() = %ld\n", shapes_lookup_.size());
//  for (int64_t i=0; i<shapes_lookup_.size(); i++)
//    printf ("shapes_lookup[%d] = %s\n", i, shapes_lookup_[i].c_str());
//  fflush(stdout);



  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- num_kmers_...\n"), "DeserializeIndex_");
  if (fread(&num_kmers_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_kmers_."));
    return 1;
  }

  if (num_kmers_ <= 0) {
    return 1;
  }

  kmer_counts_ = (int64_t *) malloc(sizeof(int64_t) * num_kmers_);
  if (fread(kmer_counts_, sizeof(int64_t), num_kmers_, fp_in) != num_kmers_) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable kmer_counts_.\n"));
    return 3;
  }

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- all_kmers_size_...\n"), "DeserializeIndex_");
  if (fread(&all_kmers_size_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable all_kmers_size_.\n"));
    return 1;
  }

  all_kmers_ = (int64_t *) malloc(sizeof(int64_t) * all_kmers_size_);
  if (fread(all_kmers_, sizeof(int64_t), all_kmers_size_, fp_in) != all_kmers_size_) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable all_kmers.\n"));
    return 3;
  }

  kmer_hash_array_ = (int64_t **) malloc(sizeof(int64_t *) * num_kmers_);
  int64_t kmer_ptr = 0;
  for (int64_t i = 0; i < num_kmers_; i++) {
    if (kmer_counts_[i] > 0)
      kmer_hash_array_[i] = (all_kmers_ + kmer_ptr);
    else
      kmer_hash_array_[i] = NULL;
    kmer_ptr += kmer_counts_[i];
  }

  return 0;
}

void IndexSpacedHashFast::Verbose(FILE* fp) const {
//  fprintf (fp, "Num sequences forward: %ld\n", num_sequences_forward_);
//  fprintf (fp, "Num sequences: %ld\n", num_sequences_);
//  fprintf (fp, "Data length forward: %ld\n", data_length_forward_);
//  fprintf (fp, "Data length: %ld\n", data_length_);
//  fprintf (fp, "\n");
//  for (uint64_t i=0; i<num_sequences_forward_; i++) {
//    fprintf (fp, "Header: '%s'\n", headers_[i].c_str());
//    fprintf (fp, "Sequence start: %ld\n", reference_starting_pos_[i]);
//    fprintf (fp, "Sequence length: %ld\n", reference_lengths_[i]);
//  }
}

std::string IndexSpacedHashFast::VerboseToString() const {
  return (std::string(""));
}

int IndexSpacedHashFast::InitShapesPredefined(uint32_t shape_type) {
//  std::string shape_temp = "11111011111";
//  std::vector<std::string> shapes_lookup_temp;
//  shapes_lookup_temp.push_back("111110111110");
//  shapes_lookup_temp.push_back("111110011111");
//  shapes_lookup_temp.push_back("111111111100");

//  std::string shape_temp = "1111110111111";
//  std::vector<std::string> shapes_lookup_temp;
//  shapes_lookup_temp.push_back("11111101111110");
//  shapes_lookup_temp.push_back("11111100111111");
//  shapes_lookup_temp.push_back("11111111111100");

//  std::string shape_temp =     "111111000111111";
//  std::vector<std::string> shapes_lookup_temp;
//  shapes_lookup_temp.push_back("11111111111100000");
//  shapes_lookup_temp.push_back("11111101111110000");
//  shapes_lookup_temp.push_back("11111100111111000");
//  shapes_lookup_temp.push_back("11111100011111100");
//  shapes_lookup_temp.push_back("11111100001111110");
//  shapes_lookup_temp.push_back("11111100000111111");
//  std::vector<std::string> shapes_lookup_temp;

  // Ovaj je radio ok na svim datasetovima, ali metagenomni je izgubio na preciznosti (oko 94.5%).
//  std::string shape_temp = "11110111101111";
//  shapes_lookup_temp.push_back("111111111111");
//  shapes_lookup_temp.push_back("1111111101111");
//  shapes_lookup_temp.push_back("11111111001111");
//  shapes_lookup_temp.push_back("1111011111111");
//  shapes_lookup_temp.push_back("11110111101111");
//  shapes_lookup_temp.push_back("111101111001111");
//  shapes_lookup_temp.push_back("11110011111111");
//  shapes_lookup_temp.push_back("111100111101111");
//  shapes_lookup_temp.push_back("1111001111001111");


//  std::string shape_temp =     "111101111101111";
//  shapes_lookup_temp.push_back("1111111111111");
//  shapes_lookup_temp.push_back("11111111101111");
//  shapes_lookup_temp.push_back("111111111001111");
//  shapes_lookup_temp.push_back("11110111111111");
//  shapes_lookup_temp.push_back("111101111101111");
//  shapes_lookup_temp.push_back("1111011111001111");
//  shapes_lookup_temp.push_back("111100111111111");
//  shapes_lookup_temp.push_back("1111001111101111");
//  shapes_lookup_temp.push_back("11110011111001111");

//  std::string shape_temp =     "1111011111101111";
//  shapes_lookup_temp.push_back("11111111111111");
//  shapes_lookup_temp.push_back("111111111101111");
//  shapes_lookup_temp.push_back("1111111111001111");
//  shapes_lookup_temp.push_back("111101111111111");
//  shapes_lookup_temp.push_back("1111011111101111");
//  shapes_lookup_temp.push_back("11110111111001111");
//  shapes_lookup_temp.push_back("1111001111111111");
//  shapes_lookup_temp.push_back("11110011111101111");
//  shapes_lookup_temp.push_back("111100111111001111");

//  std::string shape_temp =     "11101110111";
//  shapes_lookup_temp.push_back("111111111");
//  shapes_lookup_temp.push_back("1111110111");
//  shapes_lookup_temp.push_back("11111100111");
//  shapes_lookup_temp.push_back("1110111111");
//  shapes_lookup_temp.push_back("11101110111");
//  shapes_lookup_temp.push_back("111011100111");
//  shapes_lookup_temp.push_back("11100111111");
//  shapes_lookup_temp.push_back("111001110111");
//  shapes_lookup_temp.push_back("1110011100111");

  std::string shape_for_indexing = "";
  std::vector<std::string> shapes_for_search;

  if (shape_type == SHAPE_TYPE_444) {
    shape_for_indexing = "11110111101111";
    shapes_for_search = {"111111111111",   "1111111101111",    "11111111001111",
                         "1111011111111",  "11110111101111",   "111101111001111",
                         "11110011111111", "111100111101111",  "1111001111001111"};
  } else {  // SHAPE_TYPE_66
    shape_for_indexing = "1111110111111";
    shapes_for_search = {"111111111111", "1111110111111", "11111100111111"};

  }


  if (shape_for_indexing.size() == 0 || shapes_for_search.size() == 0) {
    return 1;
  }

  shape_index_length_ = shape_for_indexing.size();
  shape_index_ = (char *) malloc(sizeof(char) * (shape_index_length_ + 1));
  memmove(shape_index_, shape_for_indexing.c_str(), shape_index_length_);
  shape_index_[shape_index_length_] = '\0';

  shapes_lookup_ = shapes_for_search;

  return 0;
}



// Example usage:
//std::string shape_for_indexing_66 = "1111110111111";
//std::vector<std::string> shapes_for_search_66 = {"111111111111", "1111110111111", "11111100111111"};
//
//std::string shape_for_indexing_444 = "11110111101111";
//std::vector<std::string> shapes_for_search_444 = {"111111111111",   "1111111101111",    "11111111001111",
//                                                  "1111011111111",  "11110111101111",   "111101111001111",
//                                                  "11110011111111", "111100111101111",  "1111001111001111"};
//
//index66.InitShapes(shape_for_indexing_66, shapes_for_search_66);
//index444.InitShapes(shape_for_indexing_444, shapes_for_search_444);
int IndexSpacedHashFast::InitShapes(std::string shape_for_indexing, std::vector<std::string> &shapes_for_search) {
  shape_index_length_ = shape_for_indexing.size();
  shape_index_ = (char *) malloc(sizeof(char) * (shape_index_length_ + 1));
  memmove(shape_index_, shape_for_indexing.c_str(), shape_index_length_);
  shape_index_[shape_index_length_] = '\0';

  shapes_lookup_ = shapes_for_search;

  return 0;
}
