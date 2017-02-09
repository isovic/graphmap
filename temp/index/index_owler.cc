/*
 * index_owler.cc
 *
 *  Created on: July 11, 2015
 *      Author: ivan
 */

#include <algorithm>
#include <memory>

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "index/index_owler.h"

IndexOwler::IndexOwler() {
  data_ = NULL;
  kmer_hash_array_ = NULL;
  kmer_counts_ = NULL;
  all_kmers_ = NULL;
  shape_index_ = NULL;
  read_subindex_ = NULL;
  all_subindexes_ = NULL;
  subindex_counts_ = NULL;
  all_subindexes_size_ = 0;

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

IndexOwler::IndexOwler(uint32_t shape_type) {
  data_ = NULL;
  kmer_hash_array_ = NULL;
  kmer_counts_ = NULL;
  all_kmers_ = NULL;
  shape_index_ = NULL;
  read_subindex_ = NULL;
  all_subindexes_ = NULL;
  subindex_counts_ = NULL;
  all_subindexes_size_ = 0;

  Clear();

  InitShapesPredefined(shape_type);
}

IndexOwler::~IndexOwler() {
  Clear();

  shapes_lookup_.clear();
  if (shape_index_)
    free(shape_index_);
  shape_index_ = NULL;
  shape_index_length_ = 0;

}

void IndexOwler::Clear() {
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

  // num_sequences_ counts both forward and reverse sequences.
  all_subindexes_size_ = 0;
  if (all_subindexes_)
    free(all_subindexes_);
  all_subindexes_ = NULL;
  if (subindex_counts_)
    free(subindex_counts_);
  subindex_counts_ = NULL;
  if (read_subindex_)
    free(read_subindex_);
  read_subindex_ = NULL;

//  for (int64_t i = 0; i < num_sequences_; i++) {
//    if (read_subindex_[i].positions) {
//      free(read_subindex_[i].positions);
//      read_subindex_[i].positions = NULL;
//    }
//  }
//  if (read_subindex_)
//    free(read_subindex_);
//  read_subindex_ = NULL;
//  read_subindex_.clear();

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

//int64_t IndexOwler::GenerateHashKey(int8_t *seed, uint64_t seed_length) {
//  int64_t ret = 0;
//  int64_t current_accepted_base = 0;
//
//  return GenerateHashKeySplit(seed, seed_length, 1, 2);
//}
//
//int64_t IndexOwler::GenerateHashKeySplit(int8_t *seed, uint64_t seed_length, int num_split_spaces, int max_num_split_spaces) {
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

int64_t IndexOwler::GenerateHashKeyFromShape(int8_t *seed, const char *shape, int64_t shape_length) const {
  uint64_t ret = 0;
  uint64_t current_accepted_base = 0;

//  int64_t num_consecutive_bases = (seed_length - max_num_split_spaces) / 2;

  for (uint64_t current_base = 0; current_base < shape_length; current_base++) {
    if (!kIsBase[seed[current_base]])
      return -1;

    if (shape[current_base] != '1')
      continue;

    ret |= (((uint64_t) kBaseToBwa[seed[current_base]]) << ((current_accepted_base << 1)));
    current_accepted_base += 1;
  }

  return ((int64_t) ret);
}

//void IndexOwler::CountKmers(int8_t *sequence_data, int64_t sequence_length, int k, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) {  // std::vector<int64_t> &ret_kmer_counts) {
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
//int64_t IndexOwler::CalcNumHashKeys() {
//  // Every third base from the kmer is left out.
////  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * (2.0f * (float) k_) / 3.0f)));
////  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * (k_ - 1.0f))));
//  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * (shape_index_length_ - 2))));
//
//  return num_sparse_kmers;
//}

void IndexOwler::CountKmersFromShape(int8_t *sequence_data, int64_t sequence_length, const char *shape, int64_t shape_length, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) const {  // std::vector<int64_t> &ret_kmer_counts) {
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

int64_t IndexOwler::CalcNumHashKeysFromShape(const char *shape, int64_t shape_length) const {
  int64_t num_accepted_bases = 0;

  for (int64_t i=0; i<shape_length; i++) {
    if (shape[i] != '1')
      continue;
    num_accepted_bases += 1;
  }

  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * num_accepted_bases)));

  return num_sparse_kmers;
}

//int IndexOwler::FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) {  //, std::vector<int64_t> &return_positions) {
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

//int IndexOwler::CreateIndex_(int8_t *data, uint64_t data_length) {
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

int IndexOwler::FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const {  //, std::vector<int64_t> &return_positions) {
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

int IndexOwler::FindAllRawPositionsOfSeedKey(int64_t hash_key, int64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const {  //, std::vector<int64_t> &return_positions) {
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

int IndexOwler::CreateIndex_(int8_t *data, uint64_t data_length) {
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Creating double hashed spaced hash index.\n"), "CreateIndex_");

  if (kmer_hash_array_)
    free(kmer_hash_array_);
  kmer_hash_array_ = NULL;
  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  if (kmer_counts_)
    free(kmer_counts_);
  kmer_counts_ = NULL;

  all_subindexes_size_ = 0;
  if (all_subindexes_)
    free(all_subindexes_);
  all_subindexes_ = NULL;
  if (subindex_counts_)
    free(subindex_counts_);
  subindex_counts_ = NULL;
  if (read_subindex_)
    free(read_subindex_);
  read_subindex_ = NULL;

//  for (int64_t i = 0; i < num_sequences_; i++) {
//    if (read_subindex_[i].positions) {
//      free(read_subindex_[i].positions);
//      read_subindex_[i].positions = NULL;
//    }
//  }
//  if (read_subindex_)
//    free(read_subindex_);
//  read_subindex_ = NULL;
//  read_subindex_.clear();

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_), "CreateIndex_");

  int64_t num_kmers = 0;
  CountKmersFromShape(data_, data_length_, shape_index_, shape_index_length_, &kmer_counts_, &num_kmers);
  int64_t *kmer_countdown = (int64_t *) malloc(sizeof(int64_t) * num_kmers);
  memmove(kmer_countdown, kmer_counts_, sizeof(int64_t) * num_kmers);
  num_kmers_ = num_kmers;

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Kmer counting finished (kmer_counts.size() = %ld)\n", num_kmers_), "CreateIndex_");

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

//  read_subindex_ = (SubIndex *) malloc(sizeof(SubIndex) * num_sequences_);
//  read_subindex_.resize(num_sequences_);
  std::vector<std::vector<SubIndex> > read_subindex;

  read_subindex.resize(num_sequences_);

  int64_t kmer_hash_ptr = 0;
  for (uint64_t i = 0; i < num_kmers; i++) {
    if (kmer_counts_[i] > 0)
      kmer_hash_array_[i] = (all_kmers_ + kmer_hash_ptr);
    else
      kmer_hash_array_[i] = NULL;
    kmer_hash_ptr += kmer_counts_[i];
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Index memory allocated.\n"), "CreateIndex_");

  int64_t hash_key = -1;

  uint64_t current_ref_id = 0;
//  printf ("current_ref_id = %ld\n", current_ref_id);
//  fflush(stdout);

  for (uint64_t i = 0; i < (data_length_ - shape_index_length_ + 1); i++) {
    if (i >= (reference_starting_pos_[current_ref_id] + reference_lengths_[current_ref_id])) {
      current_ref_id += 1;
//      printf ("current_ref_id = %ld (i = %ld)\n", current_ref_id, i);
//      fflush(stdout);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\rProcessed %.2f%%", (((float) i) / ((float) (data_length_ - shape_index_length_ + 1))) * 100.0f), "[]");
    }

    int8_t *seed_start = &(data_[i]);
    hash_key = GenerateHashKeyFromShape(seed_start, shape_index_, shape_index_length_);

//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("%s = %ld = %X\n", GetSubstring((char *) seed_start, k_).c_str(), hash_key, hash_key), "[]");

    if (hash_key < 0)
      continue;

//    int64_t coded_position = (int64_t) (((uint64_t) i) & ((uint64_t) 0xFFFFFFFF)) | (((uint64_t) current_ref_id) << 32);
    uint64_t local_pos = ((uint64_t) (i - reference_starting_pos_[current_ref_id])) & ((uint64_t) 0x00000000FFFFFFFF);
    uint64_t ref_id = ((uint64_t) current_ref_id) & ((uint64_t) 0x00000000FFFFFFFF);
    int64_t coded_position = (int64_t) (local_pos << 32) | ref_id;

//    if (local_pos > reference_lengths_[current_ref_id]) {
//      printf ("Tu sam 1! i = %ld, local_pos = %ld, current_ref_id = %ld, reference_lengths_[current_ref_id] = %ld\n", i, local_pos, current_ref_id, reference_lengths_[current_ref_id]);
//      fflush(stdout);
//    }

    kmer_hash_array_[hash_key][kmer_countdown[hash_key]] = coded_position;
//    kmer_countdown[hash_key] -= 1;
    kmer_countdown[hash_key] += 1;



    /// Generate read subindex. This generates all keys present in the read, and relates them to a position.
    SubIndex subindex;
    subindex.key = (uint32_t) (((uint64_t) hash_key) & ((uint64_t) 0x00000000FFFFFFFF));
    subindex.position = (uint32_t) local_pos;
    read_subindex[current_ref_id].push_back(subindex);

//    for (int j = 0; j < shapes_lookup_.size(); j++) {
//      if (shapes_lookup_[j] != std::string(shape_index_)) {
//        int64_t alternate_hash_key = GenerateHashKeyFromShape(seed_start, shapes_lookup_[j].c_str(), shapes_lookup_[j].size());
//        if (alternate_hash_key < 0)
//          continue;
//
//        SubIndex subindex_alternate;
//        subindex_alternate.key = (uint32_t) (((uint64_t) alternate_hash_key) & ((uint64_t) 0x00000000FFFFFFFF));
//        subindex_alternate.position = (uint32_t) local_pos;
//        read_subindex[current_ref_id].push_back(subindex_alternate);
//      }
//    }
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\rProcesseg 100.00%%"), "[]");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\n"), "[]");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Converting read subindex data structures...\n"), "CreateIndex_");

  if (kmer_countdown)
    free(kmer_countdown);
  kmer_countdown = NULL;

  subindex_counts_ = (int64_t *) malloc(sizeof(int64_t) * num_sequences_);
  int64_t total_num_subindexes = 0;
  for (int64_t i = 0; i < read_subindex.size(); i++) {
    std::sort(read_subindex[i].begin(), read_subindex[i].end(), subindex_less_than_key());
    total_num_subindexes += read_subindex[i].size();
    subindex_counts_[i] = read_subindex[i].size();
  }
  all_subindexes_size_ = total_num_subindexes;
  all_subindexes_ = (SubIndex *) malloc(sizeof(SubIndex) * total_num_subindexes);
  read_subindex_ = (SubIndex **) malloc(sizeof(SubIndex *) * num_sequences_);
  int64_t current_num_subindexes = 0;
  for (int64_t i = 0; i < num_sequences_; i++) {
    if (subindex_counts_[i] > 0) {
      memmove((all_subindexes_ + current_num_subindexes), &(read_subindex[i][0]), sizeof(SubIndex) * subindex_counts_[i]);
      read_subindex[i].clear();
      read_subindex_[i] = (all_subindexes_ + current_num_subindexes);
      current_num_subindexes += subindex_counts_[i];
    } else {
      read_subindex_[i] = NULL;
    }
  }
  read_subindex.clear();

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Finished creating spaced hash index.\n"), "CreateIndex_");

  return 0;
}



int IndexOwler::SerializeIndex_(FILE* fp_out) {
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Started index serialization...\n"), "SerializeIndex_");

  int64_t vector_length = 0;
  fwrite(&shape_index_length_, sizeof(int64_t), 1, fp_out);
  fwrite(shape_index_, sizeof(char), shape_index_length_, fp_out);

  fwrite(&num_kmers_, sizeof(int64_t), 1, fp_out);

  fwrite(kmer_counts_, sizeof(int64_t), num_kmers_, fp_out);
  fwrite(&all_kmers_size_, sizeof(int64_t), 1, fp_out);
  fwrite(all_kmers_, sizeof(int64_t), all_kmers_size_, fp_out);



//  std::vector<int64_t> counts;
//  counts.resize(read_subindex_.size());
//  for (int64_t i = 0; i < read_subindex_.size(); i++) {
//    counts[i] = read_subindex_[i].size();
//  }
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Serializing the read subindex...\n"), "SerializeIndex_");

  fwrite(subindex_counts_, sizeof(int64_t), num_sequences_, fp_out);
  fwrite(&all_subindexes_size_, sizeof(int64_t), 1, fp_out);
  fwrite(all_subindexes_, sizeof(SubIndex), all_subindexes_size_, fp_out);

//  for (int64_t i = 0; i < read_subindex_.size(); i++) {
//    if (read_subindex_[i].size() > 0)
//      fwrite(&(read_subindex_[i][0]), sizeof(SubIndex), read_subindex_[i].size(), fp_out);
//  }

//
//  for (int64_t i=0; i<num_kmers_; i++) {
//    vector_length = kmer_counts_[i];
//    if (vector_length > 0)
//      fwrite(&(kmer_hash_array_[i][0]), sizeof(int64_t), vector_length, fp_out);
//  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Finished serializing the index.\n"), "SerializeIndex_");

  return 0;
}

int IndexOwler::IsManualCleanupRequired(std::string function_name) const {
  if (function_name == "FindAllRawPositionsOfSeed")
    return 0;

  return 1;
}

char* IndexOwler::get_shape_index() const {
  return shape_index_;
}

void IndexOwler::set_shape_index(char* shapeIndex) {
  shape_index_ = shapeIndex;
}

int64_t IndexOwler::get_shape_index_length() const {
  return shape_index_length_;
}

void IndexOwler::set_shape_index_length(int64_t shapeIndexLength) {
  shape_index_length_ = shapeIndexLength;
}

//int IndexOwler::get_k() const {
//  return k_;
//}
//
//void IndexOwler::set_k(int k) {
//  k_ = k;
//}

int IndexOwler::DeserializeIndex_(FILE* fp_in) {
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

  all_subindexes_size_ = 0;
  if (all_subindexes_)
    free(all_subindexes_);
  all_subindexes_ = NULL;
  if (subindex_counts_)
    free(subindex_counts_);
  subindex_counts_ = NULL;
  if (read_subindex_)
    free(read_subindex_);
  read_subindex_ = NULL;

  int64_t vector_length = 0;

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- k_...\n"), "DeserializeIndex_");
  if (fread(&shape_index_length_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable shape_index_length_."));
    return 1;
  }

  shape_index_ = (char *) malloc (sizeof(char) * (shape_index_length_ + 1));
  if (fread(shape_index_, sizeof(char), shape_index_length_, fp_in) != shape_index_length_) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable shape_index_."));
    return 1;
  }
  shape_index_[shape_index_length_] = '\0';

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_), "DeserializeIndex_");

//  printf ("shapes_lookup_.size() = %ld\n", shapes_lookup_.size());
//  for (int64_t i=0; i<shapes_lookup_.size(); i++)
//    printf ("shapes_lookup[%d] = %s\n", i, shapes_lookup_[i].c_str());
//  fflush(stdout);



  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- num_kmers_...\n"), "DeserializeIndex_");
  if (fread(&num_kmers_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable num_kmers_."));
    return 1;
  }

  if (num_kmers_ <= 0) {
    return 1;
  }

  kmer_counts_ = (int64_t *) malloc(sizeof(int64_t) * num_kmers_);
  if (fread(kmer_counts_, sizeof(int64_t), num_kmers_, fp_in) != num_kmers_) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable kmer_counts_.\n"));
    return 3;
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- all_kmers_size_...\n"), "DeserializeIndex_");
  if (fread(&all_kmers_size_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable all_kmers_size_.\n"));
    return 1;
  }

  all_kmers_ = (int64_t *) malloc(sizeof(int64_t) * all_kmers_size_);
  if (fread(all_kmers_, sizeof(int64_t), all_kmers_size_, fp_in) != all_kmers_size_) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable all_kmers.\n"));
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



//  fwrite(subindex_counts_, sizeof(int64_t), num_sequences_, fp_out);
//  fwrite(&all_subindexes_size_, sizeof(int64_t), 1, fp_out);
//  fwrite(all_subindexes_, sizeof(SubIndex), all_subindexes_size_, fp_out);

//  printf ("num_sequences_ = %ld\n", num_sequences_);
//  fflush(stdout);

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- allocating space for read subindex\n"), "DeserializeIndex_");
  read_subindex_ = (SubIndex **) malloc(sizeof(SubIndex *) * num_sequences_);
  subindex_counts_ = (int64_t *) malloc(sizeof(int64_t) * num_sequences_);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reading subindex_counts_...\n"), "DeserializeIndex_");
  if (fread(subindex_counts_, sizeof(int64_t), num_sequences_, fp_in) != num_sequences_) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable subindex_counts_!\n"));
    return 3;
  }
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reading all_subindexes_size_..."), "DeserializeIndex_");
  if (fread(&all_subindexes_size_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable all_subindexes_size_!\n"));
    return 3;
  }
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("%ld\n", all_subindexes_size_), "[]");
//  printf ("all_subindexes_size_ = %ld\n", all_subindexes_size_);
//  fflush(stdout);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- allocating space for all_subindexes_\n"), "DeserializeIndex_");
  all_subindexes_ = (SubIndex *) malloc(sizeof(SubIndex) * all_subindexes_size_);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- reading all_subindexes_...\n"), "DeserializeIndex_");
  if (fread(all_subindexes_, sizeof(SubIndex), all_subindexes_size_, fp_in) != all_subindexes_size_) {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable all_subindexes_!\n"));
    return 3;
  }
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- formatting the read subindex data structures...\n"), "DeserializeIndex_");
  int64_t current_num_subindexes = 0;
  for (int64_t i = 0; i < num_sequences_; i++) {
    if (subindex_counts_[i] > 0) {
      read_subindex_[i] = (all_subindexes_ + current_num_subindexes);
      current_num_subindexes += subindex_counts_[i];
    } else {
      read_subindex_[i] = NULL;
    }
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Finished deserializing the index!\n"), "DeserializeIndex_");

//  int64_t *counts = (int64_t *) malloc(sizeof(int64_t) * num_sequences_);
//  if (fread(counts, sizeof(int64_t), num_sequences_, fp_in) != num_sequences_) {
//    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable counts.\n"));
//    return 3;
//  }
//  read_subindex_.clear();
//  read_subindex_.resize(num_sequences_);
//  for (int64_t i = 0; i < read_subindex_.size(); i++) {
//
//    if (counts[i] > 0)
////      fwrite(&(read_subindex_[i][0]), sizeof(SubIndex), read_subindex_[i].size(), fp_out);
//
//  }
//  if (counts)
//    free(counts);

  return 0;
}

void IndexOwler::Verbose(FILE* fp) const {
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

std::string IndexOwler::VerboseToString() const {
  return (std::string(""));
}

SubIndex** IndexOwler::get_read_subindex() const {
  return read_subindex_;
}

void IndexOwler::set_read_subindex(SubIndex** readSubindex) {
  read_subindex_ = readSubindex;
}

int64_t* IndexOwler::get_subindex_counts() const {
  return subindex_counts_;
}

void IndexOwler::set_subindex_counts(int64_t* subindexCounts) {
  subindex_counts_ = subindexCounts;
}

int IndexOwler::InitShapesPredefined(uint32_t shape_type) {
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
int IndexOwler::InitShapes(std::string shape_for_indexing, std::vector<std::string> &shapes_for_search) {
  shape_index_length_ = shape_for_indexing.size();
  shape_index_ = (char *) malloc(sizeof(char) * (shape_index_length_ + 1));
  memmove(shape_index_, shape_for_indexing.c_str(), shape_index_length_);
  shape_index_[shape_index_length_] = '\0';

  shapes_lookup_ = shapes_for_search;

  return 0;
}
