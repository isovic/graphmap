/*
 * index_hash.cpp
 *
 *  Created on: Oct 21, 2014
 *      Author: ivan
 */

#include <algorithm>
#include "index_hash.h"

IndexHash::IndexHash() {
  Init();
}

IndexHash::~IndexHash() {
  Clear();
}

void IndexHash::Init() {
  k_ = 6;
//  kmer_hash_.clear();
  kmer_hash_ = NULL;
  num_threads_ = 1;
  data_ = NULL;
  suffix_array_ = NULL;
  kmer_hash_size_ = 0;
  kmer_hash_counts_size_ = 0;
  kmer_hash_counts_ = NULL;
  kmer_hash_last_key_ = 0;
  kmer_hash_last_key_initialized_ = false;

  all_kmers_ = NULL;
  all_kmers_size_ = 0;
}

void IndexHash::Clear() {

  if (data_)
    delete[] data_;
  data_ = NULL;

  if (suffix_array_)
    delete suffix_array_;
  suffix_array_ = NULL;

  if (kmer_hash_counts_)
    free(kmer_hash_counts_);
  kmer_hash_counts_ = NULL;

//  if (kmer_hash_ != NULL) {
//    for (int64_t i=0; i<kmer_hash_size_; i++) {
//      if (kmer_hash_[i] != NULL)
//        free(kmer_hash_[i]);
//      kmer_hash_[i] = NULL;
//    }
//    free(kmer_hash_);
//  }
  if (kmer_hash_)
    free(kmer_hash_);
  kmer_hash_ = NULL;

  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  all_kmers_size_ = 0;

//  kmer_hash_.clear();
  reference_starting_pos_.clear();
  reference_lengths_.clear();
  data_length_ = 0;
  data_length_forward_ = 0;
  data_ptr_ = 0;
  num_sequences_ = 0;
  kmer_hash_size_ = 0;
  kmer_hash_last_key_ = 0;
  kmer_hash_last_key_initialized_ = false;

  kmer_hash_counts_size_ = 0;
}

int64_t IndexHash::GenerateHashKey(const int8_t *seed, uint64_t seed_length) const {
  int64_t ret = 0;

//  std::string temp = "";

  int current_base_MSB = (seed_length - 1) * 2;

  for (uint64_t current_base = 0; current_base < seed_length;
      current_base++) {
    if (!kIsBase[seed[current_base]])
      return -1;

//    temp += ((char) seed[current_base]);
//    ret |= (kBaseToBwa[seed[current_base]] << (current_base * 2));
    ret |= (kBaseToBwa[seed[current_base]] << (current_base_MSB));
    current_base_MSB -= 2;

//    new_data[current_base >> 2] |=
//        kBaseToBwa[(int32_t) data_[current_base]]
//            << ((3 - (current_base & 3)) << 1);
  }
//
//
//  printf ("hash_key = %X, seed = %s\n", ret, temp.c_str());
//  printf ("%s", temp.c_str());

  return ret;
}

int64_t IndexHash::UpdateHashKey(const int8_t *seed, uint64_t seed_length, int64_t hash_key) const {
  int64_t ret = 0;

  uint64_t current_base_num = seed_length - 1;
  int8_t current_base_2bit = kBaseToBwa[seed[current_base_num]];
  if (current_base_2bit > 3)
    return -1;

//  ret = (hash_key >> 2) | (kBaseToBwa[seed[current_base]] << (current_base * 2));
  ret = (hash_key << 2) | current_base_2bit;

  uint64_t bit_mask = (uint64_t) 0;
  bit_mask = (1 << (seed_length * 2)) - 1;      // Only the seed bits are equal to 1, others are 0.
  // Seed length should never be more than 31, or we'll jump into negative numbers.
  ret &= bit_mask;

  return ret;
}

void IndexHash::CountKmers(const SingleSequence &sequence, int k, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) const { // std::vector<int64_t> &ret_kmer_counts) {
  int64_t hash_key = -1;

//  ret_kmer_counts.resize(std::pow(2, (2 * k)));
  int64_t num_kmers = std::pow(2, (2*k));
  int64_t *kmer_counts = (int64_t *) calloc(sizeof(int64_t), num_kmers);

  for (uint64_t i=0; i<(sequence.get_sequence_length() - k + 1); i++) {
    const int8_t *seed_start = &(sequence.get_data()[i]);

    if (i == 0) {
      hash_key = GenerateHashKey(seed_start, k);
    } else {
      hash_key = UpdateHashKey(seed_start, k, hash_key);
    }

    if (hash_key < 0)
      continue;
    kmer_counts[hash_key] += 1;
  }

  *ret_kmer_counts = kmer_counts;
  *ret_num_kmers = num_kmers;
}

int IndexHash::GenerateFromSingleSequenceOnlyForward(const SingleSequence& sequence) {
  Clear();

//  std::vector<int64_t> kmer_counts;
//  CountKmers(sequence, k_, kmer_counts);
  CountKmers(sequence, k_, &kmer_hash_counts_, &kmer_hash_counts_size_);
  int64_t *kmer_countdown = (int64_t *) malloc(sizeof(int64_t) * kmer_hash_counts_size_);
  memmove(kmer_countdown, kmer_hash_counts_, sizeof(int64_t) * kmer_hash_counts_size_);

//  kmer_counts.resize(std::pow(2, (2 * k_)));
//  for (uint64_t i=0; i<(sequence.get_sequence_length() - k_ + 1); i++) {
//    int8_t *seed_start = &(sequence.get_data()[i]);
//    int64_t hash_key = GenerateHashKey(seed_start, k_);
//    if (hash_key < 0)
//      continue;
//    kmer_counts[hash_key] += 1;
//  }

  if (data_)
    delete[] data_;
  data_ = NULL;
  uint64_t mem_to_alloc = (sequence.get_data_length() + 1);
  data_ = new int8_t[mem_to_alloc];
  if (data_ == NULL) {
//    ErrorReporting::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, ErrorReporting::GetInstance().GenerateErrorMessage(ERR_MEMORY, "Offending variable: data_."));
    return 1;
  }
  data_length_ = mem_to_alloc;
  data_ptr_ = 0;
  data_length_forward_ = data_length_ + 1;
  num_sequences_forward_ = 1;
  memmove(data_, sequence.get_data(), data_length_);



//  int k = 6;
//
//  k_ = k;

//  kmer_hash_.resize(std::pow(2, (2 * k_)));
//  kmer_hash_size_ = kmer_hash_.size();
  kmer_hash_size_ = std::pow(2, (2 * k_));
  kmer_hash_ = (int64_t **) malloc(sizeof(int64_t *) * kmer_hash_size_);


  all_kmers_size_ = sequence.get_data_length() - k_ + 1;
  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  all_kmers_ = (int64_t *) calloc(sizeof(int64_t), all_kmers_size_);



  int64_t kmer_hash_ptr = 0;
  for (int64_t i=0; i<kmer_hash_size_; i++) {
//    kmer_hash_[i] = (int64_t *) calloc(sizeof(int64_t), kmer_hash_counts_[i]);
    if (kmer_hash_counts_[i] > 0)
      kmer_hash_[i] = (all_kmers_ + kmer_hash_ptr);
    else
      kmer_hash_[i] = NULL;

    kmer_hash_ptr += kmer_hash_counts_[i];
    kmer_countdown[i] -= 1;
  }

  kmer_hash_last_key_ = 0;
  kmer_hash_last_key_initialized_ = false;

//  for (uint64_t i=0; i<kmer_hash_size_; i++) {
//    kmer_hash_[i].resize(kmer_hash_counts_[i]);
//    kmer_countdown[i] -= 1;
//  }

  int64_t hash_key = -1;
  bool last_skipped = true;

  for (uint64_t i=0; i<(sequence.get_sequence_length() - k_ + 1); i++) {
    int8_t *seed_start = &(data_[i]);
//    int64_t hash_key = GenerateHashKey(seed_start, k_);

    if (i == 0 || last_skipped == true) {
      hash_key = GenerateHashKey(seed_start, k_);
      last_skipped = false;
    } else {
      hash_key = UpdateHashKey(seed_start, k_, hash_key);
    }

    if (hash_key < 0) {
      last_skipped = true;
      continue;
    }

//    kmer_hash_[hash_key].push_back(((int64_t) i));
//    printf ("\nkmer_counts[%ld] = %ld\n", hash_key, kmer_counts[hash_key]);
//    fflush(stdout);
    kmer_hash_[hash_key][kmer_countdown[hash_key]] = ((int64_t) i);
    kmer_countdown[hash_key] -= 1;

//    for (int j=0; j<kmer_hash_[hash_key].size(); j++)
//    {
//      printf (", ");
//      printf ("%ld", kmer_hash_[hash_key][j]);
//    }
//    printf ("\n");
  }

  if (kmer_countdown)
    free(kmer_countdown);
  kmer_countdown = NULL;

//  for (uint64_t i=0; i<kmer_hash_.size(); i++) {
//    std::sort(kmer_hash_[i].begin(), kmer_hash_[i].end(), std::greater<int64_t>());
//  }

  return 0;
}

int IndexHash::FindAllRawPositionsOfSeed(int8_t* seed, uint64_t seed_length,
                                         uint64_t max_num_of_hits,
                                         int64_t** hits,
                                         uint64_t* start_hit,
                                         uint64_t* num_hits) const {

  int64_t hash_key = GenerateHashKey(seed, seed_length);

  if (hash_key < 0 || hash_key >= kmer_hash_size_) {
    return 3;
  }

  *hits = &(kmer_hash_[hash_key][0]);
  *start_hit = 0;
  *num_hits = kmer_hash_counts_[hash_key];

  if (kmer_hash_counts_[hash_key] == 0)
    return 1;

  if (kmer_hash_counts_[hash_key] > max_num_of_hits)
    return 2;

  return 0;
}

int IndexHash::FindAllRawPositionsOfIncrementalSeed(int8_t* seed, uint64_t seed_length,
                                         uint64_t max_num_of_hits,
                                         int64_t** hits,
                                         uint64_t* start_hit,
                                         uint64_t* num_hits) {



//  int64_t hash_key = GenerateHashKey(seed, seed_length);

  int64_t hash_key = 0;

  if (kmer_hash_last_key_initialized_ == false) {
    hash_key = GenerateHashKey(seed, seed_length);
    kmer_hash_last_key_initialized_ = true;
//    printf ("GenerateHashKey: hash_key = %ld\n", hash_key);
//    printf ("GenerateHashKey: kmer_hash_last_key_ = %ld\n", kmer_hash_last_key_);
//    printf ("seed = %s\n", GetSubstring((char *) seed, seed_length).c_str());
//    printf ("kmer_hash_size_ = %ld\n", kmer_hash_size_);
//    fflush(stdout);
  } else {
    hash_key = UpdateHashKey(seed, seed_length, kmer_hash_last_key_);
//    printf ("UpdateHashKey: hash_key = %ld\n", hash_key);
//    printf ("UpdateHashKey: kmer_hash_last_key_ = %ld\n", kmer_hash_last_key_);
//    printf ("seed = %s\n", GetSubstring((char *) seed, seed_length).c_str());
//    fflush(stdout);
//    exit(1);
  }

  if (hash_key < 0 || hash_key >= kmer_hash_size_) {
    kmer_hash_last_key_initialized_ = false;
    return 3;
  }

  kmer_hash_last_key_ = hash_key;
  int64_t kmer_hash_count = kmer_hash_counts_[hash_key];

  *hits = &(kmer_hash_[hash_key][0]);
  *start_hit = 0;
  *num_hits = kmer_hash_count;

  if (kmer_hash_count == 0)
    return 1;

  if (kmer_hash_count > max_num_of_hits)
    return 2;

  return 0;
}

//int IndexHash::FindAllRawPositionsOfSeedHash(int8_t* seed, uint64_t seed_length,
//                                         uint64_t max_num_of_hits,
//                                         HashIndexVector** hits,
//                                         uint64_t* start_hit,
//                                         uint64_t* num_hits, int64_t *kmer_hash_key) {
//  int64_t hash_key = *kmer_hash_key;
//  if (hash_key < 0)
//    hash_key = GenerateHashKey(seed, seed_length);
//  else
//    hash_key = UpdateHashKey(seed, seed_length, hash_key);
//  *kmer_hash_key = hash_key;
//
////  for (uint64_t current_base = 0; current_base < seed_length;
////      current_base++) {
////    printf ("%c", ((char) seed[current_base]));
////  }
////  printf ("\n");
//
//  if (hash_key < 0) {
//    *hits = NULL;
//    *start_hit = 0;
//    *num_hits = 0;
//    return 3;
//  }
//
//  *hits = &(kmer_hash_[hash_key]);
//  *start_hit = 0;
////  printf ("num_hits = %ld, hash_key = %ld (%X), kmer_hash_.size() = %ld\n", kmer_hash_[hash_key].size(), hash_key, hash_key, kmer_hash_.size());
//  *num_hits = kmer_hash_[hash_key].size();
//
////  printf ("hits.size() = %ld\n", (*hits)->size());
//
//  if (kmer_hash_[hash_key].size() == 0)
//    return 1;
//
//  if (kmer_hash_[hash_key].size() > max_num_of_hits)
//    return 2;
//
////  printf ("num_hits = %ld\n", *num_hits);
////  printf ("Good!\n");
//
//  return 0;
//}

const int8_t* IndexHash::get_data() const {
  return data_;
}

//void IndexHash::set_data(int8_t* data) {
//  data_ = data;
//}

uint64_t IndexHash::get_data_length() const {
  return data_length_;
}

//void IndexHash::set_data_length(uint64_t dataLength) {
//  data_length_ = dataLength;
//}

uint64_t IndexHash::get_data_length_forward() const {
  return data_length_forward_;
}

//void IndexHash::set_data_length_forward(uint64_t dataLengthForward) {
//  data_length_forward_ = dataLengthForward;
//}
//
//uint64_t IndexHash::get_data_ptr() const {
//  return data_ptr_;
//}
//
//void IndexHash::set_data_ptr(uint64_t dataPtr) {
//  data_ptr_ = dataPtr;
//}

const std::vector<std::string>& IndexHash::get_headers() const {
  return headers_;
}

//void IndexHash::set_headers(const std::vector<std::string>& headers) {
//  headers_ = headers;
//}

int IndexHash::get_k() const {
  return k_;
}

void IndexHash::set_k(int k) {
  k_ = k;
}

//const std::vector<HashIndexVector>& IndexHash::get_kmer_hash() const {
//  return kmer_hash_;
//}
//
//void IndexHash::set_kmer_hash(const std::vector<HashIndexVector>& kmerHash) {
//  kmer_hash_ = kmerHash;
//}

uint64_t IndexHash::get_num_sequences() const {
  return num_sequences_;
}

//void IndexHash::set_num_sequences(uint64_t numSequences) {
//  num_sequences_ = numSequences;
//}

uint64_t IndexHash::get_num_sequences_forward() const {
  return num_sequences_forward_;
}

//void IndexHash::set_num_sequences_forward(uint64_t numSequencesForward) {
//  num_sequences_forward_ = numSequencesForward;
//}
//
//uint64_t IndexHash::get_num_threads() const {
//  return num_threads_;
//}
//
//void IndexHash::set_num_threads(uint64_t numThreads) {
//  num_threads_ = numThreads;
//}

const std::vector<uint64_t>& IndexHash::get_reference_lengths() const {
  return reference_lengths_;
}

std::string IndexHash::VerboseToString() const {
  return (std::string(""));
}

//void IndexHash::set_reference_lengths(
//    const std::vector<uint64_t>& referenceLengths) {
//  reference_lengths_ = referenceLengths;
//}

const std::vector<uint64_t>& IndexHash::get_reference_starting_pos() const {
  return reference_starting_pos_;
}

//void IndexHash::set_reference_starting_pos(
//    const std::vector<uint64_t>& referenceStartingPos) {
//  reference_starting_pos_ = referenceStartingPos;
//}

//int64_t*& IndexHash::get_suffix_array() {
//  return suffix_array_;
//}

//void IndexHash::InitNumThreads(uint32_t num_threads) {
//}

int IndexHash::LoadFromFile(std::string index_path) {
  return 0;
}

int IndexHash::GenerateFromFile(std::string sequence_file_path) {
  return 0;
}

int IndexHash::GenerateFromSequenceFile(SequenceFile& sequence_file) {
  return 0;
}

int IndexHash::GenerateFromSingleSequence(SingleSequence& sequence) {
  return 0;
}

int IndexHash::LoadOrGenerate(std::string reference_path,
                              std::string out_index_path, bool verbose) {
  return 0;
}

int IndexHash::StoreToFile(std::string output_index_path) {
  return 0;
}

int64_t IndexHash::RawPositionConverter(
    int64_t raw_position, int64_t query_length, int64_t* ret_absolute_position,
    int64_t* ret_relative_position, SeqOrientation* ret_orientation,
    int64_t* ret_reference_index_with_reverse) const {
  return 0;
}

int64_t IndexHash::RawPositionToReferenceIndexWithReverse(int64_t raw_position) const {
  return 0;
}

//int64_t IndexHash::RawToAbsolutePositionConverter(int64_t raw_position) {
//  return raw_position;
//}

void IndexHash::Verbose(FILE* fp) const {
}

int IndexHash::SerializeIndex_(FILE *fp_out) {
  return 0;
}

int IndexHash::DeserializeIndex_(FILE *fp_in) {
  return 0;
}

int IndexHash::IsManualCleanupRequired(std::string function_name) const {
  return 1;
}

int IndexHash::CreateIndex_(int8_t* data, uint64_t data_length) {
  return 1;
}

//void IndexHash::set_suffix_array(int64_t*& suffixArray) {
//  suffix_array_ = suffixArray;
//}
