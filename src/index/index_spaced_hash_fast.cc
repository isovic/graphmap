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

CompiledSeed::CompiledSeed(std::string m_shape) {
  Generate(m_shape);
}

void CompiledSeed::Generate(std::string m_shape) {
  int32_t mask_shift = 0, zero_start = 0, zero_end = 0;
  uint64_t mask = 0;

  shape = m_shape;
  mask_ops.clear();

  for (int32_t i = 0; i <= shape.size(); i++) {
    if (i == shape.size() || (i > 0 && shape[i] == '0' && shape[i - 1] == '1')) {
      mask_shift = (zero_end - zero_start) * 2;
      if (mask != 0) {
        mask_ops.push_back(MaskOperation(mask, mask_shift));
      }

      mask = 0;
      zero_start = i;
    } else if (i > 0 && shape[i] == '1' && shape[i - 1] == '0') {
      zero_end = i;
    }

    if (shape[i] == '1') {
      mask |= (0x0000000000000003 << (i * 2));
    }
  }
}

uint64_t CompiledSeed::Apply(uint64_t full_seed) {
  uint64_t ret = 0x0;
  for (int32_t i = 0; i < mask_ops.size(); i++) {
    ret |= ((full_seed & mask_ops[i].mask) >> (mask_ops[i].shift));
  }
  return ret;
}

IndexSpacedHashFast::IndexSpacedHashFast() {
  data_ = NULL;
  kmer_hash_array_ = NULL;
  kmer_counts_ = NULL;
  all_kmers_ = NULL;
  shape_index_ = NULL;

  Clear();

  InitShapesPredefined(SHAPE_TYPE_444);
//  InitShapesPredefined(SHAPE_TYPE_66);

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

  num_kmers_ = 0;
  all_kmers_size_ = 0;
}

int64_t IndexSpacedHashFast::GenerateHashKeyFromShape(int8_t *seed, const char *shape, int64_t shape_length) const {
  uint64_t ret = 0;
  uint64_t current_accepted_base = 0;

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

  for (int64_t i = 0; i < shape_length; i++) {
    if (shape[i] != '1')
      continue;
    num_accepted_bases += 1;
  }

  int64_t num_sparse_kmers = (int64_t) ceil(std::pow(2.0f, (2.0f * num_accepted_bases)));

  return num_sparse_kmers;
}

int IndexSpacedHashFast::FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const {

  seed_length = shape_index_length_;
  *ret_hits = NULL;
  *ret_start_hit = 0;
  *ret_num_hits = 0;

  int64_t current_data_ptr = 0;
  int64_t *all_hits = NULL;

  for (int64_t i = 0; i < shapes_lookup_.size(); i++) {
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
    LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "1Offending variable: reference_index. Values: reference_index = %ld, raw_position = %ld, data_length = %ld.", reference_index, raw_position, data_length_));
    return reference_index;
  }

  int64_t relative_pos = raw_position >> 32;  // (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
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

int IndexSpacedHashFast::FindAllRawPositionsOfSeedKey(int64_t hash_key, int64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_hits, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const {
  int64_t *all_hits = NULL;
  int64_t num_hits = 0;

  if (hash_key >= 0 && hash_key < num_kmers_ && kmer_counts_[hash_key] > 0) {
    all_hits = kmer_hash_array_[hash_key];
    num_hits = kmer_counts_[hash_key];
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
  LogSystem::GetInstance().Log(
  VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG,
                                      true, FormatString("Creating spaced hash index.\n"), "CreateIndex_");

  if (kmer_hash_array_)
    free(kmer_hash_array_);
  kmer_hash_array_ = NULL;
  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  if (kmer_counts_)
    free(kmer_counts_);
  kmer_counts_ = NULL;

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_), "CreateIndex_");

  int64_t num_kmers = 0;
  CountKmersFromShape(data_, data_length_, shape_index_, shape_index_length_, &kmer_counts_, &num_kmers);
  int64_t *kmer_countdown = (int64_t *) malloc(sizeof(int64_t) * num_kmers);
  memmove(kmer_countdown, kmer_counts_, sizeof(int64_t) * num_kmers);
  num_kmers_ = num_kmers;

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Kmer counting finished (kmer_counts.size() = %ld)\n", num_kmers_), "CreateIndex_");

  int64_t total_num_kmers = 0;
  for (uint64_t i = 0; i < num_kmers; i++) {
    kmer_countdown[i] = 0;
    total_num_kmers += kmer_counts_[i];
  }

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

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("Index memory allocated.\n"), "CreateIndex_");

  int64_t hash_key = -1;

  uint64_t current_ref_id = 0;

  /// Calculate the largest gapped spaced seed length, so we don't step out of boundaries of the read.
  int64_t k = 0;
  for (int32_t i = 0; i < shape_index_length_; i++) {
    k += ((shape_index_[i] == '1') ? 1 : 2);  /// '0' can also mean an insertion, so it can occupy two bases instead of one.
  }

  for (uint64_t i = 0; i < (data_length_ - k + 1); i++) {
    if (i >= (reference_starting_pos_[current_ref_id] + reference_lengths_[current_ref_id]))
      current_ref_id += 1;

    int8_t *seed_start = &(data_[i]);
    hash_key = GenerateHashKeyFromShape(seed_start, shape_index_, shape_index_length_);

//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("%s = %ld = %X\n", GetSubstring((char *) seed_start, k_).c_str(), hash_key, hash_key), "[]");

    if (hash_key < 0)
      continue;

//    uint64_t local_pos = ((uint64_t) (i - reference_starting_pos_[current_ref_id])) & ((uint64_t) 0x00000000FFFFFFFF);
//    uint64_t ref_id = ((uint64_t) current_ref_id) & ((uint64_t) 0x00000000FFFFFFFF);
//    int64_t coded_position = (int64_t) (local_pos << 32) | ref_id;
//    int64_t reference_index = RawPositionToReferenceIndexWithReverse(i);
    int64_t local_i = ((int64_t) i) - ((int64_t) reference_starting_pos_[current_ref_id]);
    uint64_t local_pos = ((uint64_t) local_i) & ((uint64_t) 0x00000000FFFFFFFF);
    uint64_t ref_id = ((uint64_t) current_ref_id) & ((uint64_t) 0x00000000FFFFFFFF);
//    printf ("%ld\t%ld\t\tcurrent_ref_id = %ld, reference_index = %ld, i = %ld, reference_starting_pos_[reference_index] = %ld\n", ref_id, local_pos, current_ref_id, reference_index, i, reference_starting_pos_[reference_index]);
//    fflush(stdout);

//    uint64_t local_pos = ((uint64_t) (i - reference_starting_pos_[current_ref_id])) & ((uint64_t) 0x00000000FFFFFFFF);
//    uint64_t ref_id = ((uint64_t) current_ref_id) & ((uint64_t) 0x00000000FFFFFFFF);
    int64_t coded_position = (int64_t) ((ref_id << 32) | local_pos);

    kmer_hash_array_[hash_key][kmer_countdown[hash_key]] = coded_position;
    kmer_countdown[hash_key] += 1;
  }

  if (kmer_countdown)
    free(kmer_countdown);
  kmer_countdown = NULL;

  LogSystem::GetInstance().Log(
  VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG,
                                      true, FormatString("Finished creating spaced hash index.\n"), "CreateIndex_");

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

  LogSystem::GetInstance().Log(
  VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG,
                                      true, FormatString("\t- k_...\n"), "DeserializeIndex_");
  if (fread(&shape_index_length_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Error(
    SEVERITY_INT_FATAL,
                                 __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(
                                 ERR_FILE_READ_DATA,
                                                                                             "Occured when reading variable shape_index_length_."));
    return 1;
  }

  shape_index_ = (char *) malloc(sizeof(char) * (shape_index_length_ + 1));
  if (fread(shape_index_, sizeof(char), shape_index_length_, fp_in) != shape_index_length_) {
    LogSystem::GetInstance().Error(
    SEVERITY_INT_FATAL,
                                 __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(
                                 ERR_FILE_READ_DATA,
                                                                                             "Occured when reading variable shape_index_."));
    return 1;
  }
  shape_index_[shape_index_length_] = '\0';

  LogSystem::GetInstance().Log(
  VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG,
                                      true, FormatString("\t- index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_), "DeserializeIndex_");

  LogSystem::GetInstance().Log(
  VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG,
                                      true, FormatString("\t- num_kmers_...\n"), "DeserializeIndex_");
  if (fread(&num_kmers_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Error(
    SEVERITY_INT_FATAL,
                                 __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(
                                 ERR_FILE_READ_DATA,
                                                                                             "Occured when reading variable num_kmers_."));
    return 1;
  }

  if (num_kmers_ <= 0) {
    return 1;
  }

  kmer_counts_ = (int64_t *) malloc(sizeof(int64_t) * num_kmers_);
  if (fread(kmer_counts_, sizeof(int64_t), num_kmers_, fp_in) != num_kmers_) {
    LogSystem::GetInstance().Error(
    SEVERITY_INT_FATAL,
                                 __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(
                                 ERR_FILE_READ_DATA,
                                                                                             "Occured when reading variable kmer_counts_.\n"));
    return 3;
  }

  LogSystem::GetInstance().Log(
  VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG,
                                      true, FormatString("\t- all_kmers_size_...\n"), "DeserializeIndex_");
  if (fread(&all_kmers_size_, sizeof(int64_t), 1, fp_in) != 1) {
    LogSystem::GetInstance().Error(
    SEVERITY_INT_FATAL,
                                 __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(
                                 ERR_FILE_READ_DATA,
                                                                                             "Occured when reading variable all_kmers_size_.\n"));
    return 1;
  }

  all_kmers_ = (int64_t *) malloc(sizeof(int64_t) * all_kmers_size_);
  if (fread(all_kmers_, sizeof(int64_t), all_kmers_size_, fp_in) != all_kmers_size_) {
    LogSystem::GetInstance().Error(
    SEVERITY_INT_FATAL,
                                 __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(
                                 ERR_FILE_READ_DATA,
                                                                                             "Occured when reading variable all_kmers.\n"));
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
    shapes_for_search = {"111111111111", "1111111101111", "11111111001111",
      "1111011111111", "11110111101111", "111101111001111",
      "11110011111111", "111100111101111", "1111001111001111"};
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

  compiled_seeds_.clear();
  for (int32_t i = 0; i < shapes_for_search.size(); i++) {
    compiled_seeds_.push_back(CompiledSeed(shapes_for_search[i]));
  }

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

  compiled_seeds_.clear();
  for (int32_t i = 0; i < shapes_for_search.size(); i++) {
    compiled_seeds_.push_back(CompiledSeed(shapes_for_search[i]));
  }

  return 0;
}

int IndexSpacedHashFast::CalcAllKeysFromSequence(const SingleSequence *read, int64_t kmer_step, std::vector<int64_t> &ret_hash_keys, std::vector<int64_t> &ret_key_counts) {
  int64_t read_length = read->get_sequence_length();

  /// Calculate the largest gapped spaced seed length, so we don't step out of boundaries of the read.
  int64_t k = 0;
  for (int32_t i = 0; i < shape_index_length_; i++) {
    k += ((shape_index_[i] == '1') ? 1 : 2);  /// '0' can also mean an insertion, so it can occupy two bases instead of one.
  }

  /// Calculate all the hash keys.
  std::vector<unsigned __int128> hash_keys_with_pos;
  hash_keys_with_pos.clear();
  hash_keys_with_pos.reserve(read_length * shapes_lookup_.size());
  if (read->get_data_format() == kDataFormatAscii) {
    for (uint64_t i = 0; i < (read_length - k + 1); i += kmer_step) {  // i++) {
      int8_t *seed = (int8_t *) &(read->get_data()[i]);

      for (int64_t j = 0; j < shapes_lookup_.size(); j++) {
        uint64_t hash_key = (uint64_t) GenerateHashKeyFromShape(seed, shapes_lookup_[j].c_str(), shapes_lookup_[j].size());

        if (hash_key < 0 || hash_key >= num_kmers_) {
          continue;
        }

        unsigned __int128 hash_key_with_pos = (((unsigned __int128) hash_key) << 64) | i;
        hash_keys_with_pos.push_back(hash_key_with_pos);
      }
    }
  } else if (read->get_data_format() == kDataFormat2BitPacked2) {
    /// This case not implemented yet. This should be much faster than the above option.
  }

  /// Find unique hash keys.
  std::sort(hash_keys_with_pos.begin(), hash_keys_with_pos.end());
  hash_keys_with_pos.erase(std::unique(hash_keys_with_pos.begin(), hash_keys_with_pos.end()), hash_keys_with_pos.end());

  if (hash_keys_with_pos.size() == 0)
    return 1;

  /// Clear the output containers.
  ret_hash_keys.clear();
  ret_hash_keys.reserve(hash_keys_with_pos.size());
  ret_key_counts.clear();
  ret_key_counts.reserve(hash_keys_with_pos.size());

  /// Two options: a) not filtering multiple seeds in the read but handling each separately, or b) counting the same seeds and returning the number of repeats.
  /// a)
  /// This takes all keys, even the repeats in the read.
  for (int64_t i = 1; i < (hash_keys_with_pos.size()); i++) {
    int64_t hash_key = (int64_t) (hash_keys_with_pos[i] >> 64);
    int64_t pos = (int64_t) (hash_keys_with_pos[i] & MASK_32_BIT);
    ret_hash_keys.push_back(hash_key);
    ret_key_counts.push_back(pos);
  }

//  /// b)
//  /// A key can be repeated multiple times throughout the read (repeats in the sequence + tandem repeats and multiple gapped spaced seeds for same position on the read).
//  /// In this case, truncate all keys into one, and keep the counts to be used for region selection.
//  int64_t hash_key = 0, next_hash_key = 0, last_hash_key = 0;
//  int64_t key_count = 0;
//  hash_key = (int64_t) (hash_keys_with_pos[0] >> 64);
//  for (int64_t i=1; i<(hash_keys_with_pos.size()); i++) {
////    printf ("[%ld] %X %ld\n", i, (int64_t) (ret_coded_hash_keys[i] >> 64), (int64_t) (ret_coded_hash_keys[i] & 0x00000000FFFFFFFF));
////    fflush(stdout);
//    int64_t current_hash_key = (int64_t) (hash_keys_with_pos[i] >> 64);
//    if (current_hash_key == hash_key) {
//      key_count += 1;
//    } else {
//      ret_hash_keys.push_back(hash_key);
//      ret_key_counts.push_back((key_count + 1));
//      key_count = 0;
//      hash_key = current_hash_key;
//    }
//  }
//  /// Handle the last hash key.
//  ret_hash_keys.push_back(hash_key);
//  ret_key_counts.push_back((key_count + 1));

  /// Debugging.
//  for (int64_t i=1; i<(ret_hash_keys.size()); i++) {
//    printf ("[%ld] %X %ld\n", i, ret_hash_keys[i], ret_key_counts[i]);
//    fflush(stdout);
//  }
//  exit(1);

  return 0;
}

//int IndexSpacedHashFast::LookUpHashKeys(int64_t bin_size, const SingleSequence *read, const std::vector<int64_t> &hash_keys, const std::vector<int64_t> &key_counts, std::vector<int64_t> &ret_hits) {
////  std::vector<int64_t> all_hits;
//  ret_hits.reserve(hash_keys.size() * 10);     /// Just a guess at the number of hits, to preallocate the memory.
//
//  float bin_size_inverse = 1.0f / ((float) bin_size);
//
//  for (int64_t i=0; i<hash_keys.size(); i++) {
//    int64_t hash_key = hash_keys[i];
//    int64_t *hits = kmer_hash_array_[hash_key];
//    int64_t num_hits = kmer_counts_[hash_key];
//    if (num_hits == 0)
//      continue;
//    printf ("%X, %ld\n", hash_keys[i], hits[0]);
//    fflush(stdout);
//    int64_t num_hits_before = ret_hits.size();
//    ret_hits.insert(ret_hits.end(), hits, (hits + num_hits));
////    ret_x.insert(ret_x.end(), num_hits, i);
//
//    for (int64_t j=num_hits_before; j<(num_hits_before + num_hits); j++) {
//      ret_hits[j] = (ret_hits[j] & (MASK_64_BIT << 64)) | ((uint64_t) (ret_hits[j] - ((int64_t) (key_counts[i]))));
//    }
//
////    for (int32_t j=0; j<key_counts[i]; j++) {
////      ret_hits.insert(ret_hits.end(), hits, (hits + num_hits));
////    }
//  }
//
//  ret_hits.push_back(0x000000050000000A);
//  int64_t x = 0xF;
//  ret_hits[ret_hits.size()-1] = (ret_hits[ret_hits.size()-1] & (((unsigned __int128) MASK_64_BIT) << 64)) | (((int64_t) (ret_hits[ret_hits.size()-1]&((unsigned __int128) MASK_64_BIT)) - ((int64_t) (x))));
//
////  std::sort(ret_hits.begin(), ret_hits.end());
//
//  printf ("ret_hits.size() = %ld\n", ret_hits.size());
//  fflush(stdout);
//  for (int64_t i=0; i<ret_hits.size(); i++) {
//    printf ("[%ld] ref = %ld, pos = %ld\n", i, int64_t(((uint64_t) ret_hits[i])>>32), (int64_t) (((float) (ret_hits[i]&0x00000000FFFFFFFF)) * bin_size_inverse));
//    fflush(stdout);
//  }
//
//  return 0;
//}

int IndexSpacedHashFast::LookUpHashKeys(int64_t bin_size, const SingleSequence *read, const std::vector<int64_t> &hash_keys, const std::vector<int64_t> &key_counts, std::vector<SeedHit3> &ret_hits) {
//  std::vector<int64_t> all_hits;
  ret_hits.reserve(ret_hits.size() + hash_keys.size() * 10);     /// Just a guess at the number of hits, to preallocate the memory.

  float bin_size_inverse = 1.0f / ((float) bin_size);

  for (int64_t i = 0; i < hash_keys.size(); i++) {
    int64_t x = key_counts[i];
    int64_t hash_key = hash_keys[i];
    int64_t *hits = kmer_hash_array_[hash_key];
    int64_t num_hits = kmer_counts_[hash_key];
    if (num_hits == 0)
      continue;

//    printf ("hash_key = %X, ref_id = %ld, y = %ld, x = %ld\n", hash_keys[i], hits[0]>>32, hits[0]&MASK_32_BIT, x);
//    fflush(stdout);

    for (int64_t j = 0; j < num_hits; j++) {
      SeedHit3 seed_hit;
      seed_hit.ref_id = (int32_t) (((uint64_t) hits[j]) >> 32);
      seed_hit.y = (int32_t) (hits[j] & 0x00000000FFFFFFFF);
      seed_hit.x = (int32_t) x;
      ret_hits.push_back(seed_hit);
    }

//    int64_t num_hits_before = ret_hits.size();
//    ret_hits.insert(ret_hits.end(), hits, (hits + num_hits));
////    ret_x.insert(ret_x.end(), num_hits, i);
//
//    for (int64_t j=num_hits_before; j<(num_hits_before + num_hits); j++) {
//      ret_hits[j] = (ret_hits[j] & (MASK_64_BIT << 64)) | ((uint64_t) (ret_hits[j] - ((int64_t) (key_counts[i]))));
//    }
//
////    for (int32_t j=0; j<key_counts[i]; j++) {
////      ret_hits.insert(ret_hits.end(), hits, (hits + num_hits));
////    }
  }

//  ret_hits.push_back(0x000000050000000A);
//  int64_t x = 0xF;
//  ret_hits[ret_hits.size()-1] = (ret_hits[ret_hits.size()-1] & (((unsigned __int128) MASK_64_BIT) << 64)) | (((int64_t) (ret_hits[ret_hits.size()-1]&((unsigned __int128) MASK_64_BIT)) - ((int64_t) (x))));

//  std::sort(ret_hits.begin(), ret_hits.end(), seed_hit3_compare());

//  printf ("ret_hits.size() = %ld\n", ret_hits.size());
//  fflush(stdout);
//  for (int64_t i=0; i<ret_hits.size(); i++) {
////    printf ("[%ld] ref = %ld, pos = %ld\n", i, int64_t(((uint64_t) ret_hits[i])>>32), (int64_t) (((float) (ret_hits[i]&0x00000000FFFFFFFF)) * bin_size_inverse));
//    printf ("[%ld] ref_id = %ld, y = %ld, x = %ld\n", i, ret_hits[i].ref_id, ret_hits[i].y, ret_hits[i].x);
//    fflush(stdout);
//  }

  return 0;
}
