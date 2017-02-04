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
#include "log_system/log_system.h"
#include "utility/utility_general.h"

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
  is_transcriptome_ = false;

  Clear();

  InitShapesPredefined(SHAPE_TYPE_444);
//  InitShapesPredefined(SHAPE_TYPE_66);
}

IndexSpacedHashFast::IndexSpacedHashFast(uint32_t shape_type) {
  data_ = NULL;
  kmer_hash_array_ = NULL;
  kmer_counts_ = NULL;
  all_kmers_ = NULL;
  shape_index_ = NULL;
  is_transcriptome_ = false;

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

//  genome_id_to_trans_id_.clear();
//  trans_id_to_exons_.clear();
//  trans_id_to_regions_.clear();
//  is_transcriptome_ = false;
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

void IndexSpacedHashFast::CountKmersFromShape(int8_t *sequence_data, int64_t sequence_length, const char *shape, int64_t shape_length, int64_t shape_max_width, int64_t **ret_kmer_counts, int64_t *ret_num_kmers) const {  // std::vector<int64_t> &ret_kmer_counts) {
  int64_t hash_key = -1;

//  ret_kmer_counts.resize(std::pow(2, (2 * k)));
  int64_t num_kmers = CalcNumHashKeysFromShape(shape, shape_length);
  int64_t *kmer_counts = (int64_t *) calloc(sizeof(int64_t), num_kmers);

  for (uint64_t i = 0; i < (sequence_length - shape_max_width + 1); i++) {
    int8_t *seed_start = &(sequence_data[i]);
//    printf ("shape = '%s', shape_length = %ld, shape_max_width = %ld, i = %ld, sequence_length = %ld\n", shape, shape_length, shape_max_width, i, sequence_length);
//    fflush(stdout);
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

  double keygen_time = 0.0f;
  clock_t time_diff1 = clock();
  int64_t num_hits = 0;
  for (int64_t i = 0; i < shapes_lookup_.size(); i++) {
    clock_t test_keygen = clock();
    int64_t hash_key = GenerateHashKeyFromShape(seed, shapes_lookup_[i].c_str(), shapes_lookup_[i].size());
    keygen_time += ((double) clock() - test_keygen) / CLOCKS_PER_SEC;

    if (hash_key < 0 || hash_key > num_kmers_ || kmer_counts_[hash_key] <= 0) {
      continue;
    }
    num_hits += kmer_counts_[hash_key];
  }
  all_hits = (int64_t *) malloc(sizeof(int64_t) * num_hits);
//  printf ("num_hits = %ld\n", num_hits);
//  printf ("keygen_time = %f\n", keygen_time);
//  printf ("total_time = %f\n", ((double) clock() - time_diff1) / CLOCKS_PER_SEC);

//
//  if (num_hits == 134221872) {
//    printf ("WTF!?!\n");
//    printf ("Seed: '%s'\n", GetSubstring((char *) seed, seed_length).c_str());
//  }
//  fflush(stdout);

  for (int64_t i = 0; i < shapes_lookup_.size(); i++) {
    int64_t hash_key = GenerateHashKeyFromShape(seed, shapes_lookup_[i].c_str(), shapes_lookup_[i].size());
    if (hash_key < 0 || hash_key > num_kmers_ || kmer_counts_[hash_key] <= 0) {
      continue;
    }

//    if (hash_key >= 0 && hash_key < num_kmers_ && kmer_counts_[hash_key] > 0) {

//      if (all_hits == NULL)
//        all_hits = (int64_t *) malloc(sizeof(int64_t) * (current_data_ptr + kmer_counts_[hash_key]));
//      else
//        all_hits = (int64_t *) realloc(all_hits, (sizeof(int64_t) * (current_data_ptr + kmer_counts_[hash_key])));
      memmove(&(all_hits[current_data_ptr]), &(kmer_hash_array_[hash_key][0]), kmer_counts_[hash_key] * sizeof(int64_t));
      current_data_ptr += kmer_counts_[hash_key];
//    }
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

int IndexSpacedHashFast::FindAllRawPositionsOfSeedNoCopy(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, std::vector<int64_t *> &ret_hits, std::vector<uint64_t> &ret_num_hits) const {
  seed_length = shape_index_length_;
  ret_hits.clear();
  ret_num_hits.clear();
//  ret_hits.resize(shapes_lookup_.size(), NULL);
//  ret_num_hits.resize(shapes_lookup_.size(), 0);

//  printf ("seed = '%s'\n", GetSubstring((char *) seed, seed_length).c_str());
//  fflush(stdout);

  int64_t total_num_hits = 0;
  for (int64_t i = 0; i < shapes_lookup_.size(); i++) {
    int64_t hash_key = GenerateHashKeyFromShape(seed, shapes_lookup_[i].c_str(), shapes_lookup_[i].size());
    if (hash_key < 0 || hash_key > num_kmers_ || kmer_counts_[hash_key] <= 0) {
      continue;
    }

    ret_hits.push_back(&kmer_hash_array_[hash_key][0]);
    ret_num_hits.push_back(kmer_counts_[hash_key]);
    total_num_hits += ret_num_hits.back();
  }

  if (total_num_hits == 0) {
    return 1;
  }

  if ((max_num_of_hits > 0 && total_num_hits > ((int64_t) max_num_of_hits))) {
    return 2;
  }

  return 0;
}

int64_t IndexSpacedHashFast::RawPositionConverter2(int64_t raw_position, int64_t query_length, int64_t *ret_absolute_position, int64_t *ret_relative_position, SeqOrientation *ret_orientation, int64_t *ret_reference_index_with_reverse) const {
//  if (raw_position < 0 || raw_position >= data_length_)
//    return -2;

  int64_t reference_index = (int64_t) (raw_position & MASK_REF_ID);
  if (reference_index < 0) {
    ERROR_REPORT(ERR_UNEXPECTED_VALUE, "1Offending variable: reference_index. Values: reference_index = %ld, raw_position = %ld, data_length = %ld.", reference_index, raw_position, data_length_);
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

void IndexSpacedHashFast::CalcPercentileHits(double percentile, int64_t *ret_count, int64_t *ret_max_seed_count) {
  return CalcPercentileHits_(kmer_counts_, num_kmers_, percentile, ret_count, ret_max_seed_count);
}

void IndexSpacedHashFast::CalcPercentileHits_(int64_t *seed_counts, int64_t num_seeds, double percentile, int64_t *ret_count, int64_t *ret_max_seed_count) {
  // Sort the counts so that we can get an occurance histogram.
  // We will select a percentile of the data as the cutoff value.
  LOG_DEBUG("Sorting the seed counts to obtain the occurrence histogram.\n");
  std::vector<size_t> counts_indices;
  OrderedSortArray(seed_counts, num_seeds, counts_indices);
  int64_t num_nonzero = 0;
  for (int64_t i=0; i<num_seeds; i++) {
    if (seed_counts[counts_indices[i]] != 0) {
      num_nonzero = num_seeds - i;
      break;
    }
  }

  int64_t percentil_id = num_seeds - (int64_t) round(((double) 1.0 - percentile) * ((double) num_nonzero)) - 2; // The + 2 is for the case when the percentile is so low that it would amount to id = 0. This way we at least cut out two most occurring seeds, often AAAAA* and TTTTT*.
  if (percentil_id < 0) { percentil_id = 0; }
  if (percentil_id > num_seeds) { percentil_id = num_seeds - 1; }
  int64_t percentile_val = seed_counts[counts_indices[percentil_id]];

  *ret_count = percentile_val;

  if (ret_max_seed_count) {
    *ret_max_seed_count = seed_counts[counts_indices.back()];
  }
}

int IndexSpacedHashFast::CreateIndex_(int8_t *data, uint64_t data_length) {
  LOG_DEBUG_MEDHIGH("Creating spaced hash index.\n");

  if (kmer_hash_array_)
    free(kmer_hash_array_);
  kmer_hash_array_ = NULL;
  if (all_kmers_)
    free(all_kmers_);
  all_kmers_ = NULL;
  if (kmer_counts_)
    free(kmer_counts_);
  kmer_counts_ = NULL;

  LOG_DEBUG_MEDHIGH("Index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_);

  int64_t num_kmers = 0;
  CountKmersFromShape(data_, data_length_, shape_index_, shape_index_length_, shape_max_width_, &kmer_counts_, &num_kmers);
  int64_t *kmer_countdown = (int64_t *) malloc(sizeof(int64_t) * num_kmers);
  memmove(kmer_countdown, kmer_counts_, sizeof(int64_t) * num_kmers);
  num_kmers_ = num_kmers;

  LOG_DEBUG_MEDHIGH("Kmer counting finished (kmer_counts.size() = %ld)\n", num_kmers_);

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

  LOG_DEBUG_MEDHIGH("Index memory allocated.\n");

  int64_t hash_key = -1;

  uint64_t current_ref_id = 0;

  /// Calculate the largest gapped spaced seed length, so we don't step out of boundaries of the read.
  int64_t k = shape_max_width_;

  for (uint64_t i = 0; i < (data_length_ - k + 1); i++) {
    if (i >= (reference_starting_pos_[current_ref_id] + reference_lengths_[current_ref_id]))
      current_ref_id += 1;

    int8_t *seed_start = &(data_[i]);
    hash_key = GenerateHashKeyFromShape(seed_start, shape_index_, shape_index_length_);

    if (hash_key < 0)
      continue;

    int64_t local_i = ((int64_t) i) - ((int64_t) reference_starting_pos_[current_ref_id]);
    uint64_t local_pos = ((uint64_t) local_i) & ((uint64_t) 0x00000000FFFFFFFF);
    uint64_t ref_id = ((uint64_t) current_ref_id) & ((uint64_t) 0x00000000FFFFFFFF);
    int64_t coded_position = (int64_t) ((ref_id << 32) | local_pos);

//    LOG_DEBUG("%s = %ld = %X, local_pos = %lu\n", GetSubstring((char *) seed_start, k).c_str(), hash_key, hash_key, local_pos);
    if (hash_key >= num_kmers) {
      LOG_ALL("ERROR: hash_key >= num_kmers! hash_key = %ld, num_kmers = %ld\n", hash_key, num_kmers);
      LOG_ALL("  %s = %ld = %X, local_pos = %lu, coded_position = %ld\n", GetSubstring((char *) seed_start, k).c_str(), hash_key, hash_key, local_pos, coded_position);
    }
    if (kmer_countdown[hash_key] < 0 ||  kmer_countdown[hash_key] >= kmer_counts_[hash_key]) {
      LOG_ALL("ERROR: kmer_countdown[hash_key] is wrong! kmer_countdown[hash_key] = %ld, kmer_counts_[hash_key] = %ld\n", kmer_countdown[hash_key], kmer_counts_[hash_key]);
      LOG_ALL("  %s = %ld = %X, local_pos = %lu, coded_position = %ld\n", GetSubstring((char *) seed_start, k).c_str(), hash_key, hash_key, local_pos, coded_position);
    }
    kmer_hash_array_[hash_key][kmer_countdown[hash_key]] = coded_position;
    kmer_countdown[hash_key] += 1;
  }

  if (kmer_countdown)
    free(kmer_countdown);
  kmer_countdown = NULL;

  LOG_DEBUG_MEDHIGH("Finished creating spaced hash index.\n");

  return 0;
}

int IndexSpacedHashFast::SerializeIndex_(FILE* fp_out) {
  int64_t vector_length = 0;

  fwrite(&is_transcriptome_, sizeof(bool), 1, fp_out);

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

int64_t IndexSpacedHashFast::get_shape_max_width() const {
  return shape_max_width_;
}

void IndexSpacedHashFast::set_shape_max_width(int64_t shape_max_width) {
  shape_max_width_ = shape_max_width;
}

int IndexSpacedHashFast::DeserializeIndex_(FILE* fp_in) {
  if (shape_index_)
    free(shape_index_);
  shape_index_ = NULL;
  shape_index_length_ = 0;
  shape_max_width_ = 0;
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

  LOG_DEBUG_MEDHIGH("\t- is_transcriptome_...\n");
  if (fread(&is_transcriptome_, sizeof(bool), 1, fp_in) != 1) {
    FATAL_REPORT(ERR_FILE_READ_DATA, "Occured when reading variable is_transcriptome_.");
    return 1;
  }

  LOG_DEBUG_MEDHIGH("\t- k_...\n");
  if (fread(&shape_index_length_, sizeof(int64_t), 1, fp_in) != 1) {
    FATAL_REPORT(ERR_FILE_READ_DATA, "Occured when reading variable shape_index_length_.");
    return 1;
  }

  shape_index_ = (char *) malloc(sizeof(char) * (shape_index_length_ + 1));
  if (fread(shape_index_, sizeof(char), shape_index_length_, fp_in) != shape_index_length_) {
    FATAL_REPORT(ERR_FILE_READ_DATA, "Occured when reading variable shape_index_.");
    return 1;
  }
  shape_index_[shape_index_length_] = '\0';

  // Calculate the max width (each gap can consume two spaces).
  shape_max_width_ = 0;
  for (int32_t i = 0; i < shape_index_length_; i++) {
    shape_max_width_ += ((shape_index_[i] == '1') ? 1 : 2);  /// '0' can also mean an insertion, so it can occupy two bases instead of one.
  }

  LOG_DEBUG_MEDHIGH("\t- index shape: '%s', length: %ld.\n", shape_index_, shape_index_length_);
  LOG_DEBUG_MEDHIGH("\t- num_kmers_...\n");
  if (fread(&num_kmers_, sizeof(int64_t), 1, fp_in) != 1) {
    FATAL_REPORT(ERR_FILE_READ_DATA, "Occured when reading variable num_kmers_.");
    return 1;
  }

  if (num_kmers_ <= 0) {
    return 1;
  }

  kmer_counts_ = (int64_t *) malloc(sizeof(int64_t) * num_kmers_);
  if (fread(kmer_counts_, sizeof(int64_t), num_kmers_, fp_in) != num_kmers_) {
    FATAL_REPORT(ERR_FILE_READ_DATA, "Occured when reading variable kmer_counts_.\n");
    return 3;
  }

  LOG_DEBUG_MEDHIGH("\t- all_kmers_size_...\n");
  if (fread(&all_kmers_size_, sizeof(int64_t), 1, fp_in) != 1) {
    FATAL_REPORT(ERR_FILE_READ_DATA, "Occured when reading variable all_kmers_size_.\n");
    return 1;
  }

  LOG_DEBUG_MEDHIGH("\t- started allocating space for kmers...\n");
  all_kmers_ = (int64_t *) malloc(sizeof(int64_t) * all_kmers_size_);
  LOG_DEBUG_MEDHIGH("\t- started reading the kmers from file...\n");
  if (fread(all_kmers_, sizeof(int64_t), all_kmers_size_, fp_in) != all_kmers_size_) {
    FATAL_REPORT(ERR_FILE_READ_DATA, "Occured when reading variable all_kmers.\n");
    return 3;
  }

  LOG_DEBUG_MEDHIGH("\t- initializing the kmer_hash_array_...\n");
  kmer_hash_array_ = (int64_t **) malloc(sizeof(int64_t *) * num_kmers_);
  int64_t kmer_ptr = 0;
  for (int64_t i = 0; i < num_kmers_; i++) {
    if (kmer_counts_[i] > 0)
      kmer_hash_array_[i] = (all_kmers_ + kmer_ptr);
    else
      kmer_hash_array_[i] = NULL;
    kmer_ptr += kmer_counts_[i];
  }



//#ifndef RELEASE_VERSION
//  FILE *fp_debug = fopen (FormatString("temp.kmercounts.%s.csv", shape_index_).c_str(), "w");
//  std::vector<size_t> kmer_counts_indices;
//  OrderedSortArray(kmer_counts_, num_kmers_, kmer_counts_indices);
//  for (int64_t i1=(kmer_counts_indices.size()-1); i1>=0; i1--) {
//    fprintf (fp_debug, "%6X\t%ld\n", kmer_counts_indices[i1], kmer_counts_[kmer_counts_indices[i1]]);
//  }
//  fclose(fp_debug);
//#endif

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

  /// Calculate the largest gapped spaced seed length, so we don't step out of boundaries of the read.
  shape_max_width_ = 0;
  for (int32_t i = 0; i < shape_index_length_; i++) {
    shape_max_width_ += ((shape_index_[i] == '1') ? 1 : 2);  /// '0' can also mean an insertion, so it can occupy two bases instead of one.
  }

  return 0;
}

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

  /// Calculate the largest gapped spaced seed length, so we don't step out of boundaries of the read.
  shape_max_width_ = 0;
  for (int32_t i = 0; i < shape_index_length_; i++) {
    shape_max_width_ += ((shape_index_[i] == '1') ? 1 : 2);  /// '0' can also mean an insertion, so it can occupy two bases instead of one.
  }

  return 0;
}




