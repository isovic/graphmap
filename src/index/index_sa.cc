/*
 * index_sa.cc
 *
 *  Created on: Jun 4, 2014
 *      Author: ivan
 */


#include <algorithm>
#include <memory>

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "index/index_sa.h"
#include "log_system/log_system.h"

IndexSA::IndexSA() {
  data_ = NULL;
  suffix_array_ = NULL;
  Clear();
}

IndexSA::~IndexSA() {
  Clear();
}

void IndexSA::Clear() {
  if (data_)
    delete[] data_;
  data_ = NULL;

  if (suffix_array_)
    delete[] suffix_array_;
  suffix_array_ = NULL;

  reference_starting_pos_.clear();
  reference_lengths_.clear();
  data_length_ = 0;
  data_length_forward_ = 0;
  data_ptr_ = 0;
  num_sequences_ = 0;
}

int IndexSA::FindAllRawPositionsOfSeed(int8_t *seed, uint64_t seed_length, uint64_t max_num_of_hits, int64_t **ret_entire_sa, uint64_t *ret_start_hit, uint64_t *ret_num_hits) const { //, std::vector<int64_t> &return_positions) {
  if (suffix_array_ == NULL)
    return 3;

  saidx64_t sa_hit_begin = -1;
  saidx64_t num_hits = sa_search64((sauchar_t *) data_, (saidx64_t) data_length_,
                                   (sauchar_t *) seed, (saidx64_t) seed_length,
                                   (saidx64_t *) suffix_array_, (saidx64_t) data_length_,
                                   (saidx64_t *) &sa_hit_begin);

  *ret_entire_sa = suffix_array_;
  *ret_start_hit = sa_hit_begin;
  *ret_num_hits = num_hits;

  if (num_hits <= 0)
    return 1;

  if ((max_num_of_hits > 0 && num_hits > ((int64_t) max_num_of_hits)))
    return 2;

  return 0;
}

int IndexSA::CreateIndex_(int8_t *data, uint64_t data_length) {
  int64_t *sa = new int64_t[data_length];

  int ret_val = divsufsort64((sauchar_t *) data, (saidx64_t *) sa, data_length);

  if (ret_val) {
    LogSystem::GetInstance().Log(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: ret_val (divsufsort64 return)."));
    return 1;
  }

  if (suffix_array_) {
    delete[] suffix_array_;
  }

  suffix_array_ = sa;

  return 0;
}

int IndexSA::SerializeIndex_(FILE* fp_out) {
  uint64_t vector_length = 0;
  vector_length = data_length_;
  fwrite(&vector_length, sizeof(vector_length), 1, fp_out);
  fwrite(suffix_array_, sizeof(saidx64_t), vector_length, fp_out);

  return 0;
}

int IndexSA::IsManualCleanupRequired(std::string function_name) const {
  return 1;
}

int IndexSA::DeserializeIndex_(FILE* fp_in) {
  uint64_t vector_length = 0;

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t- suffix_array_...\n"), "Deserialize_");
  if (fread(&vector_length, sizeof(vector_length), 1, fp_in) != 1) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable suffix_array_->size()."));
    return 1;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t  vector_length = %ld...\n", vector_length), "Deserialize_");
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t  freeing old suffix_array_...\n", vector_length), "Deserialize_");
  if (suffix_array_)
    delete[] suffix_array_;
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t  allocating new suffix_array_...\n"), "Deserialize_");
  suffix_array_ = new int64_t[vector_length];
  if (suffix_array_ == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when allocating memory for suffix_array_."));
    return 2;
  }

  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t  reading suffix_array_...\n"), "Deserialize_");
  if (fread(suffix_array_, sizeof(saidx64_t), vector_length, fp_in) != vector_length) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FILE_READ_DATA, "Occured when reading variable suffix_array_->data()."));
    return 3;
  }
  LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_MED_DEBUG | VERBOSE_LEVEL_HIGH_DEBUG, true, FormatString("\t  done.\n"), "Deserialize_");

  return 0;
}

void IndexSA::Verbose(FILE* fp) const {
  fprintf (fp, "Num sequences forward: %ld\n", num_sequences_forward_);
  fprintf (fp, "Num sequences: %ld\n", num_sequences_);
  fprintf (fp, "Data length forward: %ld\n", data_length_forward_);
  fprintf (fp, "Data length: %ld\n", data_length_);
  fprintf (fp, "\n");
  for (uint64_t i=0; i<num_sequences_forward_; i++) {
    fprintf (fp, "Header: '%s'\n", headers_[i].c_str());
    fprintf (fp, "Sequence start: %ld\n", reference_starting_pos_[i]);
    fprintf (fp, "Sequence length: %ld\n", reference_lengths_[i]);
  }
}

std::string IndexSA::VerboseToString() const {
  return (std::string(""));
}
