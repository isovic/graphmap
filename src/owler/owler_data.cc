/*
 * owler_data.cc
 *
 *  Created on: Jul 2, 2015
 *      Author: isovic
 */

#include "owler_data.h"

PairwiseOverlapData::PairwiseOverlapData() {
  Init(0, 0, 0, 0);
}

PairwiseOverlapData::~PairwiseOverlapData() {
  seed_hits.clear();
}

void PairwiseOverlapData::Init(uint32_t query_id, int64_t query_len, uint32_t ref_id, int64_t ref_len) {
  seed_hits.clear();

  overlap_ref_.start = overlap_ref_.end = 0;
  overlap_query_.start = overlap_query_.end = 0;
  is_ref_reverse_ = false;
  is_overlap_initialized_ = false;
  query_id_ = query_id;
  query_len_ = query_len;
  ref_id_ = ref_id;
  ref_len_ = ref_len;

  viability = 0;

  num_unique_hits = 0;
  last_update = 0;
}

int PairwiseOverlapData::GetOverlap(Range &overlap_ref, Range &overlap_query, bool &is_ref_reverse) {
  if (is_overlap_initialized_ == false)
    return 1;
  overlap_ref = overlap_ref_;
  overlap_query = overlap_query_;
  is_ref_reverse = is_ref_reverse_;
  return 0;
}

void PairwiseOverlapData::SetOverlap(const Range &overlap_ref, const Range &overlap_query, const bool &is_ref_reverse) {
  overlap_ref_ = overlap_ref;
  overlap_query_ = overlap_query;
  is_ref_reverse_ = is_ref_reverse;
  is_overlap_initialized_ = true;
}

std::string PairwiseOverlapData::VerboseToString() {
  std::stringstream ss;

  ss << "query_id = " << query_id_ << ", qlen = " << query_len_ << ", ref_id = " << ref_id_ << ", rlen = " << ref_len_ << ", num_unique_hits = " << num_unique_hits << ", num_seed_hits = " << seed_hits.size() << ", viability = " << viability << ", last_update = " << last_update;

  return ss.str();
}

OwlerData::OwlerData() {
  num_seeds_over_limit = 0;
  num_seeds_with_no_hits = 0;
  num_seeds_errors = 0;
  read_ = NULL;
  indexes_ = NULL;
  seed_hits2.clear();
//  final_overlaps.clear();
  overlap_lines = "";
}

OwlerData::~OwlerData() {
  overlaps.clear();
  seed_types.clear();
  seed_hits2.clear();
//  final_overlaps.clear();
}

void OwlerData::Init(SingleSequence *read, std::vector<Index*> &indexes) {
  read_ = read;
  indexes_ = &indexes;
  overlaps.clear();
  num_unique_hits.resize((indexes[0]->get_num_sequences_forward()*2));      /// Take the reverse sequences into account too.
  last_update.resize((indexes[0]->get_num_sequences_forward()*2));      /// Take the reverse sequences into account too.
//  final_overlaps.clear();
  overlap_lines = "";

//  overlaps.resize((indexes[0]->get_num_sequences_forward()*2));      /// Take the reverse sequences into account too.
//  for (int64_t i = 0; i < overlaps.size(); i++) {
//    overlaps[i].Init(read->get_sequence_id(), read->get_sequence_length(), i, indexes[0]->get_reference_lengths()[i]);
////    overlaps[i].seed_hits.reserve(read->get_sequence_length());
//    overlaps[i].seed_hits.reserve(10);
//  }
  seed_types.clear();
  for (int64_t i = 0; i < indexes.size(); i++) {
    seed_types.push_back(((IndexSpacedHash *) indexes[i])->get_shape_index());
  }
}
