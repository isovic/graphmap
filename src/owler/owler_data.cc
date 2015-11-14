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





OverlapResult::OverlapResult() {
  read_id = ref_id = shared_minmers = read_start = read_end = read_length = ref_start = ref_end = ref_length = 0;
  jaccard_score = 0.0f;
  read_is_reverse = ref_is_reverse = false;
  front_id = back_id = 0;
  cov_bases_read = cov_bases_ref = 0;
  read_header = ref_header = "";
}

std::string OverlapResult::GenerateMHAPLine() {
  std::stringstream ret;
  ret << (read_id + 1) << " ";      /// read1_id
  ret << (ref_id + 1) << " ";      /// read2_id
  ret << jaccard_score << " ";      /// Jaccard score
  ret << shared_minmers << " ";        /// Shared minmers
  ret << (read_is_reverse ? 1 : 0) << " ";  /// A is reverse
  ret << read_start << " ";
  ret << read_end << " ";
  ret << read_length << " ";
  ret << (ref_is_reverse ? 1 : 0) << " ";
  ret << ref_start << " ";
  ret << ref_end << " ";
  ret << ref_length;
  return ret.str();
}

std::string OverlapResult::GenerateAFGLine() {
  std::stringstream ret;

//    int64_t front_id = lcskpp_indices.front();
//    int64_t back_id = lcskpp_indices.back();
//
//    int64_t A_start = seed_hits[hits_start + back_id].query_pos;
//    int64_t A_end = seed_hits[hits_start + front_id].query_pos + 12;
//    int64_t B_start = seed_hits[hits_start + back_id].ref_pos;
//    int64_t B_end = seed_hits[hits_start + front_id].ref_pos + 12;
//
//    int64_t read1_id = read_id + 1;
//    int64_t read2_id = ref_id + 1;
//    int64_t score = std::min((A_end - A_start), (B_end - B_start));
//
//    //ahg - Ahang. Length of the non-overlapping portion of the first read.
//    //bhg - Bhang. Length of the non-overlapping portion of the second read.
//    /// The position (alignment_start - clip_count_front) would be roughly where alignment of one reads starts on another.
//    /// If this value is > 0, the first read starts within read2, and thus ahang needs to be negative (hence the '-' sign).
//    int64_t ahang = 0; // - (alignment_start - clip_count_front);
//    int64_t bhang = 0; // reference_length - (alignment_end + clip_count_back);
//

//    std::string adj = (ref_is_reverse == false) ? OVERLAP_NORMAL : OVERLAP_INNIE;
//    ret << "{OVL" << "\n";
//    ret << "adj:" << adj << "\n";
//    ret << "rds:" << (read_id + 1) << "," << (ref_id + 1) << "\n";
//    ret << "scr:" << score << "\n";
//    ret << "ahg:" << ahang << "\n";
//    ret << "bhg:" << bhang << "\n";
//    ret << "}";

  return ret.str();
}

std::string OverlapResult::GeneratePAFLine() {
  std::stringstream ret;
  if (read_header == "") { ret << (read_id + 1) << "\t"; } else { ret << TrimToFirstSpace(read_header) << "\t"; }
  ret << read_length << "\t";
  ret << read_start << "\t";
  ret << read_end << "\t";

  ret << (ref_is_reverse ? "-" : "+") << "\t";  /// Or maybe this should be (read_is_reverse)? In this case, upper and lower parts need to be switched (ref and read).
  if (ref_header == "") { ret << (ref_id + 1) << "\t"; } else { ret << TrimToFirstSpace(ref_header) << "\t"; }
  ret << ref_length << "\t";
  ret << ref_start << "\t";
  ret << ref_end << "\t";

  float idty_read = ((float) cov_bases_read) / ((float) (read_end - read_start));
  float idty_ref = ((float) cov_bases_ref) / ((float) (ref_end - ref_start));
  ret << ((idty_read < idty_ref) ? cov_bases_ref : cov_bases_read) << "\t";
  ret << ((idty_read < idty_ref) ? (ref_end - ref_start) : (read_end - read_start)) << "\t";

  ret << "255" << "\t";
  ret << "cm:i:" << shared_minmers;

  return ret.str();
}
