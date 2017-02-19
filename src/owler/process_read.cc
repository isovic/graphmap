/*
 * process_read.cc
 *
 *  Created on: Mar 20, 2015
 *      Author: isovic
 */

#include <limits>
#include <algorithm>

#include "owler/owler.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"

#include "libs/edlib.h"
#include "libs/edlibcigar.h"



int Owler::ProcessRead_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data) {
  LOG_DEBUG_SPEC_NEWLINE;
  LOG_DEBUG_SPEC("Entered function.\n\n");

  // If the read length is too short, call it unmapped.
  if (read->get_sequence_length() < parameters->min_read_len) {
    LOG_DEBUG_SPEC("Read too short.\n");
    std::stringstream ss;
    ss << "Unmapped_5__readlength_too_short" << "__readlength=" << read->get_sequence_length() << "__limit=" << 80;
    owler_data.unmapped_reason += ss.str();
    return 0;
  }

  int32_t diag_epsilon = 500;     // TODO: Parametrize this.

  TicToc tt_collect;
  tt_collect.start();
  CollectHits_(index_, read, parameters, owler_data);
  tt_collect.stop();
  LOG_DEBUG_SPEC("Time collecting hits: %f sec\n", tt_collect.get_secs());

  TicToc tt_cluster;
  tt_cluster.start();
  ClusterHits2_(index, read, parameters, diag_epsilon, owler_data);
  tt_cluster.stop();
  LOG_DEBUG_SPEC("Time clustering hits: %f sec\n", tt_cluster.get_secs());

  GenerateOutput_(index_, read, parameters, owler_data);

//  owler_data->seed_hits2.reserve(500000);
//
//  /// Check if it's a case of self-overlap. In this case, overlap can be performed faster, because the index will already have pre-processed seeds of all reads.
//  if (parameters->reads_path == parameters->reference_path)
//    CollectSeedHitsExperimentalSubseededIndex(owler_data, indexes, read, parameters);
//  else
//    CollectSeedHitsExperimentalCalcSubseedsFast(owler_data, indexes, read, parameters);
//
//  ApplyLCS2(owler_data, indexes, read, parameters);

  // Just verbose.
  LOG_DEBUG_SPEC_NEWLINE;
  LOG_DEBUG_SPEC("Exiting function.\n");

  return 0;
}

int Owler::CollectHits_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data) {
  TicToc tt_collect;

  int64_t read_len = read->get_sequence_length();
  int64_t num_fwd_seqs = index->get_num_sequences_forward();
  bool is_overlapper = (parameters->reference_path == parameters->reads_path);

  owler_data.hits.clear();
  owler_data.hits.reserve(index->avg_seed_occurrence() * read_len);

  std::vector<uint128_t> seeds;       // All seeds for a given index and a given sequence.
  std::vector<int8_t> seed_key_lens;  // Lengths of seed keys. Since multiple lookup keys can be used, each can be of different length.
  std::vector<int8_t> seed_lens;    // Since there are several
  int64_t qid = read->get_sequence_id();

//  index->CollectLookupSeeds(read->get_data(), read->get_quality(), read->get_sequence_length(), 0.0f, false, false, 1, seeds);
//  index->CollectLookupSeeds(read->get_data(), read->get_quality(), read->get_sequence_length(), 0.0f, false, parameters->use_minimizers, parameters->minimizer_window, seeds);
  index->CollectIndexSeeds(read->get_data(), read->get_quality(), read->get_sequence_length(), 0.0f, false, parameters->use_minimizers, parameters->minimizer_window, seeds);

  owler_data.hits.clear();

  for (int64_t j=0; j<seeds.size(); j++) {
    AppendSeedHits_(seeds[j], index, parameters->threshold_hits, index->count_cutoff(), is_overlapper, qid, owler_data.hits);
  }

  std::sort(owler_data.hits.begin(), owler_data.hits.end());

  return 0;
}

void Owler::AppendSeedHits_(const uint128_t& seed, std::shared_ptr<is::MinimizerIndex> index, bool threshold_hits, double count_cutoff, bool is_overlapper, int64_t qid, std::vector<uint128_t> &all_hits) {
  int64_t key = is::MinimizerIndex::seed_key(seed);
  int32_t pos_read_int32 = is::MinimizerIndex::seed_position(seed);
  uint128_t pos_read = pos_read_int32;

  int32_t num_targets_fwd = index->get_num_sequences_forward();

  const uint128_t *found_seeds = NULL;
  int64_t num_found_seeds = 0;

  int lookup_ret = index->KeyLookup(key, &found_seeds, &num_found_seeds);

  // Something went wrong.
  if (lookup_ret) {
    return;
  }

  // Filter seeds with excessive hits if necessary.
  if (threshold_hits &&
      num_found_seeds >= count_cutoff) {
    return;
  }

  all_hits.insert(all_hits.end(), num_found_seeds, 0);

  // Reformat the new hits. The key part is no longer needed,
  // but source position on read is.
  auto *phits = &all_hits[0] + (all_hits.size() - num_found_seeds);
  int64_t k_added = 0;

  for (int64_t k=0; k<num_found_seeds; k++) {
    if (found_seeds[k] == kInvalidSeed) {
      continue;
    }

    int32_t seq_id_int32 = is::MinimizerIndex::seed_seq_id(found_seeds[k]);
    int32_t seq_id_int32_fwd = seq_id_int32 % num_targets_fwd;
    uint128_t seq_id = seq_id_int32;
    if (is_overlapper && seq_id_int32_fwd <= qid) {
      continue;
    }

    int32_t pos_ref_int32 = is::MinimizerIndex::seed_position(found_seeds[k]);
    uint128_t diag = ((uint128_t) (pos_ref_int32 - pos_read_int32)) & kSeedMask32_1;      // Diagonals can be negative values.
    uint128_t pos_ref = pos_ref_int32;

    phits[k_added++] = MakeHit_(seq_id, diag, pos_ref, pos_read);
  }

  all_hits.resize(all_hits.size() - num_found_seeds + k_added);
}
//
//int Owler::ClusterHits_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, int32_t diag_epsilon, OwlerData &owler_data) {
//
//  int64_t max_shape_width = index->get_shape_max_width();
//
//  int64_t begin_hit = 0;
//  int64_t total_num_hits = owler_data.hits.size();
//
//  for (int64_t end_hit=0; end_hit<total_num_hits; end_hit++) {
//    if ((end_hit + 1) == total_num_hits ||
//        HitSeqId_(owler_data.hits[end_hit+1]) != HitSeqId_(owler_data.hits[end_hit]) ||
//        (HitDiag_(owler_data.hits[end_hit+1]) - HitDiag_(owler_data.hits[end_hit])) >= diag_epsilon) {
//
//      Range rhits(begin_hit, end_hit);     // Store the current range of hits.
//      begin_hit = end_hit + 1;                // Update the beginning of the range right away in case there will be any 'continue' statements.
//
//      int64_t num_hits = rhits.end - rhits.start + 1;
//      int64_t cov_bases_estimate = num_hits * max_shape_width;
//
//      if (cov_bases_estimate < 100) { continue; }
//
//      // Add a cluster.
//      int64_t ref_id = HitSeqId_(owler_data.hits[rhits.end]);
//      int64_t ref_start = index->get_reference_starting_pos()[ref_id];
//      int64_t ref_len = index->get_reference_lengths()[ref_id];
//      int64_t read_len = read->get_sequence_length();
//
//      PairwiseOverlap overlap(read->get_sequence_id(), ref_id);
//      int rv_lcsk = LCSkFilter_(index, read, parameters, owler_data.hits, rhits.start, rhits.end, max_shape_width, overlap);
//
//      if (rv_lcsk) {      // Nothing survived the filter.
//        continue;
//      }
//
//      // Do a basic check on the sanity of an overlap.
//      if (CheckOverlapV1b_(index, read, parameters, overlap) == true) {
//        owler_data.overlaps.push_back(overlap);
//      }
//    }
//  }
//
//  std::sort(owler_data.overlaps.begin(), owler_data.overlaps.end(), [](const PairwiseOverlap& o1, const PairwiseOverlap& o2) { return o1.num_seeds > o2.num_seeds; });
//
//  return 0;
//}

int Owler::ClusterHits2_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, int32_t diag_epsilon, OwlerData &owler_data) {

  int64_t max_shape_width = index->get_shape_max_width();

  int64_t begin_hit = 0;
  int64_t total_num_hits = owler_data.hits.size();

  for (int64_t end_hit=0; end_hit<total_num_hits; end_hit++) {
    if ((end_hit + 1) == total_num_hits ||
        HitSeqId_(owler_data.hits[end_hit+1]) != HitSeqId_(owler_data.hits[end_hit])) {

      Range rhits(begin_hit, end_hit);        // Store the current range of hits.

      begin_hit = end_hit + 1;                // Update the beginning of the range right away in case there will be any 'continue' statements.

      int64_t num_hits = rhits.end - rhits.start + 1;
      int64_t cov_bases_estimate = num_hits * max_shape_width;

      if (num_hits < 4) { continue; }

      // Add a cluster.
      int64_t ref_id = HitSeqId_(owler_data.hits[rhits.end]);
      int64_t ref_start = index->get_reference_starting_pos()[ref_id];
      int64_t ref_len = index->get_reference_lengths()[ref_id];
      int64_t read_len = read->get_sequence_length();

      PairwiseOverlap overlap(read->get_sequence_id(), ref_id);
      int rv_lcsk = WrapLCSk_(index, read, parameters, owler_data.hits, rhits.start, rhits.end, max_shape_width, overlap);

      // Set the overlap bounds.
      if (overlap.lcsk_indices.size() > 0) {
        // TODO: FilterAnchorBreakpoints_ inverts the list to a normal, ascending form.
        int64_t back = rhits.start + overlap.lcsk_indices.front();
        int64_t front = rhits.start + overlap.lcsk_indices.back();

        overlap.query.start = HitPosRead_(owler_data.hits[back]);
        overlap.query.end = HitPosRead_(owler_data.hits[front]) + 1; // +1 means that the end is not inclusive. // + seed_len;  TODO: The end does not cover the last seed entirely!
        overlap.target.start = HitPosRef_(owler_data.hits[back]);
        overlap.target.end = HitPosRef_(owler_data.hits[front]) + 1; // + seed_len;
        overlap.num_seeds = overlap.lcsk_indices.size();
        overlap.cov_bases = overlap.lcsk_len;
      }

      if (rv_lcsk) {      // Nothing survived the filter.
        continue;
      }

      // Do a basic check on the sanity of an overlap.
//      if (CheckOverlapV1b_(index, read, parameters, overlap) == true) {
      if (CheckOverlapV3_(index, read, parameters, overlap) == true)
      {
        owler_data.overlaps.push_back(overlap);
      }
    }
  }

  std::sort(owler_data.overlaps.begin(), owler_data.overlaps.end(), [](const PairwiseOverlap& o1, const PairwiseOverlap& o2) { return o1.num_seeds > o2.num_seeds; });

  return 0;
}

void Owler::GenerateOutput_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data) {
  for (int64_t i=0; i<owler_data.overlaps.size(); i++) {
//    std::string overlap_line = GenerateMHAPLine_(index, read, parameters, owler_data.overlaps[i]);
    std::string overlap_line;

    if (parameters->outfmt == "mhap") {
      overlap_line = GenerateMHAPLine_(index, read, parameters, owler_data.overlaps[i]);
    } else if (parameters->outfmt == "paf") {
      overlap_line = GeneratePAFLine_(index, read, parameters, owler_data.overlaps[i]);
    } else {
      LOG_ALL("Unknown output format '%s'. Defaulting to PAF.\n", parameters->outfmt.c_str());
      overlap_line = GeneratePAFLine_(index, read, parameters, owler_data.overlaps[i]);
    }

    std::string debug_info;

    if (parameters->verbose_sam_output > 0) {
      debug_info = "\t-> " + GenerateDebugInfo_(index, read, parameters, owler_data.overlaps[i]);
    }

    owler_data.overlap_lines += overlap_line + debug_info + "\n";
  }
}



std::string Owler::GenerateMHAPLine_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap) {
  std::stringstream ret;

  int64_t read_id = overlap.qid;
  std::string read_header = TrimToFirstSpace(read->get_header());
  int64_t read_length = read->get_sequence_length();
  int64_t read_start = overlap.query.start;
  int64_t read_end = overlap.query.end;

  int64_t ref_id = overlap.tid % index->get_num_sequences_forward();
  std::string ref_header = TrimToFirstSpace(index->get_headers()[overlap.tid]);
  int64_t ref_start = overlap.target.start;
  int64_t ref_end = overlap.target.end;
  int64_t ref_length = index->get_reference_lengths()[overlap.tid];

  bool read_is_reverse = overlap.tid >= index->get_num_sequences_forward();
  bool ref_is_reverse = false;
  if (read_is_reverse) {
    ref_start = ref_length - overlap.target.end;
    ref_end = ref_length - overlap.target.start;
  }

  float jaccard_score = std::min(1.0f, ((float) overlap.cov_bases) / ((float) (overlap.target.end - overlap.target.start)));
  int64_t shared_minmers = overlap.num_seeds;

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

//  ret << "\tcov_bases = " << overlap.cov_bases;
  return ret.str();

//  std::stringstream ss;

//  float jaccard_score = std::min(1.0f, ((float) overlap.cov_bases) / ((float) (overlap.target.end - overlap.target.start)));
//  int64_t shared_minmers = overlap.num_seeds;
//
//  bool read_is_reverse = overlap.tid > index->get_num_sequences_forward();
//  bool ref_is_reverse = false;
//
//  int64_t read_id = overlap.qid;
//  int64_t read_start = overlap.query.start;
//  int64_t read_end = overlap.query.end;
//  int64_t read_length = read->get_sequence_length();
//  if (read_is_reverse) {
//    read_start = read_length - overlap.query.end - 1;
//    read_end = read_length - overlap.query.start;
//  }
//
//  int64_t ref_id = overlap.tid % index->get_num_sequences_forward();
//  int64_t ref_start = overlap.target.start;
//  int64_t ref_end = overlap.target.end;
//  int64_t ref_length = index->get_reference_lengths()[overlap.tid];
//
//  ss << (read_id + 1) << " ";      /// read1_id
//  ss << (ref_id + 1) << " ";      /// read2_id
//  ss << jaccard_score << " ";      /// Jaccard score
//  ss << shared_minmers << " ";        /// Shared minmers
//  ss << (read_is_reverse ? 1 : 0) << " ";  /// A is reverse
//  ss << read_start << " ";
//  ss << read_end << " ";
//  ss << read_length << " ";
//  ss << (ref_is_reverse ? 1 : 0) << " ";
//  ss << ref_start << " ";
//  ss << ref_end << " ";
//  ss << ref_length;
//
////  ss << "\tcov_bases = " << overlap.cov_bases;

//  return ss.str();
}

std::string Owler::GeneratePAFLine_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap) {
  std::stringstream ret;

  int64_t read_id = overlap.qid;
  std::string read_header = TrimToFirstSpace(read->get_header());
  int64_t read_length = read->get_sequence_length();
  int64_t read_start = overlap.query.start;
  int64_t read_end = overlap.query.end;

  int64_t ref_id = overlap.tid % index->get_num_sequences_forward();
  std::string ref_header = TrimToFirstSpace(index->get_headers()[overlap.tid]);
  int64_t ref_start = overlap.target.start;
  int64_t ref_end = overlap.target.end;
  int64_t ref_length = index->get_reference_lengths()[overlap.tid];


//  bool read_is_reverse = overlap.tid > index->get_num_sequences_forward();
//  bool ref_is_reverse = false;
//  if (read_is_reverse) {
//    read_start = read_length - overlap.query.end - 1;
//    read_end = read_length - overlap.query.start;
//  }
//
  bool read_is_reverse = overlap.tid >= index->get_num_sequences_forward();
  bool ref_is_reverse = false;
  if (read_is_reverse) {
    ref_start = ref_length - overlap.target.end;
    ref_end = ref_length - overlap.target.start;
  }

//  int64_t shared_minmers = overlap.num_seeds * parameters->minimizer_window;
  int64_t shared_minmers = overlap.num_seeds; // * parameters->minimizer_window;

  ret << read_header << "\t";
  ret << read_length << "\t";
  ret << read_start + 0 << "\t";      // Zero based coordinates.
  ret << read_end + 0 << "\t";

  ret << (read_is_reverse ? "-" : "+") << "\t";
  ret << ref_header << "\t";
  ret << ref_length << "\t";
  ret << ref_start + 0 << "\t";
  ret << ref_end + 0 << "\t";

//  ret << overlap.cov_bases * parameters->minimizer_window << "\t";
  ret << overlap.cov_bases << "\t";
  ret << (ref_end - ref_start) << "\t";

  ret << "255" << "\t";
//  ret << "cm:i:" << shared_minmers * parameters->minimizer_window;
  ret << "cm:i:" << shared_minmers;

  return ret.str();
}
