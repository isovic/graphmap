/*
 * owler_data.h
 *
 *  Created on: Jul 2, 2015
 *      Author: isovic
 */

#ifndef OWLER_DATA_H_
#define OWLER_DATA_H_

#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <sstream>
#include <vector>
#include "containers/range.h"
#include "sequences/single_sequence.h"
#include "index/index.h"

class SeedHit {
 public:
  SeedHit() {
    query_pos = 0;
    ref_pos = 0;
    seed_type = 0;
  }

  SeedHit(uint32_t qpos, uint32_t rpos, uint32_t stype) {
    query_pos = qpos;
    ref_pos = rpos;
    seed_type = stype;
  }

  uint32_t query_pos;
  uint32_t ref_pos;
  uint8_t seed_type;  /// ID of an enumerated seed. Since many gapped qgrams can be used, this specifies which one has been utilized.
};

class SeedHit2 {
 public:
  SeedHit2() {
    query_pos = 0;
    ref_pos = 0;
    ref_id = 0;
  }

  SeedHit2(uint32_t qpos, uint32_t rpos, uint32_t rid) {
    query_pos = qpos;
    ref_pos = rpos;
    ref_id = rid;
  }

  std::string VerboseToString() {
    std::stringstream ss;
    ss << "ref_id = " << ref_id << ", q_pos = " << query_pos << ", r_pos = " << ref_pos;
    return ss.str();
  }

  uint32_t query_pos;
  uint32_t ref_pos;
  uint32_t ref_id;
};

/// Holds all relevant info on overlaps between two sequences.
class PairwiseOverlapData {
 public:
  PairwiseOverlapData();
  ~PairwiseOverlapData();
  void Init(uint32_t query_id, int64_t query_len, uint32_t ref_id, int64_t ref_len);
  int GetOverlap(Range &overlap_ref, Range &overlap_query, bool &is_ref_reverse);
  void SetOverlap(const Range &overlap_ref, const Range &overlap_query, const bool &is_ref_reverse);
  std::string VerboseToString();

  std::vector<SeedHit> seed_hits;
  uint32_t last_update;
  int64_t num_unique_hits;
  int32_t viability;  /// If 0, the overlap is considered unviable for an overlap.

  uint32_t query_id_;
  int64_t query_len_;
  uint32_t ref_id_;
  int64_t ref_len_;

 private:
  Range overlap_ref_;
  Range overlap_query_;
  bool is_ref_reverse_;
  bool is_overlap_initialized_;
};

class SingleOverlap {
 public:
  SingleOverlap() {
    overlap_ref.start = overlap_ref.end = 0;
    overlap_query.start = overlap_query.end = 0;
    is_ref_reverse = false;
    num_kmers = 0;
  }

  SingleOverlap(int64_t qstart, int64_t qend, int64_t rstart, int64_t rend, bool rreverse, int64_t n_kmers) {
    overlap_ref.start = rstart;
    overlap_ref.end = rend;
    overlap_query.start = qstart;
    overlap_query.end = qend;
    is_ref_reverse = rreverse;
    num_kmers = n_kmers;
  }

  ~SingleOverlap() {
  }

  Range overlap_ref;
  Range overlap_query;
  bool is_ref_reverse;
  int64_t num_kmers;
};

class OwlerData {
 public:
  OwlerData();
  ~OwlerData();
  void Init(SingleSequence *read, std::vector<Index*> &indexes);

  std::vector<PairwiseOverlapData> overlaps; /// Vector is the size of number of reads in the input dataset (the number of sequences in the reference file).
  std::vector<std::string> seed_types;  /// All instances of gapped qgrams used for lookup, enumerated.
  std::string unmapped_reason;

  std::vector<SeedHit2> seed_hits2;
  std::vector<int64_t> num_unique_hits;
  std::vector<int64_t> last_update;

//  std::vector<SingleOverlap> final_overlaps;
  std::string overlap_lines;

  int64_t num_seeds_over_limit;
  int64_t num_seeds_with_no_hits;
  int64_t num_seeds_errors;

 private:
  SingleSequence *read_;
  std::vector<Index*> *indexes_;
};

class OverlapResult {
 public:
  int64_t read_id;
  int64_t ref_id;
  float jaccard_score;
  int64_t shared_minmers;
  bool read_is_reverse;   /// In the MHAP-like output, read is always considered to be forward oriented.
  int64_t read_start;
  int64_t read_end;
  int64_t read_length;
  bool ref_is_reverse;
  int64_t ref_start;
  int64_t ref_end;
  int64_t ref_length;

  int64_t front_id;
  int64_t back_id;

  std::string read_header;
  std::string ref_header;
  int64_t cov_bases_read;
  int64_t cov_bases_ref;

  OverlapResult();
  std::string GenerateMHAPLine();
  std::string GenerateAFGLine();
  std::string GeneratePAFLine();
};

struct overlapresult_sort_key
{
    inline bool operator() (const OverlapResult& op1, const OverlapResult& op2) {
      if (op1.read_id == op2.read_id) {
        if (op1.ref_id == op2.ref_id) {
          return (op1.jaccard_score > op2.jaccard_score);
        } else {
          return (op1.ref_id < op2.ref_id);
        }
      } else {
        return (op1.read_id < op2.read_id);
      }
      return false;
    }
};

struct overlaps_greater_than_key
{
    inline bool operator() (const PairwiseOverlapData& op1, const PairwiseOverlapData& op2) {
      if (op1.num_unique_hits > op2.num_unique_hits)
        return true;
      return false;
    }
};

struct seedhits2_refid_less_than_key
{
    inline bool operator() (const SeedHit2& op1, const SeedHit2& op2) {
      if (op1.ref_id < op2.ref_id)
        return true;
      else if (op1.ref_id == op2.ref_id && op1.ref_pos == op2.ref_pos) {
        return op1.query_pos < op2.query_pos;
      } else if (op1.ref_id == op2.ref_id) {
        return op1.ref_pos < op2.ref_pos;
      }

      return false;
    }
};

#endif /* OWLER_DATA_H_ */
