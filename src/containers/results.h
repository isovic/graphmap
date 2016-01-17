/*
 * results.h
 *
 *  Created on: Jan 16, 2016
 *      Author: isovic
 */

#ifndef SRC_CONTAINERS_RESULTS_H_
#define SRC_CONTAINERS_RESULTS_H_



typedef struct Cluster {
  Range query;
  Range ref;
} Cluster;

typedef struct MappingResults {
  int64_t lcs_length = 0;
  int64_t cov_bases_max = 0;
  int64_t cov_bases_query = 0;
  int64_t cov_bases_ref = 0;
  int64_t num_covering_kmers = 0;
  float deviation = 0.0f;
  Range query_coords;
  Range ref_coords;
  bool is_mapped = false;
  bool is_reverse = false;
  int64_t local_score_id = 0;
  std::vector<Cluster> clusters;

//  int64_t num_same_mappings = 0;      // How many mapping positions have exactly the same score.
} MappingResults;

typedef struct L1Results {
  int64_t l1_l = 0;
  double l1_k = 1.0f;
  int64_t l1_lmin = 0;
  int64_t l1_lmax = 0;
  double l1_confidence_abs = 0;
  double l1_std = 0;
  int64_t l1_rough_start = 0;
  int64_t l1_rough_end = 0;
} L1Results;

typedef struct AlignmentResults {
  bool is_aligned = false;
  bool is_reverse = false;
  int64_t pos_start = 0;
  int64_t pos_end = 0;
  std::string cigar = "*";
  std::vector<int8_t> alignment;
  int64_t edit_distance = 0;
  int64_t alignment_score = 0;
  int64_t mapping_quality = 0;
  double evalue = 0.0f;
  int64_t num_secondary_alns = 0;      // How many mapping positions have exactly the same score.

  int64_t num_eq_ops = 0;
  int64_t num_x_ops = 0;
  int64_t num_i_ops = 0;
  int64_t num_d_ops = 0;
  int64_t nonclipped_length = 0;

//  std::string unmapped_reason = "Not processed.";
//
//  double stats_time_region_selection = 0.0;
//  double stats_time_mapping = 0.0;
//  double stats_time_alignment = 0.0;

} AlignmentResults;



//class AlignmentResult {
// public:
//  std::vector<int8_t> query_seq;
//  std::vector<int8_t> target_seq;
//  int64_t aln_pos_start;
//  int64_t aln_pos_end;

//  int64_t edit_dist;
//  int64_t aln_score;
//  double evalue;
//
//  int64_t query_id;
//  std::string query_header;
//  int64_t ref_id;
//  std::string ref_header;
//};


typedef struct MappingMetadata {
  std::string unmapped_reason = "Not processed.";

  double stats_time_region_selection = 0.0;
  double stats_time_mapping = 0.0;
  double stats_time_alignment = 0.0;

} MappingMetadata;

#endif /* SRC_CONTAINERS_RESULTS_H_ */
