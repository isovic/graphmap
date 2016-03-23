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
  bool is_reverse = false;            // This should be deprecated and replaced with 'orientation'.
  int64_t ref_start = 0;              // Starting position of the alignment on the reference. If orientation == kReverse, this assumes that the read should be reverse complemented and the reference stays fwd. pos_start is adjusted accordingly to denote the starting position of the alignment of the reversed read.
  int64_t ref_end = 0;                // See pos_start. This is the end position of the alignment.
  int64_t query_start = 0;            // Starting position of the alignment on the read. Everything before this position should be clipped.
  int64_t query_end = 0;              // Ending position of the alignment on the read. Everything after this position should be clipped.
  std::string cigar = "*";            // In case orientation == kReverse, 'cigar' contains the reverse of the 'alignment' operations.
  std::string md = "";
  int64_t edit_distance = 0;
  int64_t alignment_score = 0;
  int64_t mapping_quality = 0;
  double evalue = 0.0f;
  int64_t num_secondary_alns = 0;      // How many mapping positions have similar score.

  int64_t raw_pos_start = 0;          // Internally, the fwd read is mapped to a reference and its reverse complement (which have been joined in a single massive sequence). The raw_pos_start then holds the absolute coordinate of the alignment in such joined sequence data.
  int64_t raw_pos_end = 0;            // See raw_pos_start. This is the end position of the alignment in global coordinates.
  std::vector<uint8_t> raw_alignment;     // Hold the alignment in the global coordinate space (between raw_pos_start and raw_pos_end). Cannot be used with pos_start and pos_end in case the read should be reverse complemented. In this case, the alignment needs to be reversed.
  std::vector<uint8_t> alignment;     // Hold the alignment in the local coordinate space (between ref_start and ref_end). If orientation == kForward, alignment == raw_alignment. Otherwise it's the reverse complement.

  SeqOrientation orientation = kForward;
  int64_t ref_id = -1;
  std::string ref_header = "*";
  int64_t ref_len = 0;
  int64_t query_id = -1;
  std::string query_header = "*";
  int64_t query_len = 0;

  int64_t num_eq_ops = 0;
  int64_t num_x_ops = 0;
  int64_t num_i_ops = 0;
  int64_t num_d_ops = 0;
  int64_t nonclipped_length = 0;

//  int8_t *ref_data = NULL;
//  int8_t *read_data = NULL;

  // These are parameters of alignment which were used to produce the results.
  int32_t aln_mode_code = 0;          // Type of alignment which was performed to produce the results stored in this structure.

  int64_t reg_pos_start = 0;          // Local coordinates of the alignment's start and end positions within the region determined by GetRegionData() function.
  int64_t reg_pos_end = 0;            // Local coordinates of the alignment's start and end positions within the region determined by GetRegionData() function.

} AlignmentResults;



typedef struct MappingMetadata {
  std::string unmapped_reason = "Not processed.";

  double time_region_selection = 0.0;
  double time_mapping = 0.0;
  double time_alignment = 0.0;
  double time_region_seed_lookup = 0.0;
  double time_region_hitsort = 0.0;
  double time_region_conversion = 0.0;
  double time_region_alloc = 0.0;
  double time_region_counting = 0.0;



} MappingMetadata;

#endif /* SRC_CONTAINERS_RESULTS_H_ */
