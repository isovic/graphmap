/*
 * raw_alignment.h
 *
 *  Created on: Nov 12, 2015
 *      Author: isovic
 */

#ifndef SRC_CONTAINERS_RAW_ALIGNMENT_H_
#define SRC_CONTAINERS_RAW_ALIGNMENT_H_

struct RawAlignment {
  int64_t aln_start = 0;
  int64_t aln_end = 0;
  std::vector<int8_t> alignment;
  std::string cigar = "*";
  std::string md = "*"; /// MD field from SAM output.
  SeqOrientation orientation = kForward;
  int64_t ref_id = 0;
  std::string ref_header = "";
  int64_t query_id = 0;
  std::string query_header = "";
  int64_t eq_ops = 0, x_ops = 0, i_ops = 0, d_ops = 0;  /// Counts of CIGAR operations.
  int64_t aligned_len = 0;  /// Number of aligned bases from the read (not counting clipped bases).
  int64_t num_clipped_front = 0;  /// Number of clipped bases at the beginning of the read.
  int64_t num_clipped_back = 0;   /// Number of clipped bases at the end of the read.
};

#endif /* SRC_CONTAINERS_RAW_ALIGNMENT_H_ */
