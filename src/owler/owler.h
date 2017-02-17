/*
 * owler.h
 *
 *  Created on: Jul 2, 2015
 *      Author: isovic
 */

#ifndef OWLER_H_
#define OWLER_H_

#include <memory>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "sequences/single_sequence.h"
#include "sequences/sequence_file.h"
#include "containers/score_registry.h"
#include "program_parameters.h"
#include "utility/utility_general.h"
#include "alignment/cigargen.h"
#include "containers/region.h"
#include "containers/mapping_data.h"
#include "utility/evalue.h"
#include "containers/vertices.h"

#include "owler/owler_data.h"

#include "minimizer_index/minimizer_index.h"

#include "utility/tictoc.h"



class Owler {
 public:
  Owler();
  ~Owler();

  // Main function for running the mapping process. It generates/loads the index, and handles batch loading of sequences from the reads file.
  void Run(ProgramParameters &parameters);

 private:
  std::shared_ptr<SequenceFile> ref_;
  std::shared_ptr<SequenceFile> reads_;
  std::shared_ptr<is::MinimizerIndex> index_;

  // Opens the output SAM file for writing if the path is specified. If the path is empty, then output is set to STDOUT.
  FILE* OpenOutFile_(std::string out_sam_path="");

  // Generates or loads the index of the reference genome.
  int BuildIndex_(ProgramParameters &parameters);

  // Process the loaded batch of reads. Uses OpenMP to do it in parallel. Calls ProcessOneRead for each read in the SequenceFile.
  int ProcessSequenceFileInParallel_(ProgramParameters &parameters, std::shared_ptr<SequenceFile> reads, TicToc &tt_all, FILE *fp_out);

  int ProcessRead_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data);

  int CollectHits_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data);

  int ClusterHits_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, int32_t diag_epsilon, OwlerData &owler_data);
  int ClusterHits2_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, int32_t diag_epsilon, OwlerData &owler_data);

  void GenerateOutput_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data);


  void AppendSeedHits_(const uint128_t& seed, std::shared_ptr<is::MinimizerIndex> index, bool threshold_hits, double count_cutoff, bool is_overlapper, int64_t qid, std::vector<uint128_t> &all_hits);

  int LCSkFilter_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int32_t seed_len, PairwiseOverlap &overlap);

  void LCSk_(std::vector<uint128_t> &events, int64_t n, int64_t k, std::vector<uint64_t> &matches_starts, std::vector<uint64_t> &matches_indices, std::vector<int32_t> &lcsk_indices, int64_t &lcsk_len);

  void FilterColinear_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters,
                       const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int64_t seed_len, const std::vector<int32_t> &raw_lcsk_indices,
                       std::vector<int32_t> &lcsk_indices, std::vector<int32_t> *cluster_ids, int32_t &num_sv);

  int PrepareEvents_(const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int64_t seed_len,
                            std::vector<uint128_t> &events, std::vector<uint64_t> &matches_starts, std::vector<uint64_t> &matches_indices, int64_t &max_seq_len);

  bool CheckOverlap_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap);


  std::string GenerateMHAPLine_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap);
  std::string GeneratePAFLine_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap);

  std::string GenerateDebugInfo_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap);
  int64_t CalcEditDist_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const PairwiseOverlap& overlap);
  double CalcRatio_(const PairwiseOverlap& overlap);

  int FilterAnchorBreakpoints_(const std::vector<int32_t> &lcskpp_indices, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length,
                                     int64_t min_cluster_length, float min_cluster_coverage, const std::vector<uint128_t> &hits,
                                     const ProgramParameters* parameters, std::vector<int32_t> &ret_filtered_lcskpp_indices,
                                     std::vector<int32_t> *ret_cluster_ids);
  bool CheckDistanceTooBig_(const std::vector<uint128_t> &hits, int64_t index_last, int64_t index_current, float error_rate);

  void WriteHits_(std::string out_path, const std::vector<uint128_t> &hits, int64_t hits_start, int64_t hits_end,
                 int64_t ref_id, std::string read_header, int64_t read_length,
                 std::string reference_header, int64_t reference_length,
                 const std::vector<int32_t> *indices_to_output, const std::vector<int32_t> *cluster_ids);

  static inline uint128_t MakeHit_(const uint128_t& seq_id, const uint128_t& diag, const uint128_t& pos_ref, const uint128_t& pos_read) {
    return ((seq_id << 96) | (diag << 64) | (pos_ref << 32) | (pos_read));
  }

  static inline int32_t HitPosRead_(const uint128_t& hit) {
    return (int32_t) (hit & kSeedMask32_1);
  }

  static inline int32_t HitPosRef_(const uint128_t& hit) {
    return (int32_t) ((hit & kSeedMask32_2) >> 32);
  }

  static inline int32_t HitDiag_(const uint128_t& hit) {
    return (int32_t) ((hit & kSeedMask32_3) >> 64);
  }

  static inline int32_t HitSeqId_(const uint128_t& hit) {
    return (int32_t) ((hit & kSeedMask32_4) >> 96);
  }
};

#endif /* OWLER_H_ */
