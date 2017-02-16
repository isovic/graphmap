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

//  std::string OverlapMHAPVerbose(OwlerData* owler_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters, int64_t ref_id, int64_t hits_start, std::vector<int> &lcskpp_indices);
//
//  OverlapResult GenerateOverlapResult(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed, std::string ref_header, int64_t read_id, int64_t read_length, bool read_reversed, std::string read_header);

 private:
  std::shared_ptr<SequenceFile> ref_;
  std::shared_ptr<SequenceFile> reads_;
//  std::vector<std::shared_ptr<is::MinimizerIndex>> indexes_;
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

  void GenerateOutput_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, OwlerData &owler_data);

  inline int32_t HitPosRef_(const uint128_t& hit);
  inline int32_t HitPosRead_(const uint128_t& hit);
  inline int32_t HitDiag_(const uint128_t& hit);
  inline int32_t HitSeqId_(const uint128_t& hit);
  void FindRegionBounds_(const std::vector<uint128_t> &all_hits, int64_t begin_hit, int64_t end_hit, int64_t &start_read, int64_t &end_read, int64_t &start, int64_t &end);
  void AppendSeedHits_(const uint128_t& seed, std::shared_ptr<is::MinimizerIndex> index, bool threshold_hits, double count_cutoff, bool is_overlapper, int64_t qid, std::vector<uint128_t> &all_hits);
  int LCSkFilter_(const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int64_t qid, int64_t tid, int32_t seed_len, PairwiseOverlap &overlap);
  void LCSk_(std::vector<uint128_t> &events, int64_t n, int64_t k, std::vector<uint64_t> &matches_starts, std::vector<uint64_t> &matches_indices, std::vector<int32_t> &lcsk_indices, int64_t &lcsk_len);
  inline uint128_t MakeHit_(const uint128_t& seq_id, const uint128_t& diag, const uint128_t& pos_ref, const uint128_t& pos_read);

  int PrepareEvents_(const std::vector<uint128_t> &hits, int64_t begin_hit, int64_t end_hit, int64_t seed_len,
                            std::vector<uint128_t> &events, std::vector<uint64_t> &matches_starts, std::vector<uint64_t> &matches_indices, int64_t &max_seq_len);
  std::string GenerateMHAPLine(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap);
  int CheckOverlap_(std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters, const PairwiseOverlap& overlap);

};

#endif /* OWLER_H_ */
