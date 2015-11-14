/*
 * owler.h
 *
 *  Created on: Jul 2, 2015
 *      Author: isovic
 */

#ifndef OWLER_H_
#define OWLER_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "index/index.h"
#include "index/index_sa.h"
#include "index/index_hash.h"
#include "index/index_spaced_hash.h"
#include "sequences/single_sequence.h"
#include "sequences/sequence_file.h"
#include "containers/score_registry.h"
#include "utility/program_parameters.h"
#include "utility/utility_general.h"
#include "alignment/cigargen.h"
#include "alignment/local_realignment.h"
#include "containers/region.h"
#include "containers/mapping_data.h"
#include "utility/evalue.h"
#include "containers/vertices.h"

#include "owler/owler_data.h"

//#include "index/index_spaced_hash_fast.h"
#include "index/index_owler.h"



class Owler {
 public:
  Owler();
  ~Owler();

  // Main function for running the mapping process. It generates/loads the index, and handles batch loading of sequences from the reads file.
  void Run(ProgramParameters &parameters);
  // Generates or loads the index of the reference genome.
  int BuildIndex(ProgramParameters &parameters);
  // Loads reads from a file in batches of given size (in MiB), or all at once.
  void ProcessReadsFromSingleFile(ProgramParameters &parameters, FILE *fp_out);
  // Process the loaded batch of reads. Uses OpenMP to do it in parallel. Calls ProcessOneRead for each read in the SequenceFile.
  int ProcessSequenceFileInParallel(ProgramParameters *parameters, SequenceFile *reads, clock_t *last_time, FILE *fp_out, int64_t *ret_num_mapped, int64_t *ret_num_unmapped);

  int ProcessRead(OwlerData *owler_data, std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params);
//  int CollectSAMLines(std::string &ret_sam_lines, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);

  int CollectAMOSLines(std::string &ret_amos_lines, OwlerData *owler_data, const SingleSequence *read, const ProgramParameters *parameters);

  int CollectSeedHits(OwlerData *owler_data, std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int FilterUnlikelyOverlaps(OwlerData *owler_data, std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int ApplyLCS(OwlerData *owler_data, std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);

  int CollectSeedHitsExperimental(OwlerData *owler_data, std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int CollectSeedHitsExperimentalCalcSubseedsFast(OwlerData *owler_data, std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int CollectSeedHitsExperimentalSubseededIndex(OwlerData *owler_data, std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int ApplyLCS2(OwlerData *owler_data, std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int CollectMHAPLines(std::string &ret_overlap_lines, OwlerData *owler_data, const SingleSequence *read, const ProgramParameters *parameters);
  std::string OverlapToMHAP(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed, int64_t read_id, int64_t read_length);
  std::string OverlapToAFG(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed, int64_t read_id, int64_t read_length);
  std::string OverlapToDot(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed, int64_t read_id, int64_t read_length);
  bool CheckDistanceTooBig(OwlerData* owler_data, int64_t index_last, int64_t index_current, float error_rate);
  bool CheckDistanceStep(OwlerData* owler_data, int64_t index_first, int64_t index_last, int64_t index_current, float max_diff);
  int FilterAnchorBreakpoints(const std::vector<int> &lcskpp_indices, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length, int64_t min_cluster_length, float min_cluster_coverage, OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters, std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids=NULL);
  int FilterAnchorBreakpointsExperimental(const std::vector<int> &lcskpp_indices, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length, int64_t min_cluster_length, float min_cluster_coverage, OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters, std::vector<int> &ret_filtered_lcskpp_indices, std::vector<int32_t> *ret_cluster_ids=NULL);

  int OverlapLength(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t *A_start, int64_t *A_end, int64_t *ret_query_length, int64_t *B_start, int64_t *B_end, int64_t *ret_ref_length);
  std::string OverlapMHAPVerbose(OwlerData* owler_data, std::vector<Index*> &indexes, const SingleSequence* read, const ProgramParameters* parameters, int64_t ref_id, int64_t hits_start, std::vector<int> &lcskpp_indices);

  OverlapResult GenerateOverlapResult(std::vector<SeedHit2> &seed_hits, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t ref_id, int64_t reference_length, bool ref_reversed, std::string ref_header, int64_t read_id, int64_t read_length, bool read_reversed, std::string read_header);

 private:
  SequenceFile *reference_;

//  Index *index_;
//  Index *index_secondary_;
  std::vector<Index *> indexes_;

  // Retrieves a file list from the given folder.
  bool GetFileList_(std::string folder, std::vector<std::string> &ret_files);
  // Check if string ends with the given suffix (parameter 'ending'), and returns true if so.
  bool StringEndsWith_(std::string const &full_string, std::string const &ending);
  // Returns only files with one of the following extensions: fasta, fastq, fa, fq, sam.
  void FilterFileList_(std::vector<std::string> &files, std::vector<std::string> &ret_read_files, std::vector<std::string> &ret_sam_files);
  // Opens the output SAM file for writing if the path is specified. If the path is empty, then output is set to STDOUT.
  FILE* OpenOutFile_(std::string out_sam_path="");

  void ClearIndexes_();

  std::string GenerateSAMHeader_(ProgramParameters &parameters, Index *index);

  void CalcLCSFromLocalScoresCacheFriendly_(OwlerData* owler_data, int64_t overlap_id, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices);

  void CalcLCSFromLocalScoresCacheFriendly2_(OwlerData* owler_data, int64_t ref_hits_start, int64_t ref_hits_end, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices);

  int CalcCoveredBases(std::vector<SeedHit2> &seed_hits, int64_t seed_length, std::vector<int> &lcskpp_indices, int64_t hits_start, int64_t hits_end, int64_t *ret_cov_A, int64_t *ret_cov_B);

  // Calculates the LCSk of the anchors using the Fenwick tree.
//  void CalcLCSFromLocalScoresCacheFriendly_(const Vertices *vertices, bool use_l1_filtering, int64_t l, int64_t allowed_dist, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices);
//
//  int HoughCollect_(int64_t bin_size, OwlerData *owler_data, const Index *index, const Index *index_secondary, const SingleSequence *read, const ProgramParameters *parameters);
//
//  int ProcessKmerCacheFriendly_(int8_t *kmer, int64_t kmer_start_position, ScoreRegistry *local_score, MappingData* mapping_data, Index *index_read, const Index* index, const Index* index_secondary, const SingleSequence* read, const ProgramParameters* parameters);
//  int PostProcessRegionWithLCS_(ScoreRegistry *local_score, MappingData *mapping_data, const Index *index, const Index *index_secondary, const SingleSequence *read, const ProgramParameters *parameters);
//  int CalculateL1ParametersWithMaximumDeviation_(ScoreRegistry *local_score, std::vector<int> &lcskpp_indices, float maximum_allowed_deviation, int64_t *ret_k, int64_t *ret_l, float *ret_sigma_L2, float *ret_confidence_L1);
//  int64_t CountBinsWithinThreshold_(const MappingData *mapping_data, float threshold);
//  Region CalcRegionFromBin_(int64_t sorted_bins_index, const MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
//  int CheckRegionSearchFinished_(int64_t current_region, float min_allowed_bin_value, float threshold_step, float *bin_value_threshold, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
//  int EvaluateMappings_(bool evaluate_edit_distance, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
//  int GenerateAlignments_(MappingData *mapping_data, const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params);
//  int CollectFinalMappingsAndMapQ_(bool generate_final_mapping_ptrs, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
//  int CheckMinimumMappingConditions_(InfoMapping *mapping_data, InfoL1 *l1_data, const Index *index, const SingleSequence *read, const ProgramParameters *parameters);
//  int VerboseLocalScoresToFile(std::string file_path, const SingleSequence *read, const ScoreRegistry *local_score, const std::vector<int> *indices, int64_t l_median, float maximum_allowed_deviation, bool check_median_filtering);
//
//  int ExperimentalPostProcessRegionWithLCS_(ScoreRegistry *local_score, MappingData *mapping_data, const Index *index, const Index *index_secondary, const SingleSequence *read, const ProgramParameters *parameters);
//  int ExperimentalPostProcessRegionWithLCS1_(ScoreRegistry *local_score, MappingData *mapping_data, const Index *index, const Index *index_secondary, const SingleSequence *read, const ProgramParameters *parameters);
};

#endif /* OWLER_H_ */
