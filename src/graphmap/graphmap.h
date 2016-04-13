/*
 * treemap_se.h
 *
 *  Created on: Jul 19, 2014
 *      Author: ivan
 */

#ifndef TREEMAP_SE_H_
#define TREEMAP_SE_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "index/index.h"
#include "index/index_sa.h"
#include "index/index_hash.h"
#include "index/index_spaced_hash.h"
#include "index/index_spaced_hash_fast.h"
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

class GraphMap {
 public:
  GraphMap();
  ~GraphMap();

  // Main function for running the mapping process. It generates/loads the index, and handles batch loading of sequences from the reads file.
  void Run(ProgramParameters &parameters);

  // Generates or loads the index of the reference genome.
  int BuildIndex(ProgramParameters &parameters);

  // Loads reads from a file in batches of given size (in MiB), or all at once.
  void ProcessReadsFromSingleFile(ProgramParameters &parameters, FILE *fp_out);

  // Process the loaded batch of reads. Uses OpenMP to do it in parallel. Calls ProcessOneRead for each read in the SequenceFile.
  int ProcessSequenceFileInParallel(ProgramParameters *parameters, SequenceFile *reads, clock_t *last_time, FILE *fp_out, int64_t *ret_num_mapped, int64_t *ret_num_unmapped);

  // Processes a single read from the batch of loaded reads.
  int ProcessRead(MappingData *mapping_data, const std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params);

  // Collects alignments from the given mapping_data and converts them into an appropriate output format (string).
  int CollectAlignments(const SingleSequence *read, const ProgramParameters *parameters, MappingData *mapping_data, std::string &ret_aln_lines);



 private:
  std::vector<Index *> indexes_;

  // Opens the output SAM file for writing if the path is specified. If the path is empty, then output is set to STDOUT.
  FILE* OpenOutSAMFile_(std::string out_sam_path="");
  // Formats the SAM header from a given index.
  std::string GenerateSAMHeader_(ProgramParameters &parameters, Index *index);
  // Generates a default SAM line for unmapped reads.
  std::string GenerateUnmappedSamLine_(MappingData *mapping_data, int64_t verbose_sam_output, const SingleSequence *read) const;

  // Calculates the LCSk of the anchors using the Fenwick tree.
  void CalcLCSFromLocalScoresCacheFriendly_(const Vertices *vertices, bool use_l1_filtering, int64_t l, int64_t allowed_dist, int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices);

  // Count gapped spaced seed hits to regions on the reference.
  // Three different implementations providing the same interface.
  int RegionSelectionNoCopy_(int64_t bin_size, MappingData *mapping_data, const std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int RegionSelectionNoBins_(int64_t bin_size, MappingData *mapping_data, const std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int RegionSelectionNoCopyWithDensehash_(int64_t bin_size, MappingData *mapping_data, const std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters);

  int GraphMap_(ScoreRegistry *local_score, Index *index_read, MappingData *mapping_data, const std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int ProcessKmerCacheFriendly_(int8_t *kmer, int64_t kmer_start_position, ScoreRegistry *local_score, MappingData* mapping_data, Index *index_read, const SingleSequence* read, const ProgramParameters* parameters);

  // Perform the LCSk calculation and simple filtering of the anchores that survived the LCSk.
  int SemiglobalPostProcessRegionWithLCS_(ScoreRegistry *local_score, MappingData *mapping_data, const std::vector<Index *> indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int AnchoredPostProcessRegionWithLCS_(ScoreRegistry *local_score, MappingData *mapping_data, const std::vector<Index *> &indexes, const SingleSequence *read, const ProgramParameters *parameters);

  //
  int EvaluateMappings_(MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int GenerateAlignments_(MappingData *mapping_data, const Index *index, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params);

  // Helper functions.
  int CalculateL1ParametersWithMaximumDeviation_(ScoreRegistry *local_score, std::vector<int> &lcskpp_indices, float maximum_allowed_deviation, int64_t *ret_k, int64_t *ret_l, float *ret_sigma_L2, float *ret_confidence_L1);
  int64_t CountBinsWithinThreshold_(const MappingData *mapping_data, float threshold);
  Region CalcRegionFromBin_(int64_t sorted_bins_index, const MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int CheckRegionSearchFinished_(int64_t current_region, float min_allowed_bin_value, float threshold_step, float *bin_value_threshold, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int CollectFinalMappingsAndMapQ_(bool generate_final_mapping_ptrs, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int CheckMinimumMappingConditions_(MappingResults *mapping_data, L1Results *l1_data, const Index *index, const SingleSequence *read, const ProgramParameters *parameters);

  // Debug.
  int VerboseLocalScoresToFile(std::string file_path, const SingleSequence *read, const ScoreRegistry *local_score, const std::vector<int> *indices, int64_t l_median, float maximum_allowed_deviation, bool check_median_filtering, std::vector<int32_t> *cluster_ids=NULL);

};

#endif /* TREEMAP_SE_H_ */
