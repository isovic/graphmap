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
#include <limits.h>

#include "minimizer_index/minimizer_index.h"
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
#include "transcriptome.h"
#include "aligner/sam_parser.h"

struct ExonMatch {
	int start;
	int stop;
	int coverage;
};

class ExonsCluster {
public:
	int cuttedRefStart;
	int cuttedRefEnd;
	std::vector<ExonMatch> exons;
	ExonsCluster(int _cuttedRefStart, int _cuttedRefEnd, std::vector<ExonMatch> _exons) {
		cuttedRefStart = _cuttedRefStart;
		cuttedRefEnd = _cuttedRefEnd;
		exons = _exons;
	}

	ExonsCluster(std::vector<ExonMatch> _exons) {
		exons = _exons;

		cuttedRefStart = 0;
		cuttedRefEnd = 0;
	}

	int avg_covereage() {
		int coverageSum = 0;
		int length = 0;
		for (int var = 0; var < exons.size(); ++var) {
			coverageSum += exons[var].coverage;
			length += exons[var].stop - exons[var].start;
		}
		double ratio = (double) coverageSum / (double) length;
		return ratio;
	}

	bool isValid() {
		return avg_covereage() > 30;
	}
};

class CigarExon {
public:
	int length;
	std::vector<is::CigarOp> cigar;
	bool isGap;
	CigarExon(int _length, std::vector<is::CigarOp> _cigar, bool _isGap) {
		length = _length;
		cigar = _cigar;
		isGap = _isGap;
	}
	void expandGap(int value) {
		cigar.front().count += value;
		length = cigar.front().count;
	}
};

class RealignmentStructure {
public:
	int order_number;
	const SingleSequence* sequence;
	int start;
	int stop;
	double score;
	bool isAligned;
	int64_t ref_number;
	std::vector<uint8_t> raw_alignment;
	SeqOrientation orientation;
	std::vector<CigarExon> previousCigarExons;

	RealignmentStructure(int _order_number, const SingleSequence* _sequence, PathGraphEntry* _entry, int64_t _ref_number, double _score, std::vector<CigarExon> _previousCigarExons) {
		order_number = _order_number;
		sequence = _sequence;
		ref_number = _ref_number;
		score = _score;
		previousCigarExons = _previousCigarExons;
		if (_entry != NULL) {
			std::vector<AlignmentResults> alignments =_entry->get_alignments();
			if (alignments.size() > 0) {
				AlignmentResults ar = alignments.back();
				raw_alignment = ar.raw_alignment;
				if (ar.orientation == kReverse) {
					reverse(raw_alignment.begin(), raw_alignment.end());
				}
				orientation = ar.orientation;
				start = ar.ref_start;
				stop = ar.ref_end;
			} else {
				start = 0;
				stop = 0;
			}
			isAligned = true;
		} else {
			isAligned = false;
		}
	}
};

class GraphMap {
 public:
  GraphMap();
  ~GraphMap();

  static bool comparePtrToNode(RealignmentStructure* a, RealignmentStructure* b);

  // Main function for running the mapping process. It generates/loads the index, and handles batch loading of sequences from the reads file.
  void Run(ProgramParameters &parameters);

  // Generates or loads the index of the reference genome.
  int BuildIndexes(ProgramParameters &parameters);

  double RealignRead(const SingleSequence *read, std::shared_ptr<is::MinimizerIndex> index, MappingData *mapping_data, const ProgramParameters *parameters, std::string cutted_reference, ExonsCluster exonsClusters, SeqOrientation orientation, int64_t ref_number, std::vector<CigarExon> *cigarExons);
  // Loads reads from a file in batches of given size (in MiB), or all at once.
  void ProcessReadsFromSingleFile(ProgramParameters &parameters, FILE *fp_out);

  // Process the loaded batch of reads. Uses OpenMP to do it in parallel. Calls ProcessOneRead for each read in the SequenceFile.
  int ProcessSequenceFileInParallel(ProgramParameters *parameters, const SequenceFile *reads, clock_t *last_time, FILE *fp_out, int64_t *ret_num_mapped, int64_t *ret_num_unmapped);

  // Processes a single read from the batch of loaded reads.
  int ProcessRead(int order_number, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params, std::vector<RealignmentStructure *> *realignment_structures);

  // Collects alignments from the given mapping_data and converts them into an appropriate output format (string).
  int CollectAlignments(const SingleSequence *read, const ProgramParameters *parameters, MappingData *mapping_data, std::string &ret_aln_lines);

  // Allows the usage of GraphMap as an API.
  int Align(const SequenceFile *ref, const SequenceFile *reads, const ProgramParameters &parameters);

  // Allows the usage of GraphMap as an API. Converts the std::string objects to SingleSequence and initializes SequenceFiles,
  // after which the above function is called. Headers for sequences are automatically generated: ref_%d and query_%d for ref_seqs and read_seqs, respectivelly.
  int Align(std::vector<std::string> ref_seqs, std::vector<std::string> read_seqs, const ProgramParameters &parameters);

  void PostprocessRNAData(std::vector<RealignmentStructure *> realignment_structures, std::vector<std::string> *sam_lines, int64_t num_threads, ProgramParameters *parameters, EValueParams *evalue_params);

 private:
  std::vector<std::shared_ptr<is::MinimizerIndex>> indexes_;
  std::shared_ptr<is::Transcriptome> transcriptome_;

  // Opens the output SAM file for writing if the path is specified. If the path is empty, then output is set to STDOUT.
  FILE* OpenOutSAMFile_(std::string out_sam_path="");
  // Formats the SAM header from a given index.

  // Generates a default SAM line for unmapped reads.
  std::string GenerateUnmappedSamLine_(const SingleSequence *read, const std::string& unmapped_reason, int64_t verbose_sam_output) const;

  // Calculates the LCSk of the anchors using the Fenwick tree.
  void CalcLCSFromLocalScoresCacheFriendly_(const Vertices *vertices, bool use_l1_filtering, int64_t l, int64_t allowed_dist,
                                            int* ret_lcskpp_length, std::vector<int> *ret_lcskpp_indices, int64_t allowed_anchor_overlap);

  std::shared_ptr<is::MinimizerIndex> SetupIndex_(std::shared_ptr<SequenceFile> ref, const std::string &index_path, const std::string &shape,
                  const ProgramParameters &parameters, int64_t num_threads) const;

  // Count gapped spaced seed hits to regions on the reference.
  int RegionSelectionNoCopy_(int64_t bin_size, MappingData *mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence *read, const ProgramParameters *parameters, std::vector<uint64_t> &keys);
  int RegionSelectionWithSort_(int64_t bin_size, MappingData *mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence *read, const ProgramParameters *parameters, std::vector<Region>& regions);
  void AppendSeedHits_(const uint128_t& seed, std::shared_ptr<is::MinimizerIndex> index, bool threshold_hits, double count_cutoff, std::vector<uint128_t> &all_hits);

  int GraphMap_(ScoreRegistry *local_score, std::shared_ptr<is::MinimizerIndex> index_read, MappingData *mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence *read, const ProgramParameters *parameters, std::vector<const uint128_t*>& hits,std::vector<int64_t> &num_hits, std::vector<uint64_t> &keys);
  int ProcessKmerCacheFriendly_(int8_t *kmer, int64_t kmer_start_position, ScoreRegistry *local_score, MappingData* mapping_data, std::shared_ptr<is::MinimizerIndex> index_read, const SingleSequence* read, const ProgramParameters* parameters, std::vector<const uint128_t*>& hits_ptr,std::vector<int64_t> &num_hits_ptr, std::vector<uint64_t> &keys, uint64_t *buffer, bool is_start);

  // Perform the LCSk calculation and simple filtering of the anchores that survived the LCSk.
  int SemiglobalPostProcessRegionWithLCS_(ScoreRegistry *local_score, MappingData *mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence *read, const ProgramParameters *parameters);
  int AnchoredPostProcessRegionWithLCS_(ScoreRegistry *local_score, MappingData *mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence *read, const ProgramParameters *parameters, int64_t allowed_anchor_overlap);

  // Performs a knapsack algorithm implementation on the set of clusters, to determine the most likely RNA-seq alignment. Clusters will be marked as valid or invalid. Valid should stay in play, but invalid filtered out.
  int RNAFilterClusters_(MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters);
  void RNAFilterClustersWithHugeGaps_(MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters);
  int CleanupIntermediateMappings_(MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters);

  //
  int EvaluateMappings_(MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int GenerateAlignments_(MappingData *mapping_data, std::shared_ptr<is::MinimizerIndex> index, std::shared_ptr<is::Transcriptome> transcriptome, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params);
  int RNAGenerateAlignments_(int order_number, MappingData *mapping_data, std::shared_ptr<is::MinimizerIndex> index, std::shared_ptr<is::Transcriptome> transcriptome, const SingleSequence *read, const ProgramParameters *parameters, const EValueParams *evalue_params, std::vector<RealignmentStructure *> *realignment_structures);

  // Helper functions.
  int CalculateL1ParametersWithMaximumDeviation_(ScoreRegistry *local_score, std::vector<int> &lcskpp_indices, float maximum_allowed_deviation, int64_t *ret_k, int64_t *ret_l, float *ret_sigma_L2, float *ret_confidence_L1);
  int64_t CountBinsWithinThreshold_(const MappingData *mapping_data, float threshold);
  Region CalcRegionFromBin_(int64_t sorted_bins_index, const MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int CheckRegionSearchFinished_(int64_t current_region, float min_allowed_bin_value, float threshold_step, float *bin_value_threshold, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int CheckRegionSearchFinished2_(int64_t current_region, float min_allowed_bin_value, float threshold_step,
                                            float *bin_value_threshold, MappingData *mapping_data, const std::vector<Region> &regions,
                                            const SingleSequence *read, const ProgramParameters *parameters);
  int CollectFinalMappingsAndMapQ_(bool generate_final_mapping_ptrs, MappingData *mapping_data, const SingleSequence *read, const ProgramParameters *parameters);
  int CheckMinimumMappingConditions_(MappingResults *mapping_data, L1Results *l1_data, std::shared_ptr<is::MinimizerIndex> index, const SingleSequence *read, const ProgramParameters *parameters);

  // Debug.
  int VerboseLocalScoresToFile(std::string file_path, const SingleSequence *read, const ScoreRegistry *local_score, const std::vector<int> *indices, int64_t l_median, float maximum_allowed_deviation, bool check_median_filtering, std::vector<int32_t> *cluster_ids=NULL);

  void VerboseRegions_(const ProgramParameters* parameters, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const std::vector<Region>& regions);
  void OpenDebugClustersFile_(const ProgramParameters* parameters, const SingleSequence* read);
};

#endif /* TREEMAP_SE_H_ */
