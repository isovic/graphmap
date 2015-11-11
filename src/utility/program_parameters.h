/*
 * program_parameters.h
 *
 *  Created on: Jul 24, 2014
 *      Author: ivan
 */

#ifndef PROGRAM_PARAMETERS_H_
#define PROGRAM_PARAMETERS_H_

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <sstream>
//#include "utility/utility_basic_defines.h"
#include "utility/program_version.h"
#include "utility/program_licence.h"

#define DEFUALT_K_REGION      (int64_t) 13
#define DEFAULT_K_GRAPH       (int64_t) 6
#define DEFAULT_NUM_LINKS     (int64_t) 9 // 18
#define DEFAULT_MIN_NUM_ANCHOR_BASES      (int64_t) 12
#define DEFAULT_ERROR_RATE    (float) 0.45f
#define DEFAULT_START_READ    (int64_t) 0
#define DEFAULT_NUM_READS_TO_PROCESS      ((int64_t) -1)
#define DEFAULT_DEBUG_READ                ((int64_t) -1)
#define DEFAULT_DEBUG_READ_BY_QNAME       std::string("");
#define DEFAULT_NUM_THREADS               ((int64_t) -1)
#define DEFAULT_REFERENCE_PATH            std::string("")
#define DEFAULT_INDEX_REFERENCE_PATH      std::string("")
#define DEFAULT_READS_PATH                std::string("")
#define DEFAULT_OUT_SAM_PATH              std::string("")
#define DEFAULT_VERBOSE_SAM_OUTPUT        (int64_t) 0
#define DEFAULT_VERBOSE_LEVEL             (int64_t) 5
#define DEFAULT_COMMAND_LINE              std::string("")
//#define DEFAULT_MAX_NUM_REGIONS_CUTOFF    100
//#define DEFAULT_MAX_NUM_REGIONS           1000
//#define DEFAULT_MAX_NUM_HITS              2000
#define DEFAULT_MAX_NUM_REGIONS_CUTOFF    (int64_t) 0
#define DEFAULT_MAX_NUM_REGIONS           (int64_t) 0
#define DEFAULT_MAX_NUM_HITS              (int64_t) 0 // 1000
#define DEFAULT_SKIP_MULTIPLE_KMERS_PER_BIN     true
#define DEFAULT_FILTER_BY_EDIT_DISTANCE         true
#define DEFAULT_OUTPUT_IN_ORIGINAL_ORDER        false
#define DEFAULT_KMER_STEP                 (int64_t) 1
#define DEFAULT_READS_FOLDER                    std::string("")
#define DEFAULT_OUTPUT_FOLDER                   std::string("")
#define DEFAULT_PROCESS_READS_FROM_FOLDER       false
#define DEFAULT_BATCH_SIZE_IN_MB                (int64_t) 200

#define DEFAULT_ALIGNMENT_ALGORITHM     std::string("myers")
#define DEFAULT_ALIGNMENT_APPROACH        std::string("sg")
#define DEFAULT_CALCULATE_ONLY_INDEX      false

#define DEFAULT_MATCH_SCORE               (int64_t) 5 // 5 // 4 // 1
#define DEFAULT_MEX_SCORE        (int64_t) 1
#define DEFAULT_MISMATCH_PENALTY          (int64_t) 4 // 4 // 4 // 3
#define DEFAULT_GAP_OPEN_PENALTY          (int64_t) 8 // 8 // 5 //1 // 6 // 5
#define DEFAULT_GAP_EXTEND_PENALTY        (int64_t) 6 // 6 // 2 // 1 // 2
// Probane kombinacije za E-value i komentari:
// 5, 4, 8, 6 - Nije skroz lose, ali se dosta poklapa s BWA-MEM-om na Salmonela Typhi datasetu, ali je BWA-MEM za nijansu bolji na koljenu.
// 1, 1, 2, 2 - Lose, koljeno je bilo na puno manjem recallu cak i od LAST-a.
// 4, 5, 6, 5 - Lose, koljeno je tocno na pola izmedju LAST-a i BWA-MEM-a.
// 5, 4, 10, 6 - kao 4-5-6-5, ali malo losija mozda
// 5, 4, 8, 6 ali sam popravio da zbrojim i gap extend za 1bp indel - najlosije do sada, manje od LAST-a.

#define DEFAULT_IS_REFERENCE_CIRCULAR     false

#define DEFAULT_COMPOSITE_PARAMETERS      std::string("")

#define DEFAULT_MARGIN_FOR_AMBIGUITY     (float)  0.05f
#define DEFAULT_OUTPUT_MULTIPLE_ALIGNMENTS  (int64_t) 0

#define DEFAULT_PARSIMONIOUS_MODE         false

#define DEFAULT_MAX_EVALUE_THRESHOLD      (double) -1.0
#define DEFAULT_MIN_MAPQ_VALUE            0

#define DEFAULT_OUTPUT_FORMAT             std::string("sam")



struct ProgramParameters {
  int64_t k_region = DEFUALT_K_REGION;                    // 'j', Kmer size for region search (binning).
  int64_t k_graph = DEFAULT_K_GRAPH;                      // 'k', Kmer size for graph building.
  int64_t num_links = DEFAULT_NUM_LINKS;                   // 'l', Number of backward edges to check.
  float error_rate = DEFAULT_ERROR_RATE;                 // 'e', Approximate error rate of the input read sequences.
  int64_t start_read = DEFAULT_START_READ;                   // 's', Start processing reads from the one specified.
  int64_t num_reads_to_process = DEFAULT_NUM_READS_TO_PROCESS;        // 'n', Number of reads to process. If equal to -1, all reads will be processed.
  int64_t debug_read = DEFAULT_DEBUG_READ;                  // 'y', Verbose output for read marked with this variable.
  std::string debug_read_by_qname = DEFAULT_DEBUG_READ_BY_QNAME;
  int64_t num_threads = DEFAULT_NUM_THREADS;                 // 't', Number of threads to use. If equal to -1, number of threads will be equal to number of processors.
  std::string reference_path = DEFAULT_REFERENCE_PATH;          // 'r', The path to the reference file.
  std::string index_reference_path = DEFAULT_INDEX_REFERENCE_PATH;    // 'i', The path to the reference file's index. If it does not exist, index will be created in this path.
  std::string reads_path = DEFAULT_READS_PATH;              // 'd', The path to the reads file, in FASTA or FASTQ format.
  std::string out_sam_path = DEFAULT_OUT_SAM_PATH;            // 'o', The output path. If left blank, all sam output will be placed to stdout.
  int64_t verbose_sam_output = DEFAULT_VERBOSE_SAM_OUTPUT;           // 'b', Helpful debug comments can be placed in SAM output lines (at the end), however, some tools (like SAMtools) don't quite like them. Comments can be turned of by setting this variable to 0. Different values increase/decrease verbosity level.
  int64_t verbose_level = DEFAULT_VERBOSE_LEVEL;                // 'v', Verbose level. If equal to 0 nothing except strict output will be placed on stdout.
  std::string command_line = DEFAULT_COMMAND_LINE;            // The actuall commandline that was used to generate the parameters.
  int64_t max_num_regions_cutoff = DEFAULT_MAX_NUM_REGIONS_CUTOFF;     // 'q' Before the read is skipped, it will be attempted to reduce the number of selected regions if their number is higher than max_num_regions_cutoff.
  int64_t max_num_regions = DEFAULT_MAX_NUM_REGIONS;            // 'g' If still more regions than this are selected, the read is too ambiguous for processing, so it will be skipped.

  // Binning parameters
  int64_t max_num_hits = DEFAULT_MAX_NUM_HITS;               // 'm' Maximum number of hits per kmer during the binning process.
  bool skip_multiple_kmers_per_bin = DEFAULT_SKIP_MULTIPLE_KMERS_PER_BIN;  // 'p' One kmer of a read can have multiple hits withing the same bin. If true, this parameter prevents this.

  // Finalization.
  bool filter_by_edit_distance = DEFAULT_FILTER_BY_EDIT_DISTANCE;      // 'f' If a read passed all the tests, it's edit distance is still measured at the end. If true, the edit distance has to be lower than a certain threshold to call a read 'mapped'.

  bool output_in_original_order = DEFAULT_OUTPUT_IN_ORIGINAL_ORDER;    // 'u' If true, SAM alignments will be output after the processing has finished, in the order of input reads.
  int64_t kmer_step = DEFAULT_KMER_STEP;              // 'w' The number of bases to skip between beginnings of every adjecent kmer.

  std::string reads_folder = DEFAULT_READS_FOLDER;            // 'D', The path to a folder that contains reads, in FASTA or FASTQ format. Intended for batch processing.
  std::string output_folder = DEFAULT_OUTPUT_FOLDER;           // 'O', The path to the output folder for batch processing.

  bool process_reads_from_folder = DEFAULT_PROCESS_READS_FROM_FOLDER;

  int64_t batch_size_in_mb = DEFAULT_BATCH_SIZE_IN_MB;             // 'B', specifies the size of a batch for sequence loading. If <= 0, all sequences will be loaded at once, otherwise the specified number of megabytes will be loaded consequentially.

  std::string alignment_algorithm = DEFAULT_ALIGNMENT_ALGORITHM;  // 'w', specifies whether EDlib or SSW or hybrid should be used for realignment in the last step.
  std::string alignment_approach = DEFAULT_ALIGNMENT_APPROACH;

  bool calc_only_index = DEFAULT_CALCULATE_ONLY_INDEX;

  int64_t match_score = DEFAULT_MATCH_SCORE;
  int64_t mex_score = DEFAULT_MEX_SCORE;
  int64_t mismatch_penalty = DEFAULT_MISMATCH_PENALTY;
  int64_t gap_open_penalty = DEFAULT_GAP_OPEN_PENALTY;
  int64_t gap_extend_penalty = DEFAULT_GAP_EXTEND_PENALTY;

  int64_t evalue_match = 5;
  int64_t evalue_mismatch = 4;
  int64_t evalue_gap_open = 8;
  int64_t evalue_gap_extend = 6;
  
  bool is_reference_circular = DEFAULT_IS_REFERENCE_CIRCULAR;     // 'C'

  std::string composite_parameters = DEFAULT_COMPOSITE_PARAMETERS;  // 'x', specifies several parameters at the same time, such as 'nanopore' and 'illumina'.

  float margin_for_ambiguity = DEFAULT_MARGIN_FOR_AMBIGUITY;  // All mapping positions within the given fraction of the top score will be counted for ambiguity (mapping quality). Value of 0.0f counts only identical mappings.
  int64_t output_multiple_alignments = DEFAULT_OUTPUT_MULTIPLE_ALIGNMENTS;  // If 0, only one best alignment will be output. Otherwise, all alignments within margin_for_ambiguity will be output to a file.

  bool parsimonious_mode = DEFAULT_PARSIMONIOUS_MODE; // If true, only one index will be used, but the memory consumption will be reduced by half. If false, the fast and memory-hungry mode will be used.
  int64_t min_num_anchor_bases = DEFAULT_MIN_NUM_ANCHOR_BASES;

  double evalue_threshold = DEFAULT_MAX_EVALUE_THRESHOLD;
  int64_t mapq_threshold = DEFAULT_MIN_MAPQ_VALUE;

  std::string outfmt = DEFAULT_OUTPUT_FORMAT;
};

int ProcessArgs(int argc, char **argv, ProgramParameters *parameters);
void VerboseProgramParameters(ProgramParameters *parameters);
void VerboseShortHelpAndExit(int argc, char **argv);
void VerboseUsageAndExit(FILE *fp=stdout);
void VerboseEula(FILE *fp=stdout);



#endif /* PROGRAM_PARAMETERS_H_ */
