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

#define SOFTWARE_NAME "GraphMap"
#define GRAPHMAP_CURRENT_VERSION "v0.3.2"
#define GRAPHMAP_CURRENT_VERSION_RELEASE_DATE (std::string(__DATE__) + std::string(" at ") + std::string(__TIME__)) // __TIMESTAMP__ // "12 October 2014"
#define COPYRIGHT "Copyright Ivan Sovic, Mile Sikic and Niranjan Nagarajan, 2015, 2016.\n" \
                  "\n" \
                  "Affiliations: Ivan Sovic (1, 3), Mile Sikic (2), Niranjan Nagarajan (3)\n" \
                  "  (1) Ruder Boskovic Institute, Zagreb, Croatia\n" \
                  "  (2) University of Zagreb, Faculty of Electrical Engineering and Computing\n" \
                  "  (3) Genome Institute of Singapore, A*STAR, Singapore\n"

#define LICENCE_INFORMATION   \
  "GraphMap (c) by Ivan Sovic, Mile Sikic and Niranjan Nagarajan\n" \
  "GraphMap is licensed under The MIT License.\n" \

#define AFFILIATIONS \
  "Affiliations: Ivan Sovic (1, 3), Mile Sikic (2), Niranjan Nagarajan (3)\n" \
  "  (1) Ruder Boskovic Institute, Zagreb, Croatia\n" \
  "  (2) University of Zagreb, Faculty of Electrical Engineering and Computing\n" \
  "  (3) Genome Institute of Singapore, A*STAR, Singapore\n"

struct ProgramParameters {
  std::string subprogram = "";

  int64_t k_region = 13;                    // 'j', Kmer size for region search (binning).
  int64_t k_graph = 6;                      // 'k', Kmer size for graph building.
  int64_t num_links = 9;                   // 'l', Number of backward edges to check.
  float error_rate = 0.45;                 // 'e', Approximate error rate of the input read sequences.
  int64_t start_read = 0;                   // 's', Start processing reads from the one specified.
  int64_t num_reads_to_process = -1;        // 'n', Number of reads to process. If equal to -1, all reads will be processed.
  int64_t debug_read = -1;                  // 'y', Verbose output for read marked with this variable.
  std::string debug_read_by_qname = "";
  int64_t num_threads = -1;                 // 't', Number of threads to use. If equal to -1, number of threads will be equal to number of processors.
  std::string reference_path = "";          // 'r', The path to the reference file.
  std::string index_file = "";    // 'i', The path to the reference file's index. If it does not exist, index will be created in this path.
  std::string reads_path = "";              // 'd', The path to the reads file, in FASTA or FASTQ format.
  std::string out_sam_path = "";            // 'o', The output path. If left blank, all sam output will be placed to stdout.
  int64_t verbose_sam_output = 0;           // 'b', Helpful debug comments can be placed in SAM output lines (at the end), however, some tools (like SAMtools) don't quite like them. Comments can be turned of by setting this variable to 0. Different values increase/decrease verbosity level.
  int64_t verbose_level = 5;                // 'v', Verbose level. If equal to 0 nothing except strict output will be placed on stdout.
  std::string command_line = "";            // The actuall commandline that was used to generate the parameters.
  int64_t max_num_regions_cutoff = 0;     // 'q' Before the read is skipped, it will be attempted to reduce the number of selected regions if their number is higher than max_num_regions_cutoff.
  int64_t max_num_regions = 0;            // 'g' If still more regions than this are selected, the read is too ambiguous for processing, so it will be skipped.

  // Binning parameters
  int64_t max_num_hits = 0;               // 'm' Maximum number of hits per kmer during the binning process.
  bool skip_multiple_kmers_per_bin = true;  // 'p' One kmer of a read can have multiple hits withing the same bin. If true, this parameter prevents this.

  bool output_in_original_order = false;    // 'u' If true, SAM alignments will be output after the processing has finished, in the order of input reads.
  int64_t kmer_step = 1;              // 'w' The number of bases to skip between beginnings of every adjecent kmer.

  std::string reads_folder = "";            // 'D', The path to a folder that contains reads, in FASTA or FASTQ format. Intended for batch processing.
  std::string output_folder = "";           // 'O', The path to the output folder for batch processing.
  bool process_reads_from_folder = false;
  int64_t batch_size_in_mb = -1;             // 'B', specifies the size of a batch for sequence loading. If <= 0, all sequences will be loaded at once, otherwise the specified number of megabytes will be loaded consequentially.
  std::string alignment_algorithm = "sg";  // 'a', specifies whether EDlib or SSW or hybrid should be used for realignment in the last step.
  std::string alignment_approach = "normal";      // 'w'
  bool calc_only_index = false;
  int64_t match_score = 5;
  int64_t mex_score = 1;
  int64_t mismatch_penalty = 4;
  int64_t gap_open_penalty = 8;
  int64_t gap_extend_penalty = 6;
  int64_t evalue_match = 5;
  int64_t evalue_mismatch = 4;
  int64_t evalue_gap_open = 8;
  int64_t evalue_gap_extend = 6;
  bool is_reference_circular = false;     // 'C'
  std::string composite_parameters = "";  // 'x', specifies several parameters at the same time, such as 'nanopore' and 'illumina'.
  float margin_for_ambiguity = 0.05;  // All mapping positions within the given fraction of the top score will be counted for ambiguity (mapping quality). Value of 0.0f counts only identical mappings.
  bool output_multiple_alignments = false;  // If 0, only one best alignment will be output. Otherwise, all alignments within margin_for_ambiguity will be output to a file.
  bool sensitive_mode = false; // If false, only one index will be used, but the memory consumption will be reduced by half. If false, sensitive and memory-hungry mode will be used.
  int64_t min_num_anchor_bases = 12;
  double evalue_threshold = -1;
  int64_t mapq_threshold = 0;
  std::string infmt = "auto";
  std::string outfmt = "sam";

//  bool extend_aln_to_end = true;

  bool use_extended_cigar = false;

  int64_t min_read_len = 80;      // If a read is shorter than this, it will be marked as unmapped.

  double min_bin_percent = 0.75f;
  double bin_threshold_step = 0.10f;

  bool use_spliced = false;
  bool use_split = false;
  bool disable_end_to_end = true;
  bool overlapper = false;
  bool no_self_hits = false;
  bool rebuild_index = false;

  double max_error_rate = 1.0f;
  double max_indel_error_rate = 1.0f;

  std::string gtf_path;
};

int ProcessArgsGraphMap(int argc, char **argv, ProgramParameters *parameters);
int ProcessArgsOwler(int argc, char **argv, ProgramParameters *parameters);
void VerboseProgramParameters(ProgramParameters *parameters);
void VerboseShortHelpAndExit(int argc, char **argv);

#endif /* PROGRAM_PARAMETERS_H_ */
