/*
 * program_parameters.cc
 *
 *  Created on: Jul 24, 2014
 *      Author: ivan
 */

#include "program_parameters.h"
#include "argparser.h"
#include "utility/utility_general.h"

int ProcessArgsGraphMap(int argc, char **argv, ProgramParameters *parameters)
{
  bool help = false;

  ArgumentParser argparser;

  argparser.AddCompositeArgument("illumina", "-a gotoh -w normal -M 5 -X 4 -G 8 -E 6");
//  argparser.AddCompositeArgument("overlap", "-a anchor -w normal --overlapper --evalue 1e0 --ambiguity 0.50 --min-bin-perc 0.10 --bin-step 0.90 --max-regions -1 --mapq -1 --secondary");
  argparser.AddCompositeArgument("overlap", "-a anchor -w normal --overlapper --evalue 1e0 --ambiguity 0.50 --min-bin-perc 0.10 --bin-step 0.90 --max-regions -1 --mapq -1 --secondary");
//  argparser.AddCompositeArgument("sensitive", "-a gotoh -w sg -M 5 -X 4 -G 8 -E 6");

  argparser.AddArgument(&parameters->reference_path, VALUE_TYPE_STRING, "r", "ref", "", "Path to the reference sequence (fastq or fasta).", 0, "Input/Output options");
  argparser.AddArgument(&parameters->index_file, VALUE_TYPE_STRING, "i", "index", "", "Path to the index of the reference sequence. If not specified, index is generated in the same folder as the reference file, with .gmidx extension. For non-parsimonious mode, secondary index .gmidxsec is also generated.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->reads_path, VALUE_TYPE_STRING, "d", "reads", "", "Path to the reads file.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->out_sam_path, VALUE_TYPE_STRING, "o", "out", "", "Path to the output file that will be generated.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->infmt, VALUE_TYPE_STRING, "K", "in-fmt", "auto", "Format in which to input reads. Options are:\n auto  - Determines the format automatically from file extension.\n fastq - Loads FASTQ or FASTA files.\n fasta - Loads FASTQ or FASTA files.\n gfa   - Graphical Fragment Assembly format.\n sam   - Sequence Alignment/Mapping format.", 0, "Input/Output options");
//  argparser.AddArgument(&parameters->outfmt, VALUE_TYPE_STRING, "L", "out-fmt", "sam", "Format in which to output results. Options are:\n sam  - Standard SAM output (in normal and '-w overlap' modes).\n m5   - BLASR M5 format.\n mhap - MHAP overlap format (use with '-w owler').\n paf  - PAF (Minimap) overlap format (use with '-w owler').", 0, "Input/Output options");
  argparser.AddArgument(&parameters->outfmt, VALUE_TYPE_STRING, "L", "out-fmt", "sam", "Format in which to output results. Options are:\n sam  - Standard SAM output (in normal and '-w overlap' modes).\n m5   - BLASR M5 format.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->calc_only_index, VALUE_TYPE_BOOL, "I", "index-only", "0", "Build only the index from the given reference and exit. If not specified, index will automatically be built if it does not exist, or loaded from file otherwise.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->rebuild_index, VALUE_TYPE_BOOL, "", "rebuild-index", "0", "Rebuild index even if it already exists in given path.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->output_in_original_order, VALUE_TYPE_BOOL, "u", "ordered", "0", "SAM alignments will be output after the processing has finished, in the order of input reads.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->batch_size_in_mb, VALUE_TYPE_INT64, "B", "batch-mb", "1024", "Reads will be loaded in batches of the size specified in megabytes. Value <= 0 loads the entire file.", 0, "Input/Output options");
  //    argparser.AddArgument(&parameters->reads_folder, VALUE_TYPE_STRING, "D", "readsfolder", "", "Path to a folder containing read files (in fastq or fasta format) to process. Cannot be used in combination with '-d' or '-o'.", 0, "Input/Output options");
  //    argparser.AddArgument(&parameters->output_folder, VALUE_TYPE_STRING, "O", "outfolder", "", "Path to a folder for placing SAM alignments. Use in combination with '-D'.", 0, "Input/Output options");

  argparser.AddArgument(&parameters->composite_parameters, VALUE_TYPE_COMPOSITE, "x", "preset", "", "Pre-set parameters to increase sensitivity for different sequencing technologies. Valid options are:\n illumina - Equivalent to: '-a gotoh -w normal -M 5 -X 4 -G 8 -E 6'\n overlap  - Equivalent to: '-a anchor -w normal --overlapper --evalue 1e0 --ambiguity 0.50 --secondary'", 0, "General-purpose pre-set options");

  argparser.AddArgument(&parameters->alignment_algorithm, VALUE_TYPE_STRING, "a", "alg", "anchor", "Specifies which algorithm should be used for alignment. Options are:\n sg       - Myers' bit-vector approach. Semiglobal. Edit dist. alignment.\n sggotoh       - Gotoh alignment with affine gaps. Semiglobal.\n anchor      - anchored alignment with end-to-end extension.\n               Uses Myers' global alignment to align between anchors.\n anchorgotoh - anchored alignment with Gotoh.\n               Uses Gotoh global alignment to align between anchors.", 0, "Alignment options");
//  argparser.AddArgument(&parameters->alignment_approach, VALUE_TYPE_STRING, "w", "appr", "sg", "Additional alignment approaches. Changes the way alignment algorithm is applied. Options are:\n sg         - Normal (default) alignment mode (non-overlapping).\n overlapper - (Experimental) Runs the entire GraphMap pipeline with small\n              modifications for better overlapping. Output in SAM format.\n              This is also a composite parameter - it changes values of other params to:\n              '-a anchor -Z -F 0.50 -z 1e0'.\n owler      - (Experimental) Runs reduced pipeline, does not produce alignments, fast.\n              Output in MHAP format.", 0, "Alignment options");
#ifndef RELEASE_VERSION
  argparser.AddArgument(&parameters->alignment_approach, VALUE_TYPE_STRING, "w", "approach", "normal", "Additional alignment approaches. Changes the way alignment algorithm is applied. Options are:\n normal         - Normal alignment of reads to the reference.\n (Currently no other options are provided. This is a placeholder for future features, such as cDNA mapping)", 0, "Alignment options");
#endif
//  argparser.AddArgument(&parameters->alignment_approach, VALUE_TYPE_STRING, "w", "appr", "normal", "Additional alignment approaches. Changes the way alignment algorithm is applied. Options are:\n normal         - Normal alignment of reads to the reference.\n overlapper - (Experimental) Runs the entire GraphMap pipeline with small\n              modifications for better overlapping. Output in SAM format.\n              This is also a composite parameter - it changes values of other params to:\n              '-a anchor -Z -F 0.50 -z 1e0'.", 0, "Alignment options");
  argparser.AddArgument(&parameters->overlapper, VALUE_TYPE_BOOL, "", "overlapper", "0", "Perform overlapping instead of mapping. Skips self-hits if reads and reference files contain same sequences, and outputs lenient secondary alignments.", 0, "Alignment options");

  argparser.AddArgument(&parameters->match_score, VALUE_TYPE_INT64, "M", "match", "5", "Match score for the DP alignment. Ignored for Myers alignment.", 0, "Alignment options");
#ifndef RELEASE_VERSION
  argparser.AddArgument(&parameters->mex_score, VALUE_TYPE_INT64, "T", "mex", "1", "Mex score.", 0, "Alignment options");
#endif
  argparser.AddArgument(&parameters->mismatch_penalty, VALUE_TYPE_INT64, "X", "mismatch", "4", "Mismatch penalty for the DP alignment. Ignored for Myers alignment.", 0, "Alignment options");
  argparser.AddArgument(&parameters->gap_open_penalty, VALUE_TYPE_INT64, "G", "gapopen", "8", "Gap open penalty for the DP alignment. Ignored for Myers alignment.", 0, "Alignment options");
  argparser.AddArgument(&parameters->gap_extend_penalty, VALUE_TYPE_INT64, "E", "gapext", "6", "Gap extend penalty for the DP alignment. Ignored for Myers alignment.", 0, "Alignment options");
  argparser.AddArgument(&parameters->evalue_threshold, VALUE_TYPE_DOUBLE, "z", "evalue", "1e0", "Threshold for E-value. If E-value > FLT, read will be called unmapped. If FLT < 0.0, thredhold not applied.", 0, "Alignment options");
  argparser.AddArgument(&parameters->mapq_threshold, VALUE_TYPE_INT64, "c", "mapq", "1", "Threshold for mapping quality. If mapq < INT, read will be called unmapped.", 0, "Alignment options");
  argparser.AddArgument(&parameters->use_extended_cigar, VALUE_TYPE_BOOL, "", "extcigar", "0", "Use the extended CIGAR format for output alignments.", 0, "Alignment options");
//#ifndef RELEASE_VERSION
//  argparser.AddArgument(&parameters->use_spliced, VALUE_TYPE_BOOL, "", "spliced", "0", "Align clusters of anchors independently and report them as separate alignments. Does not align in-between clusters. Works only for anchored alignment modes.", 0, "Alignment options");
//  argparser.AddArgument(&parameters->use_split, VALUE_TYPE_BOOL, "", "split", "0", "Align clusters of anchors independently and report them as separate alignments. Does not align in-between clusters. Works only for anchored alignment modes.", 0, "Alignment options");
//#endif
  argparser.AddArgument(&parameters->disable_end_to_end, VALUE_TYPE_BOOL, "", "no-end2end", "0", "Disables extending of the alignments to the ends of the read. Works only for anchored modes.", 0, "Alignment options");
//  argparser.AddArgument(&parameters->extend_aln_to_end, VALUE_TYPE_BOOL, "", "extend-to-end", "1", "Switches extension of anchored alignment to read/reference ends. If false, alignments will be bound by the first and last anchor. Otherwise, alignment will proceed until an end of one of the sequences is hit.", 0, "Alignment options");

  argparser.AddArgument(&parameters->k_graph, VALUE_TYPE_INT64, "k", "", "6", "Graph construction kmer size.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->num_links, VALUE_TYPE_INT64, "l", "", "9", "Number of edges per vertex.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->min_num_anchor_bases, VALUE_TYPE_INT64, "A", "minbases", "12", "Minimum number of match bases in an anchor.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->error_rate, VALUE_TYPE_FLOAT, "e", "error-rate", "0.45", "Approximate error rate of the input read sequences.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->max_num_hits, VALUE_TYPE_INT64, "", "max-hits", "-1", "Maximum allowed number of hits per seed. If 0, all seeds will be used. If < 0, threshold will be calculated automatically.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->max_num_regions, VALUE_TYPE_INT64, "g", "max-regions", "0", "If the final number of regions exceeds this amount, the read will be called unmapped. If 0, value will be dynamically determined. If < 0, no limit is set.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->max_num_regions_cutoff, VALUE_TYPE_INT64, "q", "reg-reduce", "0", "Attempt to heuristically reduce the number of regions if it exceeds this amount. Value <= 0 disables reduction but only if param -g is not 0. If -g is 0, the value of this parameter is set to 1/5 of maximum number of regions.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->is_reference_circular, VALUE_TYPE_BOOL, "C", "circular", "0", "Reference sequence is a circular genome.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->margin_for_ambiguity, VALUE_TYPE_FLOAT, "F", "ambiguity", "0.02", "All mapping positions within the given fraction of the top score will be counted for ambiguity (mapping quality). Value of 0.0 counts only identical mappings.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->output_multiple_alignments, VALUE_TYPE_BOOL, "Z", "secondary", "0", "If specified, all (secondary) alignments within (-F FLT) will be output to a file. Otherwise, only one alignment will be output.", 0, "Algorithmic options");
//  argparser.AddArgument(&parameters->output_multiple_alignments, VALUE_TYPE_BOOL, "Z", "no-secondary", "1", "If specified, all (secondary) alignments within (-F FLT) will be output to a file. Otherwise, only one alignment will be output.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->sensitive_mode, VALUE_TYPE_BOOL, "P", "sensitive", "0", "If false, only one gapped spaced index will be used in region selection. If true, two such indexes (with different shapes) will be used (2x memory-hungry but more powerful for very high error rates).", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->min_bin_percent, VALUE_TYPE_DOUBLE, "", "min-bin-perc", "0.75", "Consider only bins with counts above FLT * max_bin, where max_bin is the count of the top scoring bin.", 0, "Algorithmic options");
//  argparser.AddArgument(&parameters->bin_threshold_step, VALUE_TYPE_DOUBLE, "", "bin-step", "0.10", "After a chunk of bins with values above FLT * max_bin is processed, check if there is one extremely dominant region, and stop the search.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->bin_threshold_step, VALUE_TYPE_DOUBLE, "", "bin-step", "0.25", "After a chunk of bins with values above FLT * max_bin is processed, check if there is one extremely dominant region, and stop the search.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->min_read_len, VALUE_TYPE_INT64, "", "min-read-len", "80", "If a read is shorter than this, it will be marked as unmapped. This value can be lowered if the reads are known to be accurate.", 0, "Algorithmic options");

  argparser.AddArgument(&parameters->num_threads, VALUE_TYPE_INT64, "t", "threads", "-1", "Number of threads to use. If '-1', number of threads will be equal to min(24, num_cores/2).", 0, "Other options");
  argparser.AddArgument(&parameters->verbose_level, VALUE_TYPE_INT64, "v", "verbose", "5", "Verbose level. If equal to 0 nothing except strict output will be placed on stdout.", 0, "Other options");
  argparser.AddArgument(&parameters->start_read, VALUE_TYPE_INT64, "s", "start", "0", "Ordinal number of the read from which to start processing data.", 0, "Other options");
  argparser.AddArgument(&parameters->num_reads_to_process, VALUE_TYPE_INT64, "n", "numreads", "-1", "Number of reads to process per batch. Value of '-1' processes all reads.", 0, "Other options");
  argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0, "Other options");

  argparser.AddArgument(&parameters->debug_read, VALUE_TYPE_INT64, "y", "debug-read", "-1", "ID of the read to give the detailed verbose output.", 0, "Debug options");
  argparser.AddArgument(&parameters->debug_read_by_qname, VALUE_TYPE_STRING, "Y", "debug-qname", "", "QNAME of the read to give the detailed verbose output. Has precedence over -y. Use quotes to specify.", 0, "Debug options");
  argparser.AddArgument(&parameters->verbose_sam_output, VALUE_TYPE_INT64, "b", "verbose-sam", "0", "Helpful debug comments can be placed in SAM output lines (at the end). Comments can be turned off by setting this parameter to 0. Different values increase/decrease verbosity level.\n0 - verbose off\n1 - server mode, command line will be omitted to obfuscate paths.\n2 - umm this one was skipped by accident. The same as 0.\n>=3 - detailed verbose is added for each alignment, including timing measurements and other.\n4 - qnames and rnames will not be trimmed to the first space.\n5 - QVs will be omitted (if available).", 0, "Debug options");

  argparser.ProcessArguments(argc, argv);

  // Store the command line in parameters.
  std::stringstream ss_command_line;
  for (int i=0; i<argc; i++) {
    if (i > 0)
      ss_command_line << " ";
    ss_command_line << argv[i];
  }
  parameters->command_line = ss_command_line.str();



  /// Check if help was triggered.
  if (argparser.GetArgumentByLongName("help")->is_set == true) {
    std::stringstream ss;
    ss << SOFTWARE_NAME << " - A very accurate and sensitive long-read, high error-rate sequence mapper\n", SOFTWARE_NAME;
    ss << "Version: " <<  GRAPHMAP_CURRENT_VERSION << "\n";
    ss << "Build date: " <<  std::string(GRAPHMAP_CURRENT_VERSION_RELEASE_DATE).c_str() << "\n";
    ss << "\n";
    ss << LICENCE_INFORMATION << "\n";
    ss << AFFILIATIONS << "\n";
    ss << "\n";
    ss << "Usage:\n";
    ss << "\tgraphmap [options] -r <reference_file> -d <reads_file> -o <output_sam_path>\n";
    ss << "\n";

    fprintf (stderr, "%s\n", ss.str().c_str());
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    exit(1);
  }

  /// In case debug_read_by_qname was set, and it was specified using quote signs on the command line, remove the quotes.
  if (parameters->debug_read_by_qname.size() > 2 && ((parameters->debug_read_by_qname.front())) == '"' && ((parameters->debug_read_by_qname.back())) == '"') {
    parameters->debug_read_by_qname = parameters->debug_read_by_qname.substr(1, (parameters->debug_read_by_qname.size() - 2));
  }

//  /// For the 'overlapper' mode use anchored alignment by default, no other alignments should be possible. Except perhaps anchorgotoh...
//  if (parameters->alignment_approach == "overlapper") {
//    if (parameters->alignment_algorithm != "anchor" && parameters->alignment_algorithm != "anchorgotoh") {
//      parameters->alignment_algorithm = "anchor";
//    }
//    if (argparser.GetArgumentByLongName("evalue")->is_set == false) {
//      parameters->evalue_threshold = 1e0;
//    }
//    if (argparser.GetArgumentByLongName("ambiguity")->is_set == false) {
//      parameters->margin_for_ambiguity = 0.50f;
//    }
//    if (argparser.GetArgumentByLongName("secondary")->is_set == false) {
//      parameters->output_multiple_alignments = 1;
//    }
//  }

  // Sanity check for the reference path.
  if (argparser.GetArgumentByLongName("ref")->is_set == false) {
    fprintf (stderr, "Please specify the path to the reference file.\n");
    fprintf (stderr, "\n");
    VerboseShortHelpAndExit(argc, argv);
  }
  if (!fileExists(parameters->reference_path.c_str())) {
      fprintf (stderr, "Reference does not exist: '%s'\n\n", parameters->reference_path.c_str());
      VerboseShortHelpAndExit(argc, argv);
  }

  // Sanity check for the reads path.
  if (argparser.GetArgumentByLongName("reads")->is_set == false && parameters->calc_only_index == false) {
    fprintf (stderr, "Please specify the path to the reads file.\n");
    fprintf (stderr, "\n");
    VerboseShortHelpAndExit(argc, argv);
  }
  if (!fileExists(parameters->reads_path.c_str())) {
    fprintf (stderr, "Reads file does not exist: '%s'\n\n", parameters->reads_path.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

  // Check if the index path was specified, if not, generate it.
  if (argparser.GetArgumentByLongName("index")->is_set == false) {
    parameters->index_file = parameters->reference_path + std::string(".gmidx");
  }

  // Check if the input format is valid.
  if (parameters->infmt != "auto" && parameters->infmt != "fastq" &&
      parameters->infmt != "fasta" && parameters->infmt != "gfa" && parameters->infmt != "sam") {
    fprintf (stderr, "Unknown input format '%s'!\n\n", parameters->infmt.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

  // Check if the output format is valid.
  if (parameters->outfmt != "sam" && parameters->outfmt != "m5" &&
      parameters->outfmt != "mhap" && parameters->outfmt != "paf") {
    fprintf (stderr, "Unknown output format '%s'!\n\n", parameters->outfmt.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

  if (parameters->alignment_algorithm != "sg" && parameters->alignment_algorithm != "sggotoh" &&
      parameters->alignment_algorithm != "anchor" && parameters->alignment_algorithm != "anchorgotoh" && parameters->alignment_algorithm != "anchormex") {
    fprintf (stderr, "Unknown alignment algorithm '%s'!\n\n", parameters->alignment_algorithm.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

#ifndef RELEASE_VERSION
  if (parameters->debug_read >= 0 || parameters->debug_read_by_qname != "") {
    parameters->verbose_level = 9;
    parameters->num_threads = 1;
    parameters->num_reads_to_process = 1;
  }
#endif

  /// Write this out for every debug verbose level.
  if (parameters->verbose_level > 5) {
    fprintf (stderr, "%s\n", argparser.VerboseArguments().c_str());
  }

//  VerboseProgramParameters(parameters);

  return 0;
}

int ProcessArgsOwler(int argc, char **argv, ProgramParameters *parameters)
{
  bool help = false;

  ArgumentParser argparser;

  argparser.AddArgument(&parameters->reference_path, VALUE_TYPE_STRING, "r", "ref", "", "Path to the reference sequence (fastq or fasta).", 0, "Input/Output options");
  argparser.AddArgument(&parameters->index_file, VALUE_TYPE_STRING, "i", "index", "", "Path to the index of the reference sequence. If not specified, index is generated in the same folder as the reference file, with .gmidx extension. For non-parsimonious mode, secondary index .gmidxsec is also generated.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->reads_path, VALUE_TYPE_STRING, "d", "reads", "", "Path to the reads file.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->out_sam_path, VALUE_TYPE_STRING, "o", "out", "", "Path to the output file that will be generated.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->infmt, VALUE_TYPE_STRING, "K", "in-fmt", "auto", "Format in which to input reads. Options are:\n auto  - Determines the format automatically from file extension.\n fastq - Loads FASTQ or FASTA files.\n fasta - Loads FASTQ or FASTA files.\n gfa   - Graphical Fragment Assembly format.\n sam   - Sequence Alignment/Mapping format.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->outfmt, VALUE_TYPE_STRING, "L", "out-fmt", "mhap", "Format in which to output results. Options are:\n mhap - MHAP overlap format.\n paf  - PAF (Minimap) overlap format.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->calc_only_index, VALUE_TYPE_BOOL, "I", "index-only", "0", "Build only the index from the given reference and exit. If not specified, index will automatically be built if it does not exist, or loaded from file otherwise.", 0, "Input/Output options");
//  argparser.AddArgument(&parameters->rebuild_index, VALUE_TYPE_BOOL, "", "rebuild-index", "0", "Rebuild index even if it already exists in given path.", 0, "Input/Output options");
  argparser.AddArgument(&parameters->batch_size_in_mb, VALUE_TYPE_INT64, "B", "batch-mb", "1024", "Reads will be loaded in batches of the size specified in megabytes. Value <= 0 loads the entire file.", 0, "Input/Output options");

  argparser.AddArgument(&parameters->error_rate, VALUE_TYPE_FLOAT, "e", "error-rate", "0.45", "Approximate error rate of the input read sequences.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->max_num_hits, VALUE_TYPE_INT64, "", "max-hits", "-1", "Maximum allowed number of hits per seed. If 0, all seeds will be used. If < 0, threshold will be calculated automatically.", 0, "Algorithmic options");
  argparser.AddArgument(&parameters->min_read_len, VALUE_TYPE_INT64, "", "min-read-len", "80", "If a read is shorter than this, it will be skipped. This value can be lowered if the reads are known to be accurate.", 0, "Algorithmic options");

  argparser.AddArgument(&parameters->num_threads, VALUE_TYPE_INT64, "t", "threads", "-1", "Number of threads to use. If '-1', number of threads will be equal to min(24, num_cores/2).", 0, "Other options");
  argparser.AddArgument(&parameters->verbose_level, VALUE_TYPE_INT64, "v", "verbose", "5", "Verbose level. If equal to 0 nothing except strict output will be placed on stdout.", 0, "Other options");
  argparser.AddArgument(&parameters->start_read, VALUE_TYPE_INT64, "s", "start", "0", "Ordinal number of the read from which to start processing data.", 0, "Other options");
  argparser.AddArgument(&parameters->num_reads_to_process, VALUE_TYPE_INT64, "n", "numreads", "-1", "Number of reads to process per batch. Value of '-1' processes all reads.", 0, "Other options");
  argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0, "Other options");

  argparser.AddArgument(&parameters->debug_read, VALUE_TYPE_INT64, "y", "debug-read", "-1", "ID of the read to give the detailed verbose output.", 0, "Debug options");
  argparser.AddArgument(&parameters->debug_read_by_qname, VALUE_TYPE_STRING, "Y", "debug-qname", "", "QNAME of the read to give the detailed verbose output. Has precedence over -y. Use quotes to specify.", 0, "Debug options");

  argparser.ProcessArguments(argc, argv);

  // Store the command line in parameters.
  std::stringstream ss_command_line;
  for (int i=0; i<argc; i++) {
    if (i > 0)
      ss_command_line << " ";
    ss_command_line << argv[i];
  }
  parameters->command_line = ss_command_line.str();



  /// Check if help was triggered.
  if (argparser.GetArgumentByLongName("help")->is_set == true) {
    std::stringstream ss;
    ss << SOFTWARE_NAME << " - A very accurate and sensitive long-read, high error-rate sequence mapper\n", SOFTWARE_NAME;
    ss << "Version: " <<  GRAPHMAP_CURRENT_VERSION << "\n";
    ss << "Build date: " <<  std::string(GRAPHMAP_CURRENT_VERSION_RELEASE_DATE).c_str() << "\n";
    ss << "\n";
    ss << LICENCE_INFORMATION << "\n";
    ss << AFFILIATIONS << "\n";
    ss << "\n";
    ss << "Usage:\n";
    ss << "\tgraphmap [options] -r <reference_file> -d <reads_file> -o <output_sam_path>\n";
    ss << "\n";

    fprintf (stderr, "%s\n", ss.str().c_str());
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    exit(1);
  }

  /// In case debug_read_by_qname was set, and it was specified using quote signs on the command line, remove the quotes.
  if (parameters->debug_read_by_qname.size() > 2 && ((parameters->debug_read_by_qname.front())) == '"' && ((parameters->debug_read_by_qname.back())) == '"') {
    parameters->debug_read_by_qname = parameters->debug_read_by_qname.substr(1, (parameters->debug_read_by_qname.size() - 2));
  }

//  /// For the 'overlapper' mode use anchored alignment by default, no other alignments should be possible. Except perhaps anchorgotoh...
//  if (parameters->alignment_approach == "overlapper") {
//    if (parameters->alignment_algorithm != "anchor" && parameters->alignment_algorithm != "anchorgotoh") {
//      parameters->alignment_algorithm = "anchor";
//    }
//    if (argparser.GetArgumentByLongName("evalue")->is_set == false) {
//      parameters->evalue_threshold = 1e0;
//    }
//    if (argparser.GetArgumentByLongName("ambiguity")->is_set == false) {
//      parameters->margin_for_ambiguity = 0.50f;
//    }
//    if (argparser.GetArgumentByLongName("secondary")->is_set == false) {
//      parameters->output_multiple_alignments = 1;
//    }
//  }

  // Sanity check for the reference path.
  if (argparser.GetArgumentByLongName("ref")->is_set == false) {
    fprintf (stderr, "Please specify the path to the reference file.\n");
    fprintf (stderr, "\n");
    VerboseShortHelpAndExit(argc, argv);
  }
  if (!fileExists(parameters->reference_path.c_str())) {
      fprintf (stderr, "Reference does not exist: '%s'\n\n", parameters->reference_path.c_str());
      VerboseShortHelpAndExit(argc, argv);
  }

  // Sanity check for the reads path.
  if (argparser.GetArgumentByLongName("reads")->is_set == false && parameters->calc_only_index == false) {
    fprintf (stderr, "Please specify the path to the reads file.\n");
    fprintf (stderr, "\n");
    VerboseShortHelpAndExit(argc, argv);
  }
  if (!fileExists(parameters->reads_path.c_str())) {
    fprintf (stderr, "Reads file does not exist: '%s'\n\n", parameters->reads_path.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

  // Check if the index path was specified, if not, generate it.
  if (argparser.GetArgumentByLongName("index")->is_set == false) {
    parameters->index_file = parameters->reference_path + std::string(".gmidx");
  }

  // Check if the input format is valid.
  if (parameters->infmt != "auto" && parameters->infmt != "fastq" &&
      parameters->infmt != "fasta" && parameters->infmt != "gfa" && parameters->infmt != "sam") {
    fprintf (stderr, "Unknown input format '%s'!\n\n", parameters->infmt.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

  // Check if the output format is valid.
  if (parameters->outfmt != "sam" && parameters->outfmt != "m5" &&
      parameters->outfmt != "mhap" && parameters->outfmt != "paf") {
    fprintf (stderr, "Unknown output format '%s'!\n\n", parameters->outfmt.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

  if (parameters->alignment_algorithm != "sg" && parameters->alignment_algorithm != "sggotoh" &&
      parameters->alignment_algorithm != "anchor" && parameters->alignment_algorithm != "anchorgotoh" && parameters->alignment_algorithm != "anchormex") {
    fprintf (stderr, "Unknown alignment algorithm '%s'!\n\n", parameters->alignment_algorithm.c_str());
    VerboseShortHelpAndExit(argc, argv);
  }

#ifndef RELEASE_VERSION
  if (parameters->debug_read >= 0 || parameters->debug_read_by_qname != "") {
    parameters->verbose_level = 9;
    parameters->num_threads = 1;
    parameters->num_reads_to_process = 1;
  }
#endif

  /// Write this out for every debug verbose level.
  if (parameters->verbose_level > 5) {
    fprintf (stderr, "%s\n", argparser.VerboseArguments().c_str());
  }
//  VerboseProgramParameters(parameters);

  return 0;
}

void VerboseShortHelpAndExit(int argc, char **argv) {
  fprintf (stderr, "For detailed help, please run with -h option.\n");
//  fprintf (stderr, "  %s -h\n", argv[0]);
  fprintf (stderr, "\n");
  fprintf (stderr, "Example usage:\n");
  fprintf (stderr, "  ./graphmap align -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "%s\n", LICENCE_INFORMATION);
  fprintf (stderr, "Version: %s\n", std::string(GRAPHMAP_CURRENT_VERSION).c_str());
  fprintf (stderr, "Build date: %s\n", std::string(GRAPHMAP_CURRENT_VERSION_RELEASE_DATE).c_str());
  fprintf (stderr, "\n");

  exit(0);
}

void VerboseProgramParameters(ProgramParameters *parameters) {
  std::string line_prefix = "===| ";

  fprintf (stderr, "____________________________________\n");
  fprintf (stderr, "Program parameters:\n");
  fprintf (stderr, "Command line: %s\n", parameters->command_line.c_str());
  fprintf (stderr, "%snum_threads = %ld\n", line_prefix.c_str(), parameters->num_threads);
  fprintf (stderr, "%sk_graph = %ld\n", line_prefix.c_str(), parameters->k_graph);
  fprintf (stderr, "%snum_links = %ld\n", line_prefix.c_str(), parameters->num_links);
  fprintf (stderr, "%smin_num_anchor_bases = %ld\n", line_prefix.c_str(), parameters->min_num_anchor_bases);
  fprintf (stderr, "%skmer_step = %ld\n", line_prefix.c_str(), parameters->kmer_step);
  fprintf (stderr, "%salignment_algorithm = %s\n", line_prefix.c_str(), parameters->alignment_algorithm.c_str());
  fprintf (stderr, "%salignment_approach = %s\n", line_prefix.c_str(), parameters->alignment_approach.c_str());
  fprintf (stderr, "%smatch_score = %ld\n", line_prefix.c_str(), parameters->match_score);
  fprintf (stderr, "%smismatch_penalty = %ld\n", line_prefix.c_str(), parameters->mismatch_penalty);
  fprintf (stderr, "%sgap_open_penalty = %ld\n", line_prefix.c_str(), parameters->gap_open_penalty);
  fprintf (stderr, "%sgap_extend_penalty = %ld\n", line_prefix.c_str(), parameters->gap_extend_penalty);
  fprintf (stderr, "%sis_reference_circular = %s\n", line_prefix.c_str(), ((parameters->is_reference_circular == false) ? ("false") : ("true")));

  fprintf (stderr, "%serror_rate = %f\n", line_prefix.c_str(), parameters->error_rate);
  fprintf (stderr, "%sstart_read = %ld\n", line_prefix.c_str(), parameters->start_read);
  fprintf (stderr, "%snum_reads_to_process = %ld\n", line_prefix.c_str(), parameters->num_reads_to_process);
  fprintf (stderr, "%smax_num_hits = %ld\n", line_prefix.c_str(), parameters->max_num_hits);
  fprintf (stderr, "%smax_num_regions_cutoff = %ld\n", line_prefix.c_str(), parameters->max_num_regions_cutoff);
  fprintf (stderr, "%smax_num_regions = %ld\n", line_prefix.c_str(), parameters->max_num_regions);
  fprintf (stderr, "%sreference_path = %s\n", line_prefix.c_str(), parameters->reference_path.c_str());
  fprintf (stderr, "%sindex_reference_path = %s\n", line_prefix.c_str(), parameters->index_file.c_str());
  fprintf (stderr, "%sreads_path = %s\n", line_prefix.c_str(), parameters->reads_path.c_str());
  fprintf (stderr, "%sout_sam_path = %s\n", line_prefix.c_str(), parameters->out_sam_path.c_str());
  fprintf (stderr, "%sreads_folder = %s\n", line_prefix.c_str(), parameters->reads_folder.c_str());
  fprintf (stderr, "%soutput_folder = %s\n", line_prefix.c_str(), parameters->output_folder.c_str());
  fprintf (stderr, "%sverbose_level = %ld\n", line_prefix.c_str(), parameters->verbose_level);
  fprintf (stderr, "%soutput_in_original_order = %s\n", line_prefix.c_str(), (parameters->output_in_original_order == true)?"true":"false");
  fprintf (stderr, "%sprocess_reads_from_folder = %s\n", line_prefix.c_str(), (parameters->process_reads_from_folder == true)?"true":"false");
  fprintf (stderr, "%sbatch_size_in_mb = %ld\n", line_prefix.c_str(), (parameters->batch_size_in_mb));

  fprintf (stderr, "%sdebug_read = %ld\n", line_prefix.c_str(), parameters->debug_read);
  fprintf (stderr, "%sdebug_read_by_qname = %s\n", line_prefix.c_str(), parameters->debug_read_by_qname.c_str());
  fprintf (stderr, "%sverbose_sam_output = %ld\n", line_prefix.c_str(), parameters->verbose_sam_output);
  fprintf (stderr, "%sskip_multiple_kmers_per_bin = %s\n", line_prefix.c_str(), (parameters->skip_multiple_kmers_per_bin == true)?"true":"false");
  fprintf (stderr, "%scomposite_parameters = %s\n", line_prefix.c_str(), parameters->composite_parameters.c_str());

  fprintf (stderr, "%smargin_for_ambiguity = %f\n", line_prefix.c_str(), parameters->margin_for_ambiguity);
  fprintf (stderr, "%soutput_multiple_alignments = %s\n", line_prefix.c_str(), (parameters->output_multiple_alignments == true)?"true":"false");

  fprintf (stderr, "%ssensitive_mode = %s\n", line_prefix.c_str(), (parameters->sensitive_mode == true)?"true":"false");

  fprintf (stderr, "%sevalue_threshold = %f\n", line_prefix.c_str(), parameters->evalue_threshold);
  fprintf (stderr, "%smapq_threshold = %ld\n", line_prefix.c_str(), parameters->mapq_threshold);

  fprintf (stderr, "%soutput_format = '%s'\n", line_prefix.c_str(), parameters->outfmt.c_str());

  fprintf (stderr, "%scalc_only_index = %s\n", line_prefix.c_str(), (parameters->calc_only_index == true)?"true":"false");

  fprintf (stderr, "____________________________________\n");

  fflush(stderr);
}
