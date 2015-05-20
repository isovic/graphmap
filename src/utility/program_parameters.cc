/*
 * program_parameters.cc
 *
 *  Created on: Jul 24, 2014
 *      Author: ivan
 */

#include "utility/program_parameters.h"
#include "utility/argumentparser.h"

int ProcessArgs1(int argc, char **argv, ProgramParameters *parameters)
{
  if (argc == 1) {
    VerboseShortHelpAndExit(argc, argv);
  }

  std::stringstream ss_command_line;
  for (int i=0; i<argc; i++) {
    if (i > 0)
      ss_command_line << " ";
    ss_command_line << argv[i];
  }
  parameters->command_line = ss_command_line.str();

  ArgumentParser argparser;
//  argparser.AddArgument("o", "open", VALUE_TYPE_STRING, "", "Long description of the parameter. If too long, it will automatically be wrapped to 120 characters line width.", 0, "Basic options");
//  argparser.AddArgument("a", "", VALUE_TYPE_NONE, "", "Test for only a short argument.", 0, "Basic options");
//  argparser.AddArgument("", "threads", VALUE_TYPE_INT, "-1", "Test for only a long argument.", 0, "Basic options");
//  argparser.AddArgument("", "", VALUE_TYPE_NONE, "", "Test for specifying no arguments.", 0, "Basic options");
//  argparser.AddArgument("s", "start", VALUE_TYPE_NONE, "", "Starts something important.", 0, "Control options");
//  argparser.AddArgument("", "reads1", VALUE_TYPE_STRING, "reads.fasta", "Path to the file with read sequences.", -2, "Input/Output");
//  argparser.AddArgument("", "reads2", VALUE_TYPE_STRING, "reads2.fasta", "Path to the file with read sequences.", -1, "Input/Output");

  argparser.AddArgument("r", "reference", VALUE_TYPE_STRING,   "", "Path to the reference sequence (fastq or fasta).", 0, "Input/Output options");
  argparser.AddArgument("i", "index",     VALUE_TYPE_STRING,   "", "Path to the index of the reference sequence. If not specified, index is generated in the same path as the reference file, with .gmidx extension.", 0, "Input/Output options");
  argparser.AddArgument("d", "reads",     VALUE_TYPE_STRING,   "", "Path to the reads file (fastq or fasta).", 0, "Input/Output options");
  argparser.AddArgument("o", "outsam",    VALUE_TYPE_STRING,   "", "Path to the output SAM file that will be generated.", 0, "Input/Output options");
  argparser.AddArgument("D", "readsfldr", VALUE_TYPE_STRING,   "", "Path to a folder containing read files (in fastq or fasta format) to process. Cannot be used in combination with '-d' or '-o'.", 0, "Input/Output options");
  argparser.AddArgument("O", "outfldr",   VALUE_TYPE_STRING,   "", "Path to a folder for placing SAM alignments. Use in combination with '-D'.", 0, "Input/Output options");
  argparser.AddArgument("I", "idxonly",   VALUE_TYPE_NONE, "0", "Build only the index from the given reference and exit. If not specified, index will automatically be built if it does not exist, or loaded from file otherwise.", 0, "Input/Output options");
  argparser.AddArgument("u", "keeporder", VALUE_TYPE_NONE, "0", "SAM alignments will be output after the processing has finished, in the order of input reads.", 0, "Input/Output options");
  argparser.AddArgument("B", "batchsize", VALUE_TYPE_INT,   DEFAULT_BATCH_SIZE_IN_MB, "Reads will be loaded in batches of the size specified in megabytes. Value <= 0 loads the entire file.", 0, "Input/Output options");

  argparser.AddArgument("j", "kregion",     VALUE_TYPE_INT,   DEFUALT_K_REGION, "Region selection kmer size.", 0, "Algorithmic options");
  argparser.AddArgument("k", "kgraph",      VALUE_TYPE_INT,   DEFAULT_K_GRAPH, "Graph construction kmer size.", 0, "Algorithmic options");
  argparser.AddArgument("l", "nlinks",      VALUE_TYPE_INT,   DEFAULT_NUM_LINKS, "Number of edges per vertex.", 0, "Algorithmic options");
  argparser.AddArgument("e", "errrate",     VALUE_TYPE_FLOAT, DEFAULT_ERROR_RATE, "Approximate error rate of the input read sequences.", 0, "Algorithmic options");
  argparser.AddArgument("m", "maxhits",     VALUE_TYPE_INT,   DEFAULT_MAX_NUM_HITS, "Maximum number of hits per kmer. If 0, threshold will be estimated automatically. If < 0, all hits will be taken into account.", 0, "Algorithmic options");
  argparser.AddArgument("g", "maxreg",      VALUE_TYPE_INT,   DEFAULT_MAX_NUM_REGIONS, "If the final number of regions exceeds this amount, the read will be called unmapped. If 0, value will be dynamically determined. If < 0, no limit is set.", 0, "Algorithmic options");
  argparser.AddArgument("q", "regcutoff",   VALUE_TYPE_INT,   DEFAULT_MAX_NUM_REGIONS_CUTOFF, "Attempt to heuristically reduce the number of regions if it exceeds this amount. Value <= 0 disables reduction but only if param -g is not 0. If -g is 0, the value of this parameter is set to 1/5 of maximum number of regions.", 0, "Algorithmic options");
  argparser.AddArgument("f", "noedfilt",    VALUE_TYPE_NONE, "0", "Disables filtering by edit distance in the final step of mapping.", 0, "Algorithmic options");
  argparser.AddArgument("C", "circular",    VALUE_TYPE_NONE, "0", "Reference sequence is a circular genome.", 0, "Algorithmic options");
  argparser.AddArgument("F", "ambigdiff", VALUE_TYPE_FLOAT,   DEFAULT_MARGIN_FOR_AMBIGUITY, "All mapping positions within the given fraction of the top score will be counted for ambiguity (mapping quality). Value of 0.0f accounts only for identical mappings.", 0, "Algorithmic options");
  argparser.AddArgument("Z", "outall", VALUE_TYPE_NONE, "0", "If specified, all alignments within (-P FLT) will be output. Otherwise, single best alignment will be output.", 0, "Algorithmic options");
  argparser.AddArgument("P", "parsimon",    VALUE_TYPE_NONE, "0", "If specified, the parsimonoius memory mode will be used. If omitted, a fast and sensitive (but 2x memory-hungry) mode will be used.", 0, "Algorithmic options");

#ifndef RELEASE_VERSION
  argparser.AddArgument("M", "match", VALUE_TYPE_INT,   DEFAULT_MATCH_SCORE, "Match score for the DP alignment.", 0, "Algorithmic options");
  argparser.AddArgument("X", "mismatch", VALUE_TYPE_INT,   DEFAULT_MISMATCH_PENALTY, "Mismatch penalty for the DP alignment.", 0, "Algorithmic options");
  argparser.AddArgument("G", "gapopen", VALUE_TYPE_INT,   DEFAULT_GAP_OPEN_PENALTY, "Gap open penalty for the DP alignment.", 0, "Algorithmic options");
  argparser.AddArgument("E", "gapext", VALUE_TYPE_INT,   DEFAULT_GAP_EXTEND_PENALTY, "Gap extend penalty for the DP alignment.", 0, "Algorithmic options");
  argparser.AddArgument("S", "kstep", VALUE_TYPE_INT,   DEFAULT_KMER_STEP, "Kmer step for region selection.", 0, "Algorithmic options");
#endif

#ifndef RELEASE_VERSION
  argparser.AddArgument("y", "", VALUE_TYPE_INT,   DEFAULT_DEBUG_READ, "ID of the read to give the detailed verbose output.", 0, "Debug options");
  argparser.AddArgument("Y", "", VALUE_TYPE_STRING,   "", "QNAME of the read to give the detailed verbose output. Has precedence over -y.", 0, "Debug options");
  argparser.AddArgument("b", "verbosam", VALUE_TYPE_INT,   DEFAULT_VERBOSE_SAM_OUTPUT, "Helpful debug comments can be placed in SAM output lines (at the end). Comments can be turned of by setting this parameter to 0. Different values increase/decrease verbosity level.", 0, "Debug options");
#endif

  argparser.AddArgument("x", "preset", VALUE_TYPE_STRING,   "", "Pre-set parameters to increase sensitivity for different sequencing technologies. Valid options are: 'illumina' and 'nanopore'.", 0, "General-purpose pre-set options");

  argparser.AddArgument("t", "threads", VALUE_TYPE_INT,   DEFAULT_NUM_THREADS, "Number of threads to use. If '-1', number of threads will be equal to the number of cores.", 0, "Other options");
  argparser.AddArgument("v", "verbose", VALUE_TYPE_INT,   DEFAULT_VERBOSE_LEVEL, "Verbose level. If equal to 0 nothing except strict output will be placed on stdout.", 0, "Other options");
  argparser.AddArgument("s", "startread", VALUE_TYPE_INT,   DEFAULT_START_READ, "Ordinal number of the read from which to start processing data.", 0, "Other options");
  argparser.AddArgument("n", "numreads", VALUE_TYPE_INT,   DEFAULT_NUM_READS_TO_PROCESS, "Number of reads to process per batch. Value of '-1' processes all reads.", 0, "Other options");
  argparser.AddArgument("h", "help", VALUE_TYPE_NONE, "0", "View this help.", 0, "Other options");

  argparser.ProcessArguments(argc, argv);

  fprintf (stderr, "%s", argparser.VerboseArgumentsByGroup().c_str());
  fprintf (stderr, "\n");
  fprintf (stderr, "Results of parsing the command line:\n");

  argparser.VerboseArguments(stdout);

  exit(1);

  return 0;
}

int ProcessArgs(int argc, char **argv, ProgramParameters *parameters)
{
  if (argc == 1) {
    VerboseShortHelpAndExit(argc, argv);
  }

  std::stringstream ss_command_line;
  for (int i=0; i<argc; i++) {
    if (i > 0)
      ss_command_line << " ";
    ss_command_line << argv[i];
  }

  parameters->command_line = ss_command_line.str();

  int c;
  int count_necessary_parameters = 0;
  std::string machine_type = "";

  opterr = 0;

  bool reads_specified_by_file = false;
  bool reads_specified_by_folder = false;
  bool output_specified_by_file = false;
  bool output_specified_by_folder = false;

  while ((c = getopt (argc, argv, "k:l:e:s:n:y:Y:t:r:i:d:o:b:v:g:hx:a:w:uq:D:O:B:IG:E:M:X:CF:ZS:PA:")) != -1) {
    switch (c) {
//      case 'j':
//        sscanf (optarg, "%ld", &(parameters->k_region));
//        break;
      case 'A':
        sscanf (optarg, "%ld", &(parameters->min_num_anchor_bases));
        break;
      case 'k':
        sscanf (optarg, "%ld", &(parameters->k_graph));
        break;
      case 'l':
        sscanf (optarg, "%ld", &(parameters->num_links));
        break;
      case 'e':
        sscanf (optarg, "%f", &(parameters->error_rate));
        break;
      case 's':
        sscanf (optarg, "%ld", &(parameters->start_read));
        break;
      case 'n':
        sscanf (optarg, "%ld", &(parameters->num_reads_to_process));
        break;
      case 'y':
        sscanf (optarg, "%ld", &(parameters->debug_read));
        break;
      case 'Y':
        parameters->debug_read_by_qname = std::string(optarg);
        if (parameters->debug_read_by_qname.size() > 2 && ((parameters->debug_read_by_qname.front())) == '"' && ((parameters->debug_read_by_qname.back())) == '"') {
          parameters->debug_read_by_qname = parameters->debug_read_by_qname.substr(1, (parameters->debug_read_by_qname.size() - 2));
        }
        break;
      case 't':
        sscanf (optarg, "%ld", &(parameters->num_threads));
        break;
      case 'r':
        parameters->reference_path = std::string(optarg);
        count_necessary_parameters += 1;
        break;
      case 'i':
        parameters->index_reference_path = std::string(optarg);
        break;
      case 'I':
        parameters->calc_only_index = true;
        break;
      case 'd':
        if (reads_specified_by_folder == true) {
          fprintf (stderr, "Already specified path to the folder containing reads for batch processing. Can not also specify an exact reads file using option '-d'.\n\n");
//          VerboseUsage(stdout);
          VerboseShortHelpAndExit(argc, argv);
        }
        reads_specified_by_file = true;

        parameters->reads_path = std::string(optarg);
        count_necessary_parameters += 1;
        break;
      case 'D':
        if (reads_specified_by_file == true) {
          fprintf (stderr, "Already specified path to the exact FASTA/FASTQ input file. Can not also specify a folder for batch input file processing using option '-D'.\n\n");
//          VerboseUsage(stdout);
          VerboseShortHelpAndExit(argc, argv);
        }
        reads_specified_by_folder = true;

        parameters->reads_folder = std::string(optarg);
        count_necessary_parameters += 1;
        break;

      case 'o':
        if (output_specified_by_folder == true) {
          fprintf (stderr, "Already specified path to the output folder for batch processing. Can not also specify an exact SAM output file using option '-o'.\n\n");
//          VerboseUsage(stdout);
          VerboseShortHelpAndExit(argc, argv);
        }
        output_specified_by_file = true;

        parameters->out_sam_path = std::string(optarg);
//        count_necessary_parameters += 1;
        break;
      case 'O':
        if (output_specified_by_file == true) {
          fprintf (stderr, "Already specified path to the exact SAM output file. Can not also specify an output folder for batch file processing using option '-O'.\n\n");
//          VerboseUsage(stdout);
          VerboseShortHelpAndExit(argc, argv);
        }
        output_specified_by_folder = true;

        parameters->output_folder = std::string(optarg);
//        count_necessary_parameters += 1;
        break;

      case 'b':
        sscanf (optarg, "%ld", &(parameters->verbose_sam_output));
        break;
      case 'v':
        sscanf (optarg, "%ld", &(parameters->verbose_level));
        break;
//      case 'm':
//        sscanf (optarg, "%ld", &(parameters->max_num_hits));
//        break;
      case 'g':
        sscanf (optarg, "%ld", &(parameters->max_num_regions));
        break;

      case 'C':
        parameters->is_reference_circular = true;
        break;

      case 'F':
        sscanf (optarg, "%f", &(parameters->margin_for_ambiguity));
        break;
      case 'Z':
//        sscanf (optarg, "%d", &(parameters->output_multiple_alignments));
        parameters->output_multiple_alignments = 1;
        break;

      case 'P':
        parameters->parsimonious_mode = !DEFAULT_PARSIMONIOUS_MODE;
        break;

      case 'x':
//#ifndef RELEASE_VERSION
        if (optind != 3) {      // This points to the parameter after the optopt (the parameter after the value of -x parameter).
          fprintf (stderr, "Composite parameter '-x' needs to be specified first in order to avoid wrong settings.\n\n");
//          printf ("%ld\n", optind);
//          for (int i=0; i<argc; i++)
//            printf ("%s\n", argv[i]);

          VerboseShortHelpAndExit(argc, argv);
        }

        machine_type = ((char *) optarg);;

        if (machine_type == ((std::string) "illumina")) {
//          parameters->k_region = 13;
//          parameters->k_graph = 6;
//          parameters->num_links = 9;
//          parameters->error_rate = 0.45f;

//          parameters->start_read = 0;
//          parameters->num_reads_to_process = -1;
//          parameters->debug_read = -1;
//          parameters->debug_read_by_qname = "";

//          parameters->match_score = 1;  // 1
//          parameters->mismatch_penalty = 1; // 3;
//          parameters->gap_open_penalty = 8; // 6
//          parameters->gap_extend_penalty = 2;
          parameters->match_score = 5;  // 1
          parameters->mismatch_penalty = 4; // 3;
          parameters->gap_open_penalty = 8; // 6
          parameters->gap_extend_penalty = 6;

//          parameters->max_num_hits = 0;
//          parameters->max_num_regions = 0;
//          parameters->max_num_regions_cutoff = 0;
//          parameters->filter_by_edit_distance = true;
          parameters->alignment_algorithm = "gotoh";
          parameters->composite_parameters = machine_type;
        }
//        else if (machine_type == ((std::string) "pacbio")) {
//          parameters->k_region = 13;
//          parameters->k_graph = 6; // 7;  // 9
//          parameters->num_links = 18; // 21;  // 17
//          parameters->error_rate = 0.45f; // 0.20;  // 0.18
//          parameters->start_read = 0;
//          parameters->num_reads_to_process = -1;
//          parameters->debug_read = -1;
//          parameters->match_score = 1;
//          parameters->mismatch_penalty = 1;
//          parameters->gap_open_penalty = 1;
//          parameters->gap_extend_penalty = 1;
//          parameters->max_num_hits = 0;
//          parameters->max_num_regions = 0;
//          parameters->max_num_regions_cutoff = 0;
//          parameters->filter_by_edit_distance = true;
//        }
//        else if (machine_type == ((std::string) "nanopore")) {
////          parameters->k_region = 11;  // Commented out on 18.11.2014.
//          parameters->k_region = 13;
//          parameters->k_graph = 6;
//          parameters->num_links = 18;
//          parameters->error_rate = 0.45f;
//          parameters->start_read = 0;
//          parameters->num_reads_to_process = -1;
//          parameters->debug_read = -1;
//          parameters->debug_read_by_qname = "";
//          parameters->match_score = 2;
//          parameters->mismatch_penalty = 1;
//          parameters->gap_open_penalty = 2;
//          parameters->gap_extend_penalty = 1;
//          parameters->max_num_hits = 0;
//          parameters->max_num_regions = 0;
//          parameters->max_num_regions_cutoff = 0;
//          parameters->filter_by_edit_distance = true;
//          parameters->alignment_algorithm = "myers";
//          parameters->composite_parameters = machine_type;
//        }
        else {
//          printf ("Unknown option: '%s'. Valid options are: 'illumina', 'pacbio' and 'nanopore'. Using default parameters.\n", optarg);
          fprintf (stderr, "Unknown option: '%s'. Valid options are: 'illumina' and 'nanopore'. Skipping parameter.\n", optarg);
        }

//#endif

        break;

//      case 'p':
//        parameters->skip_multiple_kmers_per_bin = false;
//        break;
//      case 'f':
//        parameters->filter_by_edit_distance = false;
//        break;
      case 'h':
        VerboseUsageAndExit(stdout);
        break;

      case 'u':
        parameters->output_in_original_order = true;
        break;
      case 'q':
        sscanf (optarg, "%ld", &(parameters->max_num_regions_cutoff));
        break;

      case 'a':
        parameters->alignment_algorithm = std::string(optarg);
        break;
      case 'w':
        parameters->alignment_approach = std::string(optarg);
        break;
      case 'G':
        sscanf (optarg, "%ld", &(parameters->gap_open_penalty));
        break;
      case 'E':
        sscanf (optarg, "%ld", &(parameters->gap_extend_penalty));
        break;
      case 'M':
        sscanf (optarg, "%ld", &(parameters->match_score));
        break;
      case 'X':
        sscanf (optarg, "%ld", &(parameters->mismatch_penalty));
        break;

#ifndef RELEASE_VERSION
      case 'S':
        sscanf (optarg, "%d", &(parameters->kmer_step));
        break;

#endif

      case 'B':
        sscanf (optarg, "%ld", &(parameters->batch_size_in_mb));
        break;

//      case 'A':
//        VerboseEula(stderr);
//        break;

      case '?':
        fprintf (stderr, "Option -%c requires an argument or an unknown option occured!\n", optopt);

        return 1;
      default:
        fprintf (stderr, "Unknown parameter: '-%c'.\n\n", c);
        VerboseShortHelpAndExit(argc, argv);
//        continue;
//        abort();
    }
  }

  // Sanity check for the reference path.
  if (parameters->reference_path == "") {
    fprintf (stderr, "Please specify the path to the reference file.\n");
    fprintf (stderr, "\n");
    VerboseShortHelpAndExit(argc, argv);
  }

  // Check if the index path was specified, if not, generate it.
  if (parameters->index_reference_path.size() == 0)
    parameters->index_reference_path = parameters->reference_path + std::string(".gmidx");
//    parameters->index_reference_path = parameters->reference_path + std::string(".index");



  // If the -I option was not specified (calculate only index), then the mapping process
  // will be used. Check if the input paths are specified.
  if (parameters->calc_only_index == false) {
    if (parameters->reads_path == "" && parameters->reads_folder == "") {
      fprintf (stderr, "Input reads file/folder not specified!\n\n");
      VerboseShortHelpAndExit(argc, argv);
    }

    if (reads_specified_by_file == true && reads_specified_by_folder == false) // && output_specified_by_file == true)
      parameters->process_reads_from_folder = false;
    else if (reads_specified_by_folder == true && output_specified_by_folder == true && reads_specified_by_file == false) //)
      parameters->process_reads_from_folder = true;
    else {
      fprintf (stderr, "Input is not correctly specified. Either a single FASTA/FASTQ needs to be specified with the '-d' parameter,\n");
      fprintf (stderr, "or input and output folders need to be specified ('-D' and '-O' parameters), but not both.\n\n");
      VerboseShortHelpAndExit(argc, argv);
    }

    if (count_necessary_parameters != 2) {
      fprintf (stderr, "Please specify required parameters. These include: '-r' and '-d'; or '-r', '-D' and '-O'.\n");
      fprintf (stderr, "\n");
      VerboseShortHelpAndExit(argc, argv);
    }
  }

  #ifndef RELEASE_VERSION
    if (parameters->verbose_level > 0) {
      VerboseProgramParameters(parameters);
    }
  #endif

  #ifdef RELEASE_VERSION
    if (parameters->verbose_level < 0 || parameters->verbose_level > 5) {
      parameters->verbose_level = 5;
    }

    parameters->debug_read = -1;
    parameters->debug_read_by_qname = "";
//    parameters->skip_multiple_kmers_per_bin = true;

    parameters->verbose_sam_output = 0;

    parameters->kmer_step = 1;
//    parameters->realignment_algorithm = "edlib";
  #endif

  return 0;
}

void VerboseShortHelpAndExit(int argc, char **argv) {
//  fprintf (stderr, "Licence for this instance of GraphMap: free for non-commercial, non-profit, academic use only.\n");
//  fprintf (stderr, "To view the EULA, please run with the -A option.\n");
//  fprintf (stderr, "  %s -A\n", argv[0]);
//  fprintf (stderr, "\n");

  fprintf (stderr, "For detailed help, please run with -h option.\n");
  fprintf (stderr, "  %s -h\n", argv[0]);
//  fprintf (stderr, "\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Example usage:\n");
  fprintf (stderr, "  ./graphmap -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n");
  fprintf (stderr, "\n");
//  fprintf (stderr, "\n");

  fprintf (stderr, "%s\n", LICENCE_INFORMATION);

//  fprintf (stderr, "GraphMap (c) by Ivan Sovic, Mile Sikic and Niranjan Nagarajan\n");
////  fprintf (stderr, "\n");
//  fprintf (stderr, "GraphMap is licensed under a\n");
//  fprintf (stderr, "Creative Commons Attribution-NonCommercial 4.0 International License.\n");
////  fprintf (stderr, "\n");
//  fprintf (stderr, "You should have received a copy of the license along with this\n");
//  fprintf (stderr, "work. If not, see <http://creativecommons.org/licenses/by-nc/4.0/>.\n");
//  fprintf (stderr, "\n");
  exit(0);
}

void VerboseProgramParameters(ProgramParameters *parameters) {
  std::string line_prefix = "===| ";

  #ifndef RELEASE_VERSION
  printf ("____________________________________\n");
  printf ("Program parameters:\n");
  printf ("Command line: %s\n", parameters->command_line.c_str());
  printf ("%snum_threads = %ld\n", line_prefix.c_str(), parameters->num_threads);
//  printf ("%sk_region = %ld\n", line_prefix.c_str(), parameters->k_region);
  printf ("%sk_graph = %ld\n", line_prefix.c_str(), parameters->k_graph);
  printf ("%snum_links = %ld\n", line_prefix.c_str(), parameters->num_links);
  printf ("%smin_num_anchor_bases = %ld\n", line_prefix.c_str(), parameters->min_num_anchor_bases);
  printf ("%skmer_step = %ld\n", line_prefix.c_str(), parameters->kmer_step);
  printf ("%salignment_algorithm = %s\n", line_prefix.c_str(), parameters->alignment_algorithm.c_str());
  printf ("%salignment_approach = %s\n", line_prefix.c_str(), parameters->alignment_approach.c_str());
  printf ("%smatch_score = %ld\n", line_prefix.c_str(), parameters->match_score);
  printf ("%smismatch_penalty = %ld\n", line_prefix.c_str(), parameters->mismatch_penalty);
  printf ("%sgap_open_penalty = %ld\n", line_prefix.c_str(), parameters->gap_open_penalty);
  printf ("%sgap_extend_penalty = %ld\n", line_prefix.c_str(), parameters->gap_extend_penalty);
  printf ("%sis_reference_circular = %s\n", line_prefix.c_str(), ((parameters->is_reference_circular == false) ? ("false") : ("true")));

  printf ("%serror_rate = %f\n", line_prefix.c_str(), parameters->error_rate);
  printf ("%sstart_read = %ld\n", line_prefix.c_str(), parameters->start_read);
  printf ("%snum_reads_to_process = %ld\n", line_prefix.c_str(), parameters->num_reads_to_process);
  printf ("%smax_num_hits = %ld\n", line_prefix.c_str(), parameters->max_num_hits);
  printf ("%smax_num_regions_cutoff = %ld\n", line_prefix.c_str(), parameters->max_num_regions_cutoff);
  printf ("%smax_num_regions = %ld\n", line_prefix.c_str(), parameters->max_num_regions);
  printf ("%sreference_path = %s\n", line_prefix.c_str(), parameters->reference_path.c_str());
  printf ("%sindex_reference_path = %s\n", line_prefix.c_str(), parameters->index_reference_path.c_str());
  printf ("%sreads_path = %s\n", line_prefix.c_str(), parameters->reads_path.c_str());
  printf ("%sout_sam_path = %s\n", line_prefix.c_str(), parameters->out_sam_path.c_str());
  printf ("%sreads_folder = %s\n", line_prefix.c_str(), parameters->reads_folder.c_str());
  printf ("%soutput_folder = %s\n", line_prefix.c_str(), parameters->output_folder.c_str());
  printf ("%sverbose_level = %ld\n", line_prefix.c_str(), parameters->verbose_level);
  printf ("%soutput_in_original_order = %s\n", line_prefix.c_str(), (parameters->output_in_original_order == true)?"true":"false");
  printf ("%sprocess_reads_from_folder = %s\n", line_prefix.c_str(), (parameters->process_reads_from_folder == true)?"true":"false");
  printf ("%sfilter_by_edit_distance = %s\n", line_prefix.c_str(), (parameters->filter_by_edit_distance == true)?"true":"false");
  printf ("%sbatch_size_in_mb = %ld\n", line_prefix.c_str(), (parameters->batch_size_in_mb));

  printf ("%sdebug_read = %ld\n", line_prefix.c_str(), parameters->debug_read);
  printf ("%sdebug_read_by_qname = %s\n", line_prefix.c_str(), parameters->debug_read_by_qname.c_str());
  printf ("%sverbose_sam_output = %ld\n", line_prefix.c_str(), parameters->verbose_sam_output);
  printf ("%sskip_multiple_kmers_per_bin = %s\n", line_prefix.c_str(), (parameters->skip_multiple_kmers_per_bin == true)?"true":"false");
  printf ("%scomposite_parameters = %s\n", line_prefix.c_str(), parameters->composite_parameters.c_str());

  printf ("%smargin_for_ambiguity = %f\n", line_prefix.c_str(), parameters->margin_for_ambiguity);
  printf ("%soutput_multiple_alignments = %ld\n", line_prefix.c_str(), parameters->output_multiple_alignments);

  printf ("%sparsimonious_mode = %s\n", line_prefix.c_str(), (parameters->parsimonious_mode == true)?"true":"false");

  printf ("____________________________________\n");
  #endif

  fflush(stdout);
}



void VerboseUsageAndExit(FILE *fp) {
  std::stringstream ss;

  ss << SOFTWARE_NAME << " - A very accurate and sensitive long-read, high error-rate sequence mapper\n", SOFTWARE_NAME;
  ss << "Version: " <<  GRAPHMAP_CURRENT_VERSION << "\n";
  ss << "Build date: " <<  std::string(GRAPHMAP_CURRENT_VERSION_RELEASE_DATE).c_str() << "\n";
//  ss << COPYRIGHT << "\n";
  ss << "\n";

  ss << LICENCE_INFORMATION << "\n";
  ss << AFFILIATIONS << "\n";

  ss << "\n";

  ss << "Usage:\n";
  ss << "\tgraphmap [options] -r <reference_file> -d <reads_file> -o <output_sam_path>\n";
  ss << "\n";

  ss << "Input/Output options:\n";
  ss << "\t-r STR\tPath to the reference sequence (fastq or fasta).\n";
  ss << "\t-i STR\tPath to the index of the reference sequence. If not specified, index is generated\n\t\tin the same path as the reference file, with .gmidx extension.\n\t\tFor non-parsimonious mode, secondary index .gmidxsec is also generated.\n";
  ss << "\t-d STR\tPath to the reads file (fastq or fasta).\n";
  ss << "\t-o STR\tPath to the output SAM file that will be generated.\n";
  ss << "\n";
  ss << "\t-D STR\tPath to a folder containing read files (in fastq or fasta format) to process.\n\t\tCannot be used in combination with '-d' or '-o'.\n";
  ss << "\t-O STR\tPath to a folder for placing SAM alignments. Use in combination with '-D'.\n";
  ss << "\n";
  ss << "\t-I\tBuild only the index from the given reference and exit. If not specified, index will\n\t\tautomatically be built if it does not exist, or loaded from file otherwise.\n";
  ss << "\t-u\tSAM alignments will be output after the processing has finished,\n\t\tin the order of input reads.\n";
  ss << "\t-B INT\tReads will be loaded in batches of the size specified in megabytes.\n\t\tValue <= 0 loads the entire file. [" << DEFAULT_BATCH_SIZE_IN_MB << "]\n";
  ss << "\n";

  ss << "General-purpose pre-set options:\n";
//  ss << "\t-x STR\tPre-set parameters to increase sensitivity for different sequencing technologies.\n\t\tValid options are: 'illumina' and 'nanopore'.\n";
  ss << "\t-x STR\tPre-set parameters to increase sensitivity for different sequencing technologies.\n\t\tValid options are:\n";
  ss << "\t             illumina - Equivalent to: '-a gotoh -w sg -M 5 -X 4 -G 8 -E 6'\n";
  ss << "\n";

//#ifndef RELEASE_VERSION
//  ss << "General-purpose pre-set options:\n";
//  ss << "\t-x STR\tPre-set parameters for different sequencing technologies.\n\t\tValid options are: 'illumina', 'pacbio' and 'nanopore'.\n";
//  ss << "\t\t\tillumina\tEquivalent to: -j 13 -k 6 -l 18 -e 0.02 -f\n";
//  ss << "\t\t\tpacbio\t\tEquivalent to: -j 13 -k 7 -l 21 -e 0.20\n";
//  ss << "\t\t\tnanopore\tEquivalent to: -j 13 -k 6 -l 18 -e 0.45\n";
//  ss << "\n";
//#endif

  ss << "Algorithmic options:\n";
//  ss << "\t-j INT\tRegion selection kmer size. [" << DEFUALT_K_REGION << "]\n";
  ss << "\t-k INT\tGraph construction kmer size. [" <<  DEFAULT_K_GRAPH << "]\n";
  ss << "\t-l INT\tNumber of edges per vertex. [" <<  DEFAULT_NUM_LINKS << "]\n";
  ss << "\t-A INT\tMinimum number of match bases in an anchor. [" << DEFAULT_MIN_NUM_ANCHOR_BASES << "]\n";
  ss << "\t-e FLT\tApproximate error rate of the input read sequences [" <<  DEFAULT_ERROR_RATE << "]\n";
//  ss << "\t-m INT\tMaximum number of hits per kmer. If 0, threshold will be estimated automatically.\n\t\tIf < 0, all hits will be taken into account [" <<  DEFAULT_MAX_NUM_HITS << "]\n";
  ss << "\t-g INT\tIf the final number of regions exceeds this amount, the read will be called unmapped.\n\t\tIf 0, value will be dynamically determined. If < 0, no limit is set. [" << DEFAULT_MAX_NUM_REGIONS << "]\n";
  ss << "\t-q INT\tAttempt to heuristically reduce the number of regions if it exceeds this amount.\n\t\tValue <= 0 disables reduction but only if param -g is not 0.\n\t\tIf -g is 0, the value of this parameter is set to 1/5 of maximum number of regions. [" << DEFAULT_MAX_NUM_REGIONS_CUTOFF << "]\n";
//  ss << "\t-f\tDisables filtering by edit distance in the final step of mapping\n";
  ss << "\t-C\tReference sequence is a circular genome. [" << ((DEFAULT_IS_REFERENCE_CIRCULAR == false) ? ("false") : ("true")) << "]\n";
  ss << "\n";
  ss << "\t-P FLT\tAll mapping positions within the given fraction of the top score will be counted for\n\t\tambiguity (mapping quality). Value of 0.0 counts only identical mappings. [" << DEFAULT_MARGIN_FOR_AMBIGUITY << "]\n";
  ss << "\t-Z\tIf specified, all alignments within (-P FLT) will be output to a file.\n\t\tOtherwise, only one alignment will be output.\n";
  ss << "\n";
  ss << "\t-F\tIf specified, the parsimonious memory mode will be used. If omitted, a fast and sensitive\n\t\t(but 2x memory-hungry) mode will be used.\n";
  ss << "\n";

#ifndef RELEASE_VERSION
    ss << "\t-S INT\tKmer step for region selection. [" << DEFAULT_KMER_STEP << "]\n";
    ss << "\n";
#endif

  ss << "Alignment options:\n";
//    ss << "\t-w INT\tKmer step during graph construction (the number of bases to skip between\n\t\tbeginnings of every adjacent kmer) [" << DEFAULT_KMER_GRAPH_STEP << "]\n";
//    ss << "\t-p\tOne kmer of a read can have multiple hits within the same region.\n\t\tThis parameter enables using multiple hits per region\n";
  ss << "\t-a STR\tSpecifies which algorithm should be used for alignment. Options are: [" << DEFAULT_ALIGNMENT_ALGORITHM << "]\n";
  ss << "\t             myers     - Myers' bit-vector approach. Edit distance alignment.\n";
  ss << "\t             gotoh     - Gotoh alignment with affine gaps.\n";
  ss << "\t-w STR\tSpecifies the alignment strategy. Options are: [" << DEFAULT_ALIGNMENT_APPROACH << "]\n";
  ss << "\t             sg     - semiglobal alignment over best region. Can be used with both Myers and Gotoh.\n";
  ss << "\t             anchor - anchored alignment with end-to-end extension. Uses Myers alignment only.\n";
//  ss << "\t             overlap - anchored alignment with end-to-end extension\n";
//  ss << "\t             splice - spliced alignment (each anchor chain output separately)";

  ss << "\t-M INT\tMatch score for the DP alignment. Ignored for Myers alignment. [" << DEFAULT_MATCH_SCORE << "]\n";
  ss << "\t-X INT\tMismatch penalty for the DP alignment. Ignored for Myers alignment. [" << DEFAULT_MISMATCH_PENALTY << "]\n";
  ss << "\t-G INT\tGap open penalty for the DP alignment. Ignored for Myers alignment. [" << DEFAULT_GAP_OPEN_PENALTY << "]\n";
  ss << "\t-E INT\tGap extend penalty for the DP alignment. Ignored for Myers alignment. [" << DEFAULT_GAP_EXTEND_PENALTY << "]\n";

  ss << "\n";

  ss << "Other options:\n";
  ss << "\t-t INT\tNumber of threads to use. If '-1', number of threads will be equal\n\t\tto min(24, num_cores/2). [" << DEFAULT_NUM_THREADS << "]\n";
  ss << "\t-v INT\tVerbose level. If equal to 0 nothing except strict output will be\n\t\tplaced on stdout. [" << DEFAULT_VERBOSE_LEVEL  << "]\n";
//  ss << "\t-A    \tShow the End User Licence Agreement.\n";
//  #ifndef RELEASE_VERSION
  ss << "\t-s INT\tOrdinal number of the read from which to start processing data [" << DEFAULT_START_READ << "]\n";
  ss << "\t-n INT\tNumber of reads to process per batch. Value of '-1' processes all reads [" << DEFAULT_NUM_READS_TO_PROCESS << "]\n";
  ss << "\t-h    \tView this help.\n";
//  #endif
  ss << "\n";

#ifndef RELEASE_VERSION
    ss << "Debug options:\n";
    ss << "\t-y INT\tID of the read to give the detailed verbose output [" << DEFAULT_DEBUG_READ << "]\n";
    ss << "\t-Y STR\tQNAME of the read to give the detailed verbose output. Has precedence over -y.\n";
    ss << "\t-b INT\tHelpful debug comments can be placed in SAM output lines (at the end).\n\t\tComments can be turned off by setting this parameter to 0.\n\t\tDifferent values increase/decrease verbosity level [" << DEFAULT_VERBOSE_SAM_OUTPUT << "]\n";
    ss << "\n";
#endif

  ss << "Example usage:\n";
  ss << "\t# Process all reads from a given FASTA/FASTQ file:\n";
  ss << "\t./graphmap -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n";
  ss << "\n";
//  ss << "\t# Process reads using more sensitive parameters for Illumina and nanopore data:\n";
  ss << "\t# Process reads using more sensitive alignment parameters for Illumina data:\n";
//  ss << "\t./graphmap -x nanopore -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n";
  ss << "\t./graphmap -x illumina -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n";
  ss << "\n";
  ss << "\t# Process reads from a circular genome:\n";
  ss << "\t./graphmap -C -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n";
//  ss << "\t./graphmap -x nanopore -C -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n";
  ss << "\n";
  ss << "\t# Process all reads from a given folder.\n"; //    \t# Process nanopore reads.\n";
  ss << "\t./graphmap -r escherichia_coli.fa -D reads_folder -O alignments_folder\n";
  ss << "\n";
  ss << "\t# Generate only the index.\n";
  ss << "\t./graphmap -I -r escherichia_coli.fa\n";
  ss << "\n";

//  ss << "Example usage:\n";
//  ss << "\t# Process nanopore reads.\n"; //    \t# Process nanopore reads.\n";
//  ss << "\t./graphmap -x nanopore -r escherichia_coli.fa -d reads.fastq -o alignments.sam\n"; //    \t# Process nanopore reads.\n";
//  ss << "\n";
//  ss << "\t# Process all reads from a given folder.\n"; //    \t# Process nanopore reads.\n";
//  ss << "\t./graphmap -x nanopore -r escherichia_coli.fa -D reads_folder -O alignments_folder\n"; // \t# Process all reads from a given folder.\n";
//  ss << "\n";
//  ss << "\t# Generate only the index.\n";
//  ss << "\t./graphmap -I -r escherichia_coli.fa\n";
//  ss << "\n";

  ss << "Contact:\n";
  ss << "\tivan.sovic@irb.hr, mile.sikic@fer.hr, nagarajann@gis.a-star.edu.sg\n";
  ss << "\n";

  fprintf (fp, "%s", ss.str().c_str());

//  fprintf (fp, "%s - A very accurate and sensitive long-read, high error-rate sequence mapper\n", SOFTWARE_NAME);
//  fprintf (fp, "Version: %s\n", GRAPHMAP_CURRENT_VERSION);
//  fprintf (fp, "Build date: %s\n", std::string(GRAPHMAP_CURRENT_VERSION_RELEASE_DATE).c_str());
//  fprintf (fp, "%s\n", COPYRIGHT);
//
//  fprintf (fp, "%s\n", LICENCE_INFORMATION);
//
//  fprintf (fp, "Usage:\n");
//  fprintf (fp, "\tgraphmap [options] -r <reference_file> -d <reads_file> -o <output_sam_path>\n");
//  fprintf (fp, "\n");
//
//  fprintf (fp, "Input/Output options:\n");
//  fprintf (fp, "\t-r STR\tPath to the reference sequence (fastq or fasta)\n");
//  fprintf (fp, "\t-i STR\tPath to the index of the reference sequence. If the index does not exist, it will be automatically generated, otherwise it will be loaded from disk.\n");
//  fprintf (fp, "\t-d STR\tPath to the reads file (fastq or fasta)\n");
//  fprintf (fp, "\t-o STR\tPath to the output SAM file that will be generated\n");
//  fprintf (fp, "\n");
//  fprintf (fp, "\t-D STR\tPath to a folder containing read files (in fastq or fasta format) to process. Cannot be used in combination with '-d' or '-o'.\n");
//  fprintf (fp, "\t-O STR\tPath to a folder for placing SAM alignments. Use in combination with '-D'\n");
//  fprintf (fp, "\n");
//  fprintf (fp, "\t-u\tSAM alignments will be output after the processing has finished, in the order of input reads\n");
//  fprintf (fp, "\t-B INT\tReads will be loaded in batches of the size specified in megabytes. Value <= 0 loads the entire file [%ld]\n", DEFAULT_BATCH_SIZE_IN_MB);
//  fprintf (fp, "\n");
//
//  fprintf (fp, "General-purpose pre-set options:\n");
//  fprintf (fp, "\t-x STR\tPre-set parameters for different sequencing technologies. Valid options are: 'ngs', 'pacbio' and 'nanopore'.\n");
//  fprintf (fp, "\t\t\tngs\t\tEquivalent to: -j 11 -k 9 -l 17 -e 0.02 -f\n");
//  fprintf (fp, "\t\t\tpacbio\t\tEquivalent to: -j 11 -k 7 -l 21 -e 0.20\n");
//  fprintf (fp, "\t\t\tnanopore\tEquivalent to: -j %ld -k %ld -l %ld -e %.2f\n", 11, 6, 18, 0.45f);
//  fprintf (fp, "\n");
//
//  fprintf (fp, "Algorithmic options:\n");
//  fprintf (fp, "\t-j INT\tInitial kmer size. [%ld]\n", DEFUALT_K_REGION);
//  fprintf (fp, "\t-k INT\tSecond kmer size. [%ld]\n", DEFAULT_K_GRAPH);
//  fprintf (fp, "\t-l INT\tNumber of links [%ld]\n", DEFAULT_NUM_LINKS);
//  fprintf (fp, "\t-e FLT\tApproximate error rate of the input read sequences [%.2f]\n", DEFAULT_ERROR_RATE);
//  fprintf (fp, "\t-m INT\tMaximum number of hits per kmer [%ld]\n", DEFAULT_MAX_NUM_HITS);
//  fprintf (fp, "\t-f\tThis parameter disables filtering by edit distance in the final step of mapping\n");
//  #ifndef RELEASE_VERSION
//  fprintf (fp, "\t-q INT\tThe number of selected regions if exceeded, attempt to filter them [%ld]\n", DEFAULT_MAX_NUM_REGIONS_CUTOFF);
//  fprintf (fp, "\t-g INT\tIf the final number of regions (after filtering) exceeds this amount, read will be called unmapped [%ld]\n", DEFAULT_MAX_NUM_REGIONS);
//    fprintf (fp, "\t-w INT\tKmer step during graph construction (the number of bases to skip between beginnings of every adjacent kmer) [%ld]\n", DEFAULT_KMER_GRAPH_STEP);
//    fprintf (fp, "\t-p\tOne kmer of a read can have multiple hits within the same region. This parameter enables using multiple hits per region\n");
//  #endif
//  fprintf (fp, "\n");
//
//  fprintf (fp, "Other options:\n");
//  fprintf (fp, "\t-t INT\tNumber of threads to use. Value of '-1' utilizes number of threads equal to the number of cores [%d]\n", DEFAULT_NUM_THREADS);
//  fprintf (fp, "\t-v INT\tVerbose level. If equal to 0 nothing except strict output will be placed on stdout [%ld]\n", DEFAULT_VERBOSE_LEVEL);
//  fprintf (fp, "\t-h    \tView this help.\n", DEFAULT_VERBOSE_LEVEL);
//  fprintf (fp, "\t-E    \tShow the End User Licence Agreement.\n", DEFAULT_VERBOSE_LEVEL);
//  #ifndef RELEASE_VERSION
//    fprintf (fp, "\t-s INT\tID of the read from the input reads file from which to start processing data (i.e. its ordinal number) [%ld]\n", DEFAULT_START_READ);
//    fprintf (fp, "\t-n INT\tNumber of reads to process. Value of '-1' processes all reads [%ld]\n", DEFAULT_NUM_READS_TO_PROCESS);
//  #endif
//  fprintf (fp, "\n");
//
//  #ifndef RELEASE_VERSION
//    fprintf (fp, "Debug options:\n");
//    fprintf (fp, "\t-y INT\tID of the read to give the detailed verbose output [%ld]\n", DEFAULT_DEBUG_READ);
//    fprintf (fp, "\t-b INT\tHelpful debug comments can be placed in SAM output lines (at the end). Comments can be turned of by setting this variable to 0. Different values increase/decrease verbosity level [%ld]\n", DEFAULT_VERBOSE_SAM_OUTPUT);
//    fprintf (fp, "\n");
//  #endif
//
//  fprintf (fp, "Example usage:\n");
//  fprintf (fp, "\t./graphmap -x nanopore -r escherichia_coli.fa -d reads.fastq -o alignments.sam    \t# Process nanopore reads.\n");
//  fprintf (fp, "\t./graphmap -x nanopore -r escherichia_coli.fa -D reads_folder -O alignments_folder\t# Process all reads from a given folder.\n");
//  fprintf (fp, "\n");
//
//  fprintf (fp, "Contact:\n");
//  fprintf (fp, "\tivan.sovic@irb.hr, mile.sikic@fer.hr, nagarajann@gis.a-star.edu.sg\n");
//  fprintf (fp, "\n");

  exit(0);
}

void VerboseEula(FILE *fp) {
  fprintf (fp, "%s", PROGRAM_LICENCE);
  exit(0);
}
