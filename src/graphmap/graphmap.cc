/*
 * graphmap_se.cc
 *
 *  Created on: Aug 5, 2014
 *      Author: Ivan Sovic
 */

#include <omp.h>
#include <algorithm>
#include "libs/libdivsufsort-2.0.1-64bit/divsufsort64.h"
#include "graphmap/graphmap.h"
//#include "index/index_hash.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include "transcriptome.h"
#include "index/index_util.h"



GraphMap::GraphMap() : transcriptome_(nullptr) {
  indexes_.clear();
}

GraphMap::~GraphMap() {

}

void GraphMap::Run(ProgramParameters& parameters) {
  clock_t time_start = clock();
  clock_t last_time = time_start;

  // Set the verbose level for the execution of this program.
  LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters.verbose_level);

  // Check if the index exists, and build it if it doesn't.
  BuildIndex(parameters);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n\n", FormatMemoryConsumptionAsString().c_str()), "Index");
  last_time = clock();

  if (parameters.calc_only_index == true) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Finished generating index. Note: only index was generated due to selected program arguments.\n\n", FormatMemoryConsumptionAsString().c_str()), "Index");
    return;
  }

  if (parameters.threshold_hits) {
    LOG_MEDHIGH("Hits will be thresholded at the percentil value (percentil: %f%%, frequency: %.0f).\n", parameters.frequency_percentil*100.0, indexes_[0]->count_cutoff());
  } else {
    LOG_MEDHIGH("No thresholding will be applied during seed lookup.\n");
  }

  if (parameters.use_minimizers) {
    LOG_MEDHIGH("Minimizers will be used. Minimizer window length: %ld\n", parameters.minimizer_window);
  } else {
    LOG_MEDHIGH("Minimizers will not be used.\n");
  }

  // Dynamic calculation of the number of allowed regions. This should be relative to the genome size.
  // The following formula has been chosen arbitrarily.
  // The dynamic calculation can be overridden by explicitly stating the max_num_regions in the arguments passed to the binary.
  if (parameters.max_num_regions == 0) {
    if (this->indexes_[0]->get_data_length_forward() < 5000000){
      parameters.max_num_regions = 500;          // Limit the number of allowed regions, because log10 will drop rapidly after this point.
    } else {
      float M10 = 1000;                          // Baseline number of allowed regions. M10 is the number of allowed regions for 10Mbp reference size.
      float factor = log10(((float) this->indexes_[0]->get_data_length()) / 1000000.0f);     // How many powers of 10 above 1 million?
      parameters.max_num_regions = (int64_t) (M10 * factor);
    }
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Automatically setting the maximum allowed number of regions: max. %ld, attempt to reduce after %ld\n", parameters.max_num_regions, parameters.max_num_regions_cutoff), "Run");

  } else if (parameters.max_num_regions < 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("No limit to the maximum allowed number of regions will be set.\n"), "Run");
  }

  // The dynamic calculation can be overridden by explicitly stating the max_num_regions_cutoff in the arguments passed to the binary.
  if (parameters.max_num_regions_cutoff == 0) {
    parameters.max_num_regions_cutoff = (parameters.max_num_regions < 0) ? (parameters.max_num_regions) : (parameters.max_num_regions / 5);
  }

//  // Dynamic calculation of the number of allowed kmer hits for region selection.
//  // The following formula has been chosen arbitrarily.
//  // The correct value would be the one that calculates the mean (or median) of the kspectra and its standard deviation
//  // to detect outliers, but calculating the kspectra could be time and memory consuming for larger genomes. That is why
//  // we employ this simple heuristic.
//  // The dynamic calculation can be overridden by explicitly stating the max_num_hits in the arguments passed to the binary.
//  if (parameters.max_num_hits < 0) {
//    // This is how it was done previously.
////    int64_t num_kmers = (1 << (parameters.k_region * 2));
////    int64_t num_kmers_in_genome = (this->indexes_[0]->get_data_length_forward() * 2) - parameters.k_region + 1;
////    double average_num_kmers = ((double) num_kmers_in_genome) / ((double) num_kmers);
////    parameters.max_num_hits = (int64_t) ceil(average_num_kmers) * 500;
//    int64_t max_seed_count = 0;
////    ((IndexSpacedHashFast *) this->indexes_[0])->CalcPercentileHits(0.9999, &parameters.max_num_hits, &max_seed_count);
/////// TODO: Removed on 07.02.2017.    ((IndexSpacedHashFast *) this->indexes_[0])->CalcPercentileHits(0.9999, &parameters.max_num_hits, &max_seed_count);
//    LOG_ALL("Automatically setting the maximum number of seed hits to: %ld. Maximum seed occurrence in index: %ld.\n", parameters.max_num_hits, max_seed_count);
//
////    LogSystem::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true, FormatString("Automatically setting the maximum number of kmer hits: %ld\n", parameters.max_num_hits), "Run");
////    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true, FormatString("\tmax_num_hits = %ld\n", parameters.max_num_hits), "Run");
//  } else if (parameters.max_num_hits == 0) {
//    LOG_ALL("No limit to the maximum number of seed hits will be set in region selection.\n");
//  }

  if (parameters.is_reference_circular == false)
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Reference genome is assumed to be linear.\n"), "Run");
  else
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Reference genome is assumed to be circular.\n"), "Run");

  if (parameters.output_multiple_alignments == false)
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Only one alignment will be reported per mapped read.\n"), "Run");
  else
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("One or more similarly good alignments will be output per mapped read. Will be marked secondary.\n"), "Run");

  if (parameters.outfmt != "sam" &&
      parameters.outfmt != "afg" &&
      parameters.outfmt != "m5" &&
      parameters.outfmt != "mhap") {
    LogSystem::GetInstance().Error(SEVERITY_INT_WARNING, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_WRONG_FILE_TYPE, "Unknown output format specified: '%s'. Defaulting to SAM output.", parameters.outfmt.c_str()));
  }


  // Processing reads.
  // Reads can either be processed from a single file, or they can be processed from several files in a given folder.
  if (parameters.process_reads_from_folder == false) {    // This part processes a single given input file.
    last_time = clock();
    FILE *fp_out = OpenOutSAMFile_(parameters.out_sam_path); // Checks if the output SAM file is specified. If it is not, then output to STDOUT.

    // Do the actual work.
    ProcessReadsFromSingleFile(parameters, fp_out);
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("\n"), "[]");
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("All reads processed in %.2f sec (or %.2f CPU min).\n", (((float) (clock() - last_time))/CLOCKS_PER_SEC), ((((float) (clock() - last_time))/CLOCKS_PER_SEC) / 60.0f)), "ProcessReads");

    if (fp_out != stdout)
      fclose(fp_out);

  } else {    // This part processes all files in a specified folder.
    std::vector<std::string> file_list, file_list_out, read_files, sam_files;

    // The GetFileList_ functions also checks if the folder exists. If it doesn't exist, a fatal error is reported.
    if (GetFileList(parameters.reads_folder, file_list)) {
      // Sanity check for the output folder also. Function returns false if the folder does not exist.
      if (GetFileList(parameters.output_folder, file_list_out) == true) {
        FilterFileList(file_list, read_files, sam_files);

        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Loading reads from input folder. In total, %ld files need to be processed.\n", read_files.size()), "Run");

        clock_t all_reads_time = clock();

        for (int64_t i=0; i<((int64_t) read_files.size()); i++) {
          last_time = clock();
          parameters.reads_path = parameters.reads_folder + "/" + read_files.at(i);
          parameters.out_sam_path = parameters.output_folder + "/graphmap-" + sam_files.at(i);
          FILE *fp_out = OpenOutSAMFile_(parameters.out_sam_path); // Checks if the output SAM file is specified. If it is not, then output to STDOUT.

          // Do the actual work.
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Starting to process read file %ld/%ld ('%s').\n", (i + 1), read_files.size(), parameters.reads_path.c_str()), "ProcessReads");
          ProcessReadsFromSingleFile(parameters, fp_out);
          LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Finished processing read file %ld/%ld ('%s').\n\n", (i + 1), read_files.size(), parameters.reads_path.c_str()), "ProcessReads");

          if (fp_out != stdout)
            fclose(fp_out);
        }

        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("\n"), "[]");
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("All reads processed in %.2f sec (or %.2f CPU min). =====\n", (((float) (clock() - all_reads_time))/CLOCKS_PER_SEC), ((((float) (clock() - all_reads_time))/CLOCKS_PER_SEC) / 60.0f)), "ProcessReads");
      }
    }

    if (read_files.size() == 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("No read files found in path '%s'. Exiting.\n\n", parameters.reads_folder.c_str()), "Run");
    }
  }
}

bool FileExists(const std::string &path) {
  FILE *fp = fopen(path.c_str(), "r");
  if (fp == NULL) {
    return false;
  }
  fclose (fp);
  return true;
}

int GraphMap::BuildIndex(ProgramParameters &parameters) {
  indexes_.clear();
  transcriptome_ = is::createTranscriptome();

  if (parameters.is_transcriptome) {
    LOG_ALL("Loading GTF annotations.\n");
    transcriptome_->LoadGTF(parameters.gtf_path);
  }

  if (parameters.index_on_the_fly) {
    LOG_ALL("Index will be generated on the fly (won't be stored to disk). However, if it already exists it will be loaded from disk.\n");
  }

  std::string prim_index_path = parameters.index_file;
  std::string sec_index_path = parameters.index_file + std::string("sec");

  bool exists = FileExists(prim_index_path) && (!parameters.use_double_index || FileExists(sec_index_path));
  bool rebuild_index = !exists || parameters.rebuild_index || parameters.calc_only_index;

  if (rebuild_index) {
    LOG_ALL("Building the index.\n");

    // Either load genomic sequence or generate a transcriptome.
    std::shared_ptr<SequenceFile> refs = nullptr;
    if (!parameters.is_transcriptome) {
      LOG_ALL("Loading reference sequences.\n");
      refs = std::shared_ptr<SequenceFile>(new SequenceFile(parameters.reference_path));
    } else {
      LOG_ALL("Loading genomic sequences.\n");
      auto genomic = std::shared_ptr<SequenceFile>(new SequenceFile(parameters.reference_path));
      LOG_ALL("Generating the transcriptome.\n");
      refs = transcriptome_->GenerateTranscriptomeSeqs(genomic);
    }

    LOG_ALL("Constructing the primary index.\n");
    std::vector<std::string> shapes_prim = {"11110111101111"};
    auto index_prim = is::createMinimizerIndex(shapes_prim, parameters.frequency_percentil);

    index_prim->Create(*refs, 0.0f, true, parameters.use_minimizers, parameters.minimizer_window, parameters.num_threads, true);

    LOG_ALL("Storing the primary index to file: '%s'.\n", prim_index_path.c_str());
    if (parameters.index_on_the_fly == false) {
      index_prim->Store(prim_index_path);
    }
    indexes_.push_back(index_prim);

    if (parameters.use_double_index) {
      LOG_ALL("Sensitive mode selected.\n");
      LOG_ALL("Constructing the secondary index.\n");
      std::vector<std::string> shapes_sec = {"1111110111111"};
      auto index_sec = is::createMinimizerIndex(shapes_sec, parameters.frequency_percentil);

      index_sec->Create(*refs, 0.0f, true, parameters.use_minimizers, parameters.minimizer_window, parameters.num_threads, true);

      LOG_ALL("Storing the secondary index to file: '%s'.\n", sec_index_path.c_str());
      if (parameters.index_on_the_fly == false) {
        index_sec->Store(sec_index_path);
      }
      indexes_.push_back(index_sec);
    }

  } else {      // Load the index.
    LOG_ALL("Loading the primary index.\n");
    auto index_prim = is::createMinimizerIndex(prim_index_path, parameters.frequency_percentil);
    indexes_.push_back(index_prim);

    if (parameters.use_double_index) {
      LOG_ALL("Loading the secondary index.\n");
      auto index_sec = is::createMinimizerIndex(sec_index_path, parameters.frequency_percentil);
      indexes_.push_back(index_sec);
    }
  }

  return 0;
}

//void DeprecatedBuildIndex(ProgramParameters &parameters) {
//  // Run away, you are free now!
//  for (int32_t i=0; i<indexes_.size(); i++) {
//    if (indexes_[i]) { delete indexes_[i]; }
//    indexes_[i] = NULL;
//  }
//  indexes_.clear();
//
//  IndexSpacedHashFast *index_prim = new IndexSpacedHashFast(SHAPE_TYPE_444);
//  IndexSpacedHashFast *index_sec = NULL;
//  indexes_.push_back(index_prim);
//
//  if (parameters.sensitive_mode == false) {
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Running in normal (parsimonious) mode. Only one index will be used.\n"), "Index");
//  } else {
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Running in sensitive mode. Two indexes will be used (double memory consumption).\n"), "Index");
//    index_sec = new IndexSpacedHashFast(SHAPE_TYPE_66);
//    indexes_.push_back(index_sec);
//  }
//
//  clock_t last_time = clock();
//
//  if (parameters.calc_only_index == false) {
//    // Check if index already exists, if not generate it.
//    FILE *fp = fopen(parameters.index_file.c_str(), "r");
//    if (fp == NULL) {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Index is not prebuilt. Generating index.\n"), "Index");
//    } else {
//      fclose (fp);
//      if (parameters.rebuild_index == false) {
//        LOG_ALL("Index already exists. Loading from file.\n");
//      } else {
//        LOG_ALL("Index already exists, but will be rebuilt.\n");
//      }
//    }
//
//    // Check whether the index needs to be rebuilt, or if it can only be loaded.
//    transcriptome_ = is::createTranscriptome();
//    if (parameters.rebuild_index == false) {
//      int prim_index_loaded = 0;
//      if (parameters.gtf_path == "") {
//        index_prim->LoadOrGenerate(parameters.reference_path, parameters.index_file, (parameters.verbose_level > 0));
//      } else {
//        index_prim->LoadOrGenerateTranscriptome(parameters.reference_path, parameters.gtf_path, parameters.index_file, (parameters.verbose_level > 0));
//      }
//
//      if (prim_index_loaded) { return 1; }
//    } else {
//      int prim_index_generated = 0;
//      if (parameters.gtf_path == "") {
//        index_prim->GenerateFromFile(parameters.reference_path);
//      } else {
//        index_prim->GenerateTranscriptomeFromFile(parameters.reference_path, parameters.gtf_path);
//      }
//
//      int prim_index_stored = index_prim->StoreToFile(parameters.index_file);
//      if (prim_index_generated || prim_index_stored) { return 1; }
//    }
//
//    if (parameters.sensitive_mode == true ) {
//      fp = fopen((parameters.index_file + std::string("sec")).c_str(), "r");
//      if (fp == NULL) {
//        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Secondary index is not prebuilt. Generating index.\n"), "Index");
//      } else {
//        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Secondary index already exists. Loading from file.\n"), "Index");
//        fclose (fp);
//      }
//
//      if (parameters.rebuild_index == false) {
//        int sec_index_loaded = 0;
//        if (parameters.gtf_path == "") {
//          index_sec->LoadOrGenerate(parameters.reference_path, parameters.index_file + std::string("sec"), (parameters.verbose_level > 0));
//        } else {
//          index_sec->LoadOrGenerateTranscriptome(parameters.reference_path, parameters.gtf_path, parameters.index_file + std::string("sec"), (parameters.verbose_level > 0));
//        }
//        if (sec_index_loaded) { return 1; }
//      } else {
//        int sec_index_generated = 0;
//        if (parameters.gtf_path == "") {
//          index_sec->GenerateFromFile(parameters.reference_path);
//        } else {
//          index_sec->GenerateTranscriptomeFromFile(parameters.reference_path, parameters.gtf_path);
//        }
//        int sec_index_stored = index_sec->StoreToFile(parameters.index_file + std::string("sec"));
//        if (sec_index_generated || sec_index_stored) { return 1; }
//      }
//    }
//
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Index loaded in %.2f sec.\n", (((float) (clock() - last_time))/CLOCKS_PER_SEC)), "Index");
//    return 0;
//
//  } else {
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Generating index.\n"), "Index");
//
//    if (parameters.gtf_path == "") {
//      index_prim->GenerateFromFile(parameters.reference_path);
//    } else {
//      index_prim->GenerateTranscriptomeFromFile(parameters.reference_path, parameters.gtf_path);
//    }
//
//    index_prim->StoreToFile(parameters.index_file);
//
//    if (parameters.sensitive_mode == true) {
//      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Generating secondary index.\n"), "Index");
//      index_sec->GenerateFromFile(parameters.reference_path);
//      index_sec->StoreToFile(parameters.index_file + std::string("sec"));
//    }
//    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Index generated in %.2f sec.\n", (((float) (clock() - last_time))/CLOCKS_PER_SEC)), "Index");
//  }
//
//  return 0;
//}

void GraphMap::ProcessReadsFromSingleFile(ProgramParameters &parameters, FILE *fp_out) {
  // Write out the SAM header in fp_out.
  if (parameters.outfmt == "sam") {
    std::string sam_header;
    if (parameters.gtf_path.size() == 0) {
      sam_header = is::GenerateSAMHeader(indexes_[0], parameters);
    } else {
      sam_header = is::GenerateSAMHeader(transcriptome_);
    }

    if (sam_header.size() > 0)
      fprintf (fp_out, "%s\n", sam_header.c_str());
  }

  // Check whether to load in batches or to load all the data at once.
  if (parameters.batch_size_in_mb <= 0) {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("All reads will be loaded in memory.\n"), "ProcessReads");
  } else {
    LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Reads will be loaded in batches of up to %ld MB in size.\n", parameters.batch_size_in_mb), "ProcessReads");
  }

  SequenceFile reads;
  reads.OpenFileForBatchLoading(parameters.reads_path);

  clock_t absolute_time = clock();
  clock_t last_batch_loading_time = clock();

  int64_t num_mapped = 0;
  int64_t num_unmapped = 0;

  // Load sequences in batch (if requested), or all at once.
  while ((parameters.batch_size_in_mb <= 0 && !reads.LoadAllAsBatch(SeqFmtToString(parameters.infmt), false)) || (parameters.batch_size_in_mb > 0 && !reads.LoadNextBatchInMegabytes(SeqFmtToString(parameters.infmt), parameters.batch_size_in_mb, false))) {
    if (parameters.batch_size_in_mb <= 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("All reads loaded in %.2f sec (size around %ld MB). (%ld bases)\n", (((float) (clock() - last_batch_loading_time))/CLOCKS_PER_SEC), reads.CalculateTotalSize(MEMORY_UNIT_MEGABYTE), reads.GetNumberOfBases()), "ProcessReads");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), "ProcessReads");
    }
    else {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Batch of %ld reads (%ld MiB) loaded in %.2f sec. (%ld bases)\n", reads.get_sequences().size(), reads.CalculateTotalSize(MEMORY_UNIT_MEGABYTE), parameters.reads_path.c_str(), (((float) (clock() - last_batch_loading_time))/CLOCKS_PER_SEC), reads.GetNumberOfBases()), "ProcessReads");
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), "ProcessReads");
    }

    // This line actually does all the work.
    ProcessSequenceFileInParallel(&parameters, &reads, &absolute_time, fp_out, &num_mapped, &num_unmapped);

    if (parameters.batch_size_in_mb > 0) {
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("\n"), "[]");
    }

    last_batch_loading_time = clock();
  }

  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), "ProcessReads");

  reads.CloseFileAfterBatchLoading();
}

int GraphMap::ProcessSequenceFileInParallel(ProgramParameters *parameters, SequenceFile *reads, clock_t *last_time, FILE *fp_out, int64_t *ret_num_mapped, int64_t *ret_num_unmapped) {
  int64_t num_reads = reads->get_sequences().size();
  std::vector<std::string> sam_lines;

  if (parameters->output_in_original_order == true) {
    sam_lines.resize(num_reads, std::string(""));
  }

  // Division by to to avoid hyperthreading cores, and limit on 24 to avoid clogging a shared SMP.
//  int64_t num_threads = std::min(24, ((int) omp_get_num_procs()) / 2);
//
//  if (parameters->num_threads > 0)
  int64_t num_threads = (int64_t) parameters->num_threads;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Using %ld threads.\n", num_threads), "ProcessReads");

  // Set up the starting and ending read index.
  int64_t start_i = (parameters->start_read >= 0)?((int64_t) parameters->start_read):0;

  #ifndef RELEASE_VERSION
    if (parameters->debug_read >= 0)
      start_i = parameters->debug_read;

    if (parameters->debug_read_by_qname != "") {
      for (int64_t i=0; i<num_reads; i++) {
        if (std::string(reads->get_sequences().at(i)->get_header()).compare(0, parameters->debug_read_by_qname.size(), parameters->debug_read_by_qname) == 0) {
          start_i = i;
          parameters->debug_read = i;
          break;
        }
      }
    }
  #endif

  int64_t max_i = (parameters->num_reads_to_process >= 0) ? (start_i + (int64_t) parameters->num_reads_to_process) : num_reads;

  // Initialize the counters.
  int64_t num_mapped=0, num_unmapped=0, num_ambiguous=0, num_errors=0;
  int64_t num_reads_processed_in_thread_0 = 0;

  EValueParams *evalue_params;
  SetupScorer((char *) "EDNA_FULL_5_4", indexes_[0]->get_data_length_forward(), -parameters->evalue_gap_open, -parameters->evalue_gap_extend, &evalue_params);

  // Process all reads in parallel.
  #pragma omp parallel for num_threads(num_threads) firstprivate(num_reads_processed_in_thread_0, evalue_params) shared(reads, parameters, last_time, sam_lines, num_mapped, num_unmapped, num_ambiguous, num_errors, fp_out) schedule(dynamic, 1)
  for (int64_t i=start_i; i<max_i; i++) {
    uint32_t thread_id = omp_get_thread_num();

    // Verbose the currently processed read. If the verbose frequency is low, only output to STDOUT every 100th read.
    // If medium verbose frequency is set, every 10th read will be output, while for high every read will be reported.
    if (thread_id == 0 && parameters->verbose_level > 0) {
      if (((!(LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_ALL) ||
            (LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_LOW)) && (num_reads_processed_in_thread_0 % 100) == 0) ||
          ((LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_MED) && (num_reads_processed_in_thread_0 % 10) == 0) ||
          ((LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_HIGH))) {

        std::stringstream ss;
        ss << FormatString("\r[CPU time: %.2f sec, RSS: %ld MB] Read: %lu/%lu (%.2f%%) [m: %ld, u: %ld], length = %ld, qname: ",
                           (((float) (clock() - (*last_time)))/CLOCKS_PER_SEC), getCurrentRSS()/(1024*1024),
                           i, reads->get_sequences().size(), ((float) i) / ((float) reads->get_sequences().size()) * 100.0f,
                           num_mapped, num_unmapped,
                           reads->get_sequences()[i]->get_data_length()) << reads->get_sequences()[i]->get_header();
        std::string string_buffer = FormatStringToLength(ss.str(), 140);
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, string_buffer, "ProcessReads");

        if (parameters->verbose_level > 6 && parameters->num_threads == 1)
              ss << "\n";
      }

      #pragma omp critical
      {
        num_reads_processed_in_thread_0 += 1;
      }
    }

    // The actual interesting part.
    std::string sam_line = "";
    MappingData mapping_data;
    ProcessRead(&mapping_data, indexes_, transcriptome_, reads->get_sequences()[i], parameters, evalue_params);
//    ProcessRead2(&mapping_data, indexes_, transcriptome_, reads->get_sequences()[i], parameters, evalue_params);

    // Generate the output.
    int mapped_state = STATE_UNMAPPED;
    mapped_state = CollectAlignments(reads->get_sequences()[i], parameters, &mapping_data, sam_line);

    // Keep the counts.
    if (mapped_state == STATE_MAPPED) {
      #pragma omp critical
      num_mapped += 1;
    }
    else if (mapped_state == STATE_UNMAPPED) {
      #pragma omp critical
      num_unmapped += 1;
    }
    else if (mapped_state == STATE_AMBIGUOUS) {
      #pragma omp critical
      num_ambiguous += 1;
    }
    else {
      #pragma omp critical
      num_errors += 1;
    }

    // If the order of the reads should be kept, store them in a vector, otherwise output the alignment to file.
    if (parameters->output_in_original_order == false) {
      if (sam_line.size() > 0) {
        #pragma omp critical
        fprintf (fp_out, "%s\n", sam_line.c_str());
      }
    }
    else {
      #pragma omp critical
      sam_lines[i] = sam_line;
    }
  }

  (*ret_num_mapped) = num_mapped;
  (*ret_num_unmapped) = num_unmapped;

  if (evalue_params) {
    DeleteEValueParams(evalue_params);
  }

  // Verbose the final processing info.
  std::string string_buffer = FormatString("\r[CPU time: %.2f sec, RSS: %ld MB] Read: %lu/%lu (%.2f%%) [m: %ld, u: %ld]",
                               (((float) (clock() - (*last_time)))/CLOCKS_PER_SEC), getCurrentRSS()/(1024*1024),
                               reads->get_sequences().size(), reads->get_sequences().size(), 100.0f,
                               num_mapped, num_unmapped);
  string_buffer = FormatStringToLength(string_buffer, 140);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, string_buffer, "ProcessReads");
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, "\n", "[]");

  // Output the results to the SAM file in the exact ordering of the input file (if it was requested by the specified parameter).
  if (parameters->output_in_original_order == true) {
    for (int64_t i=0; i<num_reads; i++) {
      if (sam_lines[i].size() > 0) {
        fprintf (fp_out, "%s\n", sam_lines[i].c_str());
      }
    }
  }

  return 0;
}

FILE* GraphMap::OpenOutSAMFile_(std::string out_sam_path) {
  // Check if the output SAM file is specified. If it is not, then output to STDOUT.
  FILE *fp_out = stdout;

  if (out_sam_path.size() > 0) {
    fp_out = fopen(out_sam_path.c_str(), "w");

    if (fp_out == NULL) {
      LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_OPENING_FILE, "File path: '%s'.", out_sam_path.c_str()));
      return NULL;
    }
  }

  return fp_out;
}
