/*
 * owler_data.cc
 *
 *  Created on: Jul 2, 2015
 *      Author: isovic
 */

#include "owler/owler.h"

#include <omp.h>
#include <algorithm>
#include "log_system/log_system.h"
#include "utility/utility_general.h"

#include "utility/tictoc.h"



Owler::Owler() {

}

Owler::~Owler() {
}

void Owler::Run(ProgramParameters& parameters) {
  // Set the verbose level for the execution of this program.
  LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters.verbose_level);

  TicToc tt_all;
  tt_all.start();

  LOG_ALL("Loading genomic sequences.\n");
  TicToc tt_load;
  tt_load.start();
  ref_ = std::shared_ptr<SequenceFile>(new SequenceFile(parameters.reference_path));
  if (parameters.reads_path == parameters.reference_path) {
    reads_ = ref_;
  } else {
    reads_ = std::shared_ptr<SequenceFile>(new SequenceFile(parameters.reads_path));
  }
  tt_load.stop();
  LOG_ALL("All sequences loaded in %.2f sec (size of reads file around %ld MB). (%ld bases)\n", tt_load.get_secs(), reads_->CalculateTotalSize(MEMORY_UNIT_MEGABYTE), reads_->GetNumberOfBases());
  LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());

  // Construct the index.
  TicToc tt_index;
  tt_index.start();
  BuildIndex_(parameters);
  tt_index.stop();
  LOG_MEDHIGH("Memory consumption: %s\n\n", FormatMemoryConsumptionAsString().c_str());

  // Processing reads.
  TicToc tt_processing;
  tt_processing.start();

  FILE *fp_out = OpenOutFile_(parameters.out_sam_path); // Checks if the output file is specified. If it is not, then output to STDOUT.

  // Do the actual work.
  ProcessSequenceFileInParallel_(parameters, reads_, tt_all, fp_out);

  if (fp_out != stdout) {
    fclose(fp_out);
  }

  tt_processing.stop();

  LOG_NEWLINE;
  LOG_ALL("All reads processed in %.2f sec (or %.2f CPU min).\n", tt_processing.get_secs(), tt_processing.get_secs() / 60.0f);


}

//int Owler::BuildIndex_(ProgramParameters &parameters) {
//  LOG_ALL("Building the index.\n");
//  // Division by 2 to to avoid hyperthreading cores, and limit
//  // to 24 to avoid clogging a shared SMP.
//  int64_t num_threads = (parameters.num_threads > 0) ?
//                              parameters.num_threads :
//                              std::min(24, ((int) omp_get_num_procs()) / 2);
//
//  std::vector<std::string> shapes_prim = {"1111110111111"};
//  auto index_prim = is::createMinimizerIndex(shapes_prim, parameters.frequency_percentil);
//  index_prim->Create(*ref_, 0.0f, true, parameters.use_minimizers, parameters.minimizer_window, num_threads, true);
//
//  if (index_prim == nullptr) {
//    FATAL_REPORT(ERR_UNEXPECTED_VALUE, "No index was generated! Exiting.");
//  }
//
//  index_ = index_prim;
//
//  return 0;
//}

int Owler::BuildIndex_(ProgramParameters &parameters) {
  // Division by 2 to to avoid hyperthreading cores, and limit
  // to 24 to avoid clogging a shared SMP.
  int64_t num_threads = (parameters.num_threads > 0) ?
                              parameters.num_threads :
                              std::min(24, ((int) omp_get_num_procs()) / 2);

  std::vector<std::string> shapes_prim = {"1111110111111"};
  index_ = is::createMinimizerIndex(shapes_prim, parameters.frequency_percentil);

  std::string index_path = parameters.index_file + "owl";

  bool load = !parameters.index_on_the_fly && (!parameters.rebuild_index && FileExists(index_path));

  if (load) {
    LOG_ALL("Loading the index from file.\n");

    index_->Load(index_path);

  } else {
    LOG_ALL("Building the index.\n");

    index_->Create(*ref_, 0.0f, true, parameters.use_minimizers, parameters.minimizer_window, num_threads, true);

    if (!parameters.index_on_the_fly) {
      index_->Store(index_path);
    }
  }

  if (index_ == nullptr) {
    FATAL_REPORT(ERR_UNEXPECTED_VALUE, "No index was generated! Exiting.");
  }

  return 0;
}

int Owler::ProcessSequenceFileInParallel_(ProgramParameters &parameters, std::shared_ptr<SequenceFile> reads, TicToc &tt_all, FILE *fp_out) {
  if (parameters.outfmt == "dot") {
    fprintf (fp_out, "digraph overlaps {\n");
    fprintf (fp_out, "\tnode [shape=circle,fixedsize=true,style=filled,fillcolor=grey];\n");
  }

  int64_t num_reads = reads->get_sequences().size();
  // Division by 2 to to avoid hyperthreading cores, and limit
  // to 24 to avoid clogging a shared SMP.
  int64_t num_threads = (parameters.num_threads > 0) ?
                              parameters.num_threads :
                              std::min(24, ((int) omp_get_num_procs()) / 2);

  LOG_MEDHIGH("Using %ld threads.\n", num_threads);

  // Set up the starting and ending read index.
  int64_t start_i = (parameters.start_read >= 0)?((int64_t) parameters.start_read):0;

  #ifndef RELEASE_VERSION
    if (parameters.debug_read >= 0)
      start_i = parameters.debug_read;

    if (parameters.debug_read_by_qname != "") {
      for (int64_t i=0; i<num_reads; i++) {
        if (std::string(reads->get_sequences().at(i)->get_header()).compare(0, parameters.debug_read_by_qname.size(), parameters.debug_read_by_qname) == 0) {
          start_i = i;
          parameters.debug_read = i;
          break;
        }
      }
    }
  #endif

  int64_t max_i = (parameters.num_reads_to_process >= 0) ? (start_i + (int64_t) parameters.num_reads_to_process) : num_reads;

  // Initialize the counters.
  int64_t num_mapped=0, num_unmapped=0, num_ambiguous=0, num_errors=0;
  int64_t num_reads_processed_in_thread_0 = 0;

  // Process all reads in parallel.
  #pragma omp parallel for num_threads(num_threads) firstprivate(num_reads_processed_in_thread_0) shared(reads, parameters, num_mapped, num_unmapped, num_ambiguous, num_errors, fp_out) schedule(dynamic, 1)
  for (int64_t i=start_i; i<max_i; i++) {
    uint32_t thread_id = omp_get_thread_num();

    // Verbose the currently processed read. If the verbose frequency is low, only output to STDOUT every 100th read.
    // If medium verbose frequency is set, every 10th read will be output, while for high every read will be reported.
    if (thread_id == 0 && parameters.verbose_level > 0) {
      if (((!(LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_ALL) ||
            (LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_LOW)) && (num_reads_processed_in_thread_0 % 100) == 0) ||
          ((LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_MED) && (num_reads_processed_in_thread_0 % 10) == 0) ||
          ((LogSystem::GetInstance().PROGRAM_VERBOSE_LEVEL & VERBOSE_FREQ_HIGH))) {

        std::stringstream ss;
        if (parameters.verbose_level > 6 && parameters.num_threads == 1)
              ss << "\n";
        ss << FormatString("\r[CPU time: %.2f sec, RSS: %ld MB] Read: %lu/%lu (%.2f%%) [m: %ld, u: %ld], length = %ld, qname: ",
                           tt_all.get_secs_current(), getCurrentRSS()/(1024*1024),
                           i, reads->get_sequences().size(), ((float) i) / ((float) reads->get_sequences().size()) * 100.0f,
                           num_mapped, num_unmapped,
                           reads->get_sequences()[i]->get_data_length()) << reads->get_sequences()[i]->get_header();
        std::string string_buffer = FormatStringToLength(ss.str(), 140);
        LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, string_buffer, "ProcessReads");

        if (parameters.verbose_level > 6 && parameters.num_threads == 1)
              ss << "\n";
      }

      #pragma omp critical
      {
        num_reads_processed_in_thread_0 += 1;
      }
    }

    // The actual interesting part.

    OwlerData owler_data;
    ProcessRead_(index_, reads->get_sequences()[i], &parameters, owler_data);

    int mapped_state = STATE_UNMAPPED;
    if (owler_data.overlaps.size() > 0) {
      mapped_state = STATE_MAPPED;

      #pragma omp critical
      {
        fprintf (fp_out, "%s", owler_data.overlap_lines.c_str());
      }
    }

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
  }

  // Verbose the final processing info.
  std::string string_buffer = FormatString("\r[CPU time: %.2f sec, RSS: %ld MB] Read: %lu/%lu (%.2f%%) [m: %ld, u: %ld]",
                               tt_all.get_secs_current(), getCurrentRSS()/(1024*1024),
                               reads->get_sequences().size(), reads->get_sequences().size(), 100.0f,
                               num_mapped, num_unmapped);
  string_buffer = FormatStringToLength(string_buffer, 140);
  LOG_ALL("%s\n", string_buffer.c_str());

  if (parameters.outfmt == "dot") {
    fprintf (fp_out, "}\n");
  }

  return 0;
}

FILE* Owler::OpenOutFile_(std::string out_path) {
  // Check if the output SAM file is specified. If it is not, then output to STDOUT.
  FILE *fp_out = stdout;

  if (out_path.size() > 0) {
    fp_out = fopen(out_path.c_str(), "w");
    if (fp_out == NULL) {
      FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", out_path.c_str());
      return NULL;
    }
  }

  return fp_out;
}

