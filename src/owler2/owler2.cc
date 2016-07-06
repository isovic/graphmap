/*
 * owler_data.cc
 *
 *  Created on: Jul 2, 2015
 *      Author: isovic
 */

#include "owler2/owler2.h"

#include <omp.h>
#include <algorithm>
#include "libs/libdivsufsort-2.0.1-64bit/divsufsort64.h"
#include "index/index_hash.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include "index_gapped_minimizer.h"



Owler2::Owler2() {
  reference_ = NULL;
  indexes_.clear();
}

Owler2::~Owler2() {
  reference_ = NULL;
  ClearIndexes_();
}

void Owler2::ClearIndexes_() {
  for (int i=0; i<indexes_.size(); i++) {
    if (indexes_[i])
      delete indexes_[i];
    indexes_[i] = NULL;
  }
  indexes_.clear();
}

void Owler2::Run(ProgramParameters& parameters) {
  clock_t time_start = clock();
  clock_t last_time = time_start;

  if (parameters.num_threads <= 0) {
    parameters.num_threads = std::min(24, ((int) omp_get_num_procs()) / 2);
  }

  // Set the verbose level for the execution of this program.
  LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters.verbose_level);

  LOG_ALL("Loading reference sequences in memory.\n");
  SequenceFile ref_seqs(parameters.reference_path);
  ref_seqs.ConvertDataFormat(kDataFormat2BitSparse);
  // If the reads path is the same as the reference path, or reads path is not assigned, reference seqs will be used as query also.
  SequenceFile *query_seqs = &ref_seqs;
  LOG_ALL("Loaded %ld seqs.\n", ref_seqs.get_sequences().size());

  // Build the index on the fly.
  LOG_ALL("Building the index.\n");
  std::vector<std::string> index_shapes = {"1111110111111"};
  std::vector<CompiledShape> compiled_shapes = CompileShapes(index_shapes);
  float min_qv = -1.0f;
//  IndexBrute index;
  IndexGappedMinimizer index;
  index.CreateFromSequenceFile(ref_seqs, compiled_shapes, min_qv, true, true, 5, parameters.num_threads);
  LOG_ALL("Finished building the index.\n");

  SequenceFile read_seqs;
  if (parameters.reads_path != parameters.reference_path) {
    LOG_ALL("Loading query sequences in memory.\n");
    read_seqs.LoadAll(SEQ_FORMAT_AUTO, parameters.reads_path, false);
    query_seqs = &read_seqs;
    LOG_ALL("Loaded %ld seqs.\n", query_seqs->get_sequences().size());
  }

  LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
  last_time = clock();

  // Set up the debug options.
  if (parameters.debug_read_by_qname != "") {
    for (int64_t i=0; i<query_seqs->get_sequences().size(); i++) {
      if (std::string(query_seqs->get_sequences().at(i)->get_header()).compare(0, parameters.debug_read_by_qname.size(), parameters.debug_read_by_qname) == 0) {
        parameters.debug_read = i; break;
      }
    }
  }

  // Do the actual work.
  last_time = clock();
  FILE *fp_out = OpenOutFile_(parameters.out_sam_path); // Checks if the output file is specified. If it is not, then output to STDOUT.

  std::vector<CompiledShape> lookup_shapes;
  ProcessReads(parameters, &index, query_seqs, lookup_shapes, fp_out);

  LOG_ALL("All reads processed in %.2f sec (or %.2f CPU min).\n", (((float) (clock() - last_time))/CLOCKS_PER_SEC), ((((float) (clock() - last_time))/CLOCKS_PER_SEC) / 60.0f));
  if (fp_out != stdout)
    fclose(fp_out);
}

void Owler2::ProcessReads(const ProgramParameters& parameters, const IndexGappedMinimizer* index, const SequenceFile* reads, const std::vector<CompiledShape> &lookup_shapes, FILE* fp_out) {
  clock_t time_start = clock();
  clock_t last_time = time_start;
  int32_t num_threads = parameters.num_threads;

  // Set up the starting and ending read index.
  int64_t start_i = (parameters.debug_read >= 0) ? ((int64_t) parameters.debug_read) : (parameters.start_read >= 0)?((int64_t) parameters.start_read) : 0;
  int64_t end_i = (parameters.num_reads_to_process >= 0) ? (start_i + (int64_t) parameters.num_reads_to_process) : reads->get_sequences().size();

  // Initialize the counters.
  int64_t num_mapped = 0, num_unmapped = 0;

  LOG_ALL("Starting to process reads from %ld to %ld.\n", start_i, end_i);

  // Process all reads in parallel.
  #pragma omp parallel for num_threads(num_threads) shared(reads, parameters, last_time, num_mapped, num_unmapped, fp_out) schedule(dynamic, 1)
  for (int64_t i=start_i; i<end_i; i++) {
    uint32_t thread_id = omp_get_thread_num();

    if (thread_id == 0 && parameters.verbose_level > 5) {
      std::stringstream ss;
      ss << FormatString("\r[CPU time: %.2f sec, RSS: %ld MB] Read: %lu/%lu (%.2f%%) [m: %ld, u: %ld], length = %ld, qname: ",
                         (((float) (clock() - (last_time)))/CLOCKS_PER_SEC), getCurrentRSS()/(1024*1024),
                         i, reads->get_sequences().size(), ((float) i) / ((float) reads->get_sequences().size()) * 100.0f,
                         num_mapped, num_unmapped,
                         reads->get_sequences()[i]->get_data_length()) << reads->get_sequences()[i]->get_header();
      std::string string_buffer = FormatStringToLength(ss.str(), 140);
      LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, string_buffer, "ProcessReads");
    }

//    std::string mapping_out_str;
    OwlerResult owler_result;
    ProcessRead(parameters, index, reads->get_sequences()[i], lookup_shapes, &owler_result);
  }

}

FILE* Owler2::OpenOutFile_(std::string out_sam_path) {
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

bool Owler2::GetFileList_(std::string folder, std::vector<std::string> &ret_files) {
  ret_files.clear();

  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(folder.c_str())) != NULL) {
    // Get the list of file and folder names.
    while ((ent = readdir(dir)) != NULL) {
      ret_files.push_back(std::string(ent->d_name));
    }
    closedir (dir);

  } else {
    LogSystem::GetInstance().Error(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_FOLDER_NOT_FOUND, "Folder path: '%s'.", folder.c_str()));
    return false;
  }

  return true;
}


bool Owler2::StringEndsWith_(std::string const &full_string, std::string const &ending) {
  if (full_string.length() >= ending.length()) {
    return (full_string.compare(full_string.length() - ending.length(), ending.length(), ending) == 0);
  } else {
    return false;
  }
}

void Owler2::FilterFileList_(std::vector<std::string> &files, std::vector<std::string> &ret_read_files, std::vector<std::string> &ret_sam_files) {
  std::string ext_fasta = "fasta";
  std::string ext_fastq = "fastq";
  std::string ext_fa = "fa";
  std::string ext_fq = "fq";
  std::string ext_sam = "sam";
  std::string sam_file = "";

  ret_read_files.clear();
  ret_sam_files.clear();

  for (int64_t i=0; i<((int64_t) files.size()); i++) {
    if (StringEndsWith_(files.at(i), ext_fastq) || StringEndsWith_(files.at(i), ext_fasta)) {
      sam_file = files.at(i).substr(0, files.at(i).size() - 5) + ext_sam;
      ret_read_files.push_back(files.at(i));
      ret_sam_files.push_back(sam_file);
    }
    else if (StringEndsWith_(files.at(i), ext_fq) || StringEndsWith_(files.at(i), ext_fa)) {
      sam_file = files.at(i).substr(0, files.at(i).size() - 2) + ext_sam;
      ret_read_files.push_back(files.at(i));
      ret_sam_files.push_back(sam_file);
    }
  }
}
