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
#include <iostream>

#include "aligner/aligner_containers.h"
#include "aligner/aligner_ksw2.h"
#include "aligner/anchor_aligner.h"
#include "aligner/aligner_util.hpp"
#include "aligner/pairwise_penalties.h"

GraphMap::GraphMap() : transcriptome_(nullptr) {
  indexes_.clear();
}

GraphMap::~GraphMap() {

}

void GraphMap::Run(ProgramParameters& parameters) {
  clock_t time_start = clock();
  clock_t last_time = time_start;

  if (parameters.composite_parameters == "rnaseq") {
	  parameters.output_in_original_order = true;
	  parameters.batch_size_in_mb = -1;
  }

  // Set the verbose level for the execution of this program.
  LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters.verbose_level);

  // Check if the index exists, and build it if it doesn't.
  BuildIndexes(parameters);
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str()), "Index");
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

std::shared_ptr<is::MinimizerIndex> GraphMap::SetupIndex_(std::shared_ptr<SequenceFile> ref, const std::string &index_path, const std::string &shape,
                          const ProgramParameters &parameters, int64_t num_threads) const {

  std::vector<std::string> shapes = {shape};
  auto index = is::createMinimizerIndex(shapes, parameters.frequency_percentil);

  bool load = !parameters.rebuild_index && FileExists(index_path);
  bool store = !parameters.index_on_the_fly;

  if (load) {
    LOG_ALL("Loading index from file: '%s'.\n", index_path.c_str());

    int ret_load = index->Load(index_path);

    if (ret_load) {
      LOG_DEBUG("Problems loading index: index is of wrong version or index file corrupt.\n");

      if (parameters.auto_rebuild_index) {
        LOG_ALL("Rebuilding the index.\n");
        load = false;
      } else {
        FATAL_REPORT(ERR_GENERIC, "Not rebuilding the index automatically (specify --auto-rebuild-index).");
      }
    }
  }

  // Separate 'if' because there can be a fallthrough case when index needs to be rebuilt.
  if (!load) {
    LOG_ALL("Building the index for shape: '%s'.\n", shape.c_str());

    index->Create(*ref, 0.0f, true, parameters.use_minimizers, parameters.minimizer_window, num_threads, true);
    LOG_ALL("Finished building index.\n");

    if (store) {
      LOG_ALL("Storing the index to file: '%s'.\n", index_path.c_str());
      index->Store(index_path);
    }
  }

  if (index == nullptr) {
    FATAL_REPORT(ERR_GENERIC, "No index was generated! Exiting.");
  }

  return index;
}

int GraphMap::BuildIndexes(ProgramParameters &parameters) {

  // Division by 2 to to avoid hyperthreading cores, and limit
  // to 24 to avoid clogging a shared SMP.
  int64_t num_threads = (parameters.num_threads > 0) ?
                              parameters.num_threads :
                              std::min(24, ((int) omp_get_num_procs()) / 2);

  indexes_.clear();
  transcriptome_ = is::createTranscriptome();

  if (parameters.index_on_the_fly) {
    LOG_ALL("Index will be generated on the fly (won't be stored to disk). However, if it already exists it will be loaded from disk.\n");
  }

  if (parameters.is_transcriptome) {
    LOG_ALL("Loading GTF annotations.\n");
    transcriptome_->LoadGTF(parameters.gtf_path);
  }

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

  // Construct the primary index.
  std::string prim_index_path = parameters.index_file;
  auto index_prim = SetupIndex_(refs, prim_index_path, "11110111101111", parameters, num_threads);
  indexes_.push_back(index_prim);

  // Construct the secondary index.
  if (parameters.use_double_index) {
    LOG_ALL("Sensitive mode selected.\n");
    LOG_ALL("Constructing the secondary index.\n");
    std::string sec_index_path = parameters.index_file + std::string("sec");
    auto index_sec = SetupIndex_(refs, sec_index_path, "1111110111111", parameters, num_threads);
    indexes_.push_back(index_sec);
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

bool GraphMap::comparePtrToNode(RealignmentStructure* a, RealignmentStructure* b) { return ((*a).start < (*b).start); }

int get_read_for_reference_begining(std::vector<is::CigarOp> cigar, int len) {
	int reference_covered = 0;
	int read_covered = 0;

	for(is::CigarOp &c: cigar) {
		int count = c.count;

		if(c.op == '=' || c.op == 'X') {
			if(reference_covered + count > len) {
				read_covered += len - reference_covered;
				break;
			} else if(reference_covered + count == len) {
				read_covered += count;
				break;
			} else {
				read_covered += count;
				reference_covered += count;
			}
		}

		if(c.op == 'D') {
			if(reference_covered + count >= len) {
				break;
			} else {
				reference_covered += count;
			}
		}

		if(c.op == 'I') {
			read_covered += count;
		}
	}

	return read_covered;
}

std::vector<is::CigarOp> GraphMap::ProcessReadExons(std::vector<ExonInfo> &exonsInfos, const char *ref_data) {
	if (exonsInfos.size() <= 0) {
		std::vector<is::CigarOp> temp;
//		std::cout << "bad exit" << std::endl;
		return temp;
	}

	std::vector<is::CigarOp> complete_cigar1;

//	std::cout << exonsInfos[0].read_id << std::endl;

	std::sort(std::begin(exonsInfos), std::end(exonsInfos), [](ExonInfo a, ExonInfo b) {return a.order_number < b.order_number; });

	std::vector<ExonInfo> merged_exons;

	ExonInfo previous_exon = exonsInfos[0];
	int index = 1;
	while(index != exonsInfos.size()) {
		ExonInfo current_exon = exonsInfos[index];
		if (abs((previous_exon.start-previous_exon.leftOffset)-(current_exon.start-current_exon.leftOffset)) <  5 ||
				abs((previous_exon.stop-previous_exon.rightOffset)-(current_exon.stop-current_exon.rightOffset)) <  5) {
			ExonInfo new_exon = ExonInfo(previous_exon, current_exon);
			new_exon.reference = previous_exon.reference;
			for (long var = previous_exon.stop; var < current_exon.start; ++var) {
				new_exon.reference.push_back(ref_data[var]);
			}
			new_exon.reference += current_exon.reference;
			previous_exon = new_exon;
		} else {
			merged_exons.push_back(previous_exon);
			previous_exon = current_exon;
		}
		index += 1;
	}
	merged_exons.push_back(previous_exon);

	previous_exon = merged_exons[0];
	index = 1;

	int base_offset = 0;

	int front_clipping = 0;
	if (previous_exon.cigar.size() > 0 && previous_exon.cigar[0].op == 'S') {
		front_clipping = previous_exon.cigar[0].count;
	}

	int first_offset_left = base_offset + std::max(-previous_exon.rightOffset, 0);
	int first_offset_right = base_offset + std::max(previous_exon.rightOffset, 0);

	int first_ref_len = ((previous_exon.reference.size() - first_offset_left) + first_offset_right);

	std::string ref;

	for (int var = previous_exon.start; var < previous_exon.stop - first_offset_left; ++var) {
		ref.push_back(ref_data[var]);
	}

	if(first_offset_right > 0) {
		for (int var = previous_exon.stop - first_offset_left; var < previous_exon.stop + first_offset_right; ++var) {
			ref.push_back(ref_data[var]);
		}
	}

	std::string read = previous_exon.content.substr(front_clipping, previous_exon.content.size()-front_clipping);

	bool shouldAling = previous_exon.rightOffset != 0;

	is::PiecewisePenalties p(2, -4, std::vector<is::AffinePiece>{is::AffinePiece(-2, -4), is::AffinePiece(-1, -13)});
	is::AlignmentOptions aln_opt;

	auto aligner = is::createAlignerKSW2(p, aln_opt);

	std::vector<is::CigarOp> complete_cigar;

	bool did_adjust_read = false;

	std::vector<is::CigarOp> container_vector;

//	for(ExonInfo ei: merged_exons) {
//		std::cout << ei.start << "-" << ei.stop << "  " << ei.leftOffset << " " << ei.rightOffset << " (" << ei.invalid << ")" << std::endl;
//	}

	while(index != merged_exons.size()) {
		ExonInfo current_exon = merged_exons[index];

		if (current_exon.invalid) {
			shouldAling = true;
			read += current_exon.content;
		} else {
			if (current_exon.leftOffset != 0 || shouldAling) {

				int second_offset_left = base_offset + std::max(current_exon.leftOffset, 0);
				int second_offset_right = base_offset + std::max(-current_exon.leftOffset, 0);

				std::vector<is::CigarOp> first_vector;
				std::vector<is::CigarOp> second_vector;

				int back_clipping = 0;
				if (current_exon.cigar.size() > 0 && current_exon.cigar[current_exon.cigar.size()-1].op == 'S') {
					back_clipping = current_exon.cigar[current_exon.cigar.size()-1].count;
				}

				for (int var = current_exon.start - second_offset_left; var < current_exon.start; ++var) {
					ref.push_back(ref_data[var]);
				}

				for (int var = current_exon.start + second_offset_right; var < current_exon.stop; ++var) {
					ref.push_back(ref_data[var]);
				}

				std::string total_read = read + current_exon.content.substr(0, current_exon.content.size()-back_clipping);

				did_adjust_read = true;

				aligner->Global(total_read.c_str(), total_read.size(), ref.c_str(), ref.size(), true);

				auto aln_result = aligner->getResults();

				int reference_covered = 0;
				int read_covered = 0;
				bool found_first_vector = false;

				for (auto& c: aln_result->cigar) {
					int count = c.count;

					if (found_first_vector) {
						second_vector.push_back(c);
						continue;
					}

					if(c.op == '=' || c.op == 'X') {
						if (	reference_covered + count == first_ref_len) {
							is::CigarOp left_op = is::CigarOp(c.op, first_ref_len-reference_covered);
							first_vector.push_back(left_op);
							reference_covered += count;
							read_covered += count;
							found_first_vector = true;

//							std::cout << "Ags: " << ref[reference_covered-4] << ref[reference_covered-3] << ref[reference_covered-2] << ref[reference_covered-1] << std::endl;
//							std::cout << ref[reference_covered] << ref[reference_covered+1] << ref[reference_covered+2] << ref[reference_covered+3] << ref[reference_covered+4] << ref[reference_covered+5] << ref[reference_covered+6] << ref[reference_covered+7] <<  ref[reference_covered+8] <<  std::endl;
						}
						else if(reference_covered + count > first_ref_len) {
							is::CigarOp left_op = is::CigarOp(c.op, first_ref_len-reference_covered);
							first_vector.push_back(left_op);
							is::CigarOp right_op = is::CigarOp(c.op, (reference_covered + count)-first_ref_len);
							second_vector.push_back(right_op);

							read_covered += (first_ref_len-reference_covered);
							reference_covered += (first_ref_len-reference_covered);
							found_first_vector = true;
//							std::cout << "Ags: " << ref[reference_covered-4] << ref[reference_covered-3] << ref[reference_covered-2] << ref[reference_covered-1] << std::endl;
//							std::cout << ref[reference_covered] << ref[reference_covered+1] << ref[reference_covered+2] << ref[reference_covered+3] << ref[reference_covered+4] << ref[reference_covered+5] << ref[reference_covered+6] << ref[reference_covered+7] <<  ref[reference_covered+8] <<  std::endl;
						} else {
							first_vector.push_back(c);
							reference_covered += count;
							read_covered += count;
						}
					}

					if(c.op == 'D') {
						if (	reference_covered + count == first_ref_len) {
							is::CigarOp left_op = is::CigarOp(c.op, first_ref_len-reference_covered);
							first_vector.push_back(left_op);
							reference_covered += (first_ref_len-reference_covered);

//							std::cout << "Ags: " << ref[reference_covered-4] << ref[reference_covered-3] << ref[reference_covered-2] << ref[reference_covered-1] << std::endl;
//							std::cout << ref[reference_covered] << ref[reference_covered+1] << ref[reference_covered+2] << ref[reference_covered+3] << ref[reference_covered+4] << ref[reference_covered+5] << ref[reference_covered+6] << ref[reference_covered+7] <<  ref[reference_covered+8] <<  std::endl;
							found_first_vector = true;
						}
						else if(reference_covered + count > first_ref_len) {
							is::CigarOp left_op = is::CigarOp(c.op, first_ref_len-reference_covered);
							first_vector.push_back(left_op);
							is::CigarOp right_op = is::CigarOp(c.op, (reference_covered + count)-first_ref_len);
							second_vector.push_back(right_op);

							reference_covered += (first_ref_len-reference_covered);

							found_first_vector = true;
//							std::cout << "Ags: " << ref[reference_covered-4] << ref[reference_covered-3] << ref[reference_covered-2] << ref[reference_covered-1] << std::endl;
//							std::cout << ref[reference_covered] << ref[reference_covered+1] << ref[reference_covered+2] << ref[reference_covered+3] << ref[reference_covered+4] << ref[reference_covered+5] << ref[reference_covered+6] << ref[reference_covered+7] <<  ref[reference_covered+8] <<  std::endl;
						} else {
							first_vector.push_back(c);
							reference_covered += count;
						}
					}

					if(c.op == 'I') {
						read_covered += count;
						first_vector.push_back(c);
					}
				}

				if(front_clipping > 0 && complete_cigar.size() == 0) {
					complete_cigar.push_back(is::CigarOp('S', front_clipping));
				}

				int size_temp = 0;
				for (auto& c: first_vector) {
					complete_cigar.push_back(c);
					if(c.op == '=' || c.op == 'X' || c.op == 'I') {
						size_temp += c.count;
					}
				}
				int gap_length = ((((current_exon.start - previous_exon.stop) + first_offset_left) - first_offset_right) - second_offset_left) + second_offset_right;

				complete_cigar.push_back(is::CigarOp('N', gap_length));

				container_vector = second_vector;

				if (index == merged_exons.size() - 1) {

					int size_temp = 0;
					for (auto& c: second_vector) {
						if(c.op == '=' || c.op == 'X' || c.op == 'I') {
							size_temp += c.count;
						}
						complete_cigar.push_back(c);
					}

					if (back_clipping > 0) {
						complete_cigar.push_back(is::CigarOp('S', back_clipping));
					}

					index += 1;
					continue;
				}

				previous_exon = current_exon;

				ref = "";

				if (read_covered > total_read.size()) {
					std::vector<is::CigarOp> temp;
//					std::cout << "bad exit 2" << std::endl;
					return temp;
				}

				read = total_read.substr(read_covered, total_read.size() - read_covered);

				first_offset_left = base_offset + std::max(-previous_exon.rightOffset, 0);
				first_offset_right = base_offset + std::max(previous_exon.rightOffset, 0);


				for (int var = previous_exon.start - second_offset_left; var < previous_exon.start; ++var) {
					ref.push_back(ref_data[var]);
				}

				for (int var = previous_exon.start + second_offset_right; var < previous_exon.stop - first_offset_left; ++var) {
					ref.push_back(ref_data[var]);
				}

				if(first_offset_right > 0) {
					for (int var = previous_exon.stop - first_offset_left; var < previous_exon.stop + first_offset_right; ++var) {
						ref.push_back(ref_data[var]);
					}
				}

				first_ref_len = (((previous_exon.reference.size() - first_offset_left) + first_offset_right) + second_offset_left) - second_offset_right;

				shouldAling = previous_exon.rightOffset != 0;
			} else {
				if (container_vector.size() > 0) {
					int size_temp = 0;
					for (auto& c: container_vector) {
						complete_cigar.push_back(c);
						if(c.op == '=' || c.op == 'X' || c.op == 'I') {
							size_temp += c.count;
						}
					}
				} else {
					int size_temp = 0;
					for (auto& c: previous_exon.cigar) {
						complete_cigar.push_back(c);
						if(c.op == '=' || c.op == 'X' || c.op == 'I') {
							size_temp += c.count;
						}
					}
				}

				int gap_length = current_exon.start - previous_exon.stop;
				complete_cigar.push_back(is::CigarOp('N', gap_length));

				container_vector.clear();

				if (index == merged_exons.size() - 1) {
					int size_temp = 0;
					for (auto& c: current_exon.cigar) {
						complete_cigar.push_back(c);
						if(c.op == '=' || c.op == 'X' || c.op == 'I') {
							size_temp += c.count;
						}
					}
					index += 1;
					continue;
				}

				previous_exon = current_exon;

				ref = "";

				read = previous_exon.content;

				first_offset_left = base_offset + std::max(-previous_exon.rightOffset, 0);
				first_offset_right = base_offset + std::max(previous_exon.rightOffset, 0);

				for (int var = previous_exon.start; var < previous_exon.stop - first_offset_left; ++var) {
					ref.push_back(ref_data[var]);
				}

				if(first_offset_right > 0) {
					for (int var = previous_exon.stop - first_offset_left; var < previous_exon.stop + first_offset_right; ++var) {
						ref.push_back(ref_data[var]);
					}
				}

				first_ref_len = (previous_exon.reference.size() - first_offset_left) + first_offset_right;

				shouldAling = previous_exon.rightOffset != 0;
			}
		}
		index += 1;
	}

	if (did_adjust_read) {
		return complete_cigar;
	} else {
		std::vector<is::CigarOp> temp;
		return temp;
	}
}

void FindEquivavlencyClasses(std::vector<std::vector<int> >* classes, std::vector<std::vector<int>> *equivalency_matrix, std::vector<int>* visited_nodes, int k, int n) {
	for(int i = k+1; i < n; i++) {
//		std::cout << "i " << i << std::endl;
		if((*equivalency_matrix)[k][i] == 1) {
			if((*visited_nodes)[i] != 1) {
				(*classes)[classes->size()-1].push_back(i);
				(*visited_nodes)[i] = 1;
				FindEquivavlencyClasses(classes, equivalency_matrix, visited_nodes, i, n);
			}
		}
	}

//	std::cout << "returning FindEquivavlencyClasses" << std::endl;
}

void GraphMap::PostprocessRNAData(std::vector<RealignmentStructure *> realignment_structures, std::vector<std::string> *sam_lines, int64_t num_threads, ProgramParameters *parameters, EValueParams *evalue_params) {
	std::shared_ptr<is::MinimizerIndex> first_index = indexes_[0];
	std::vector<int64_t*> coverages_array;

	int number_of_refs = first_index->get_reference_lengths().size() / 2;

	for (int i = 0; i < number_of_refs; ++i) {
	  int64_t ref_data_start = first_index->get_reference_starting_pos()[i];
	  int64_t ref_data_len = first_index->get_reference_lengths()[i];

	  int64_t *coverage_array = new int64_t[ref_data_len];
	  for (int var = 0; var < ref_data_len; ++var) {
		  coverage_array[var] = 0;
	  }
	  coverages_array.push_back(coverage_array);
	}

	std::vector<std::string> ref_names;
	std::vector<std::vector<RealignmentStructure*>> reads_by_chromosome;

	for (int i = 0; i < number_of_refs; ++i) {
	  std::string ref_name = indexes_.front()->get_headers()[i];
	  ref_names.push_back(ref_name);
	  std::vector<RealignmentStructure *> vector;
	  reads_by_chromosome.push_back(vector);
	}

	for (int var = 0; var < realignment_structures.size(); ++var) {
	  RealignmentStructure *rs = realignment_structures[var];
	  if (rs->ref_number >= 0) {
		  reads_by_chromosome[rs->ref_number].push_back(rs);
	  }
	}

	for (int var = 0; var < reads_by_chromosome.size(); ++var) {
	  std::sort(reads_by_chromosome[var].begin(), reads_by_chromosome[var].end(), GraphMap::comparePtrToNode);
	}

	std::vector<std::vector<RealignmentStructure *>> realignment_clusters;
	std::vector<RealignmentStructure *> current_realignment_cluster;

	for (int var = 0; var < reads_by_chromosome.size(); ++var) {
		std::vector<RealignmentStructure *> current_realignment_chromosome = reads_by_chromosome[var];

		if (current_realignment_chromosome.size() <= 0) {
			continue;
		}

		current_realignment_cluster.push_back(current_realignment_chromosome[0]);

		int current_structure_index = 0;
		int next_structure_index = 1;

		while(next_structure_index < current_realignment_chromosome.size()) {
		  RealignmentStructure* current_structure = current_realignment_chromosome[current_structure_index];
		  RealignmentStructure* next_structure = current_realignment_chromosome[next_structure_index];

		  if (next_structure_index == current_structure_index) {
			  realignment_clusters.push_back(current_realignment_cluster);
			  current_realignment_cluster.clear();
			  current_realignment_cluster.push_back(current_structure);
			  next_structure_index += 1;
		  } else if (current_structure->stop < next_structure->start) {
			  current_structure_index += 1;
		  } else {
			  current_realignment_cluster.push_back(next_structure);
			  next_structure_index += 1;
		  }
		}
		realignment_clusters.push_back(current_realignment_cluster);
		current_realignment_cluster.clear();
	}

	for (int var = 0; var < realignment_clusters.size(); ++var) {
	  std::vector<RealignmentStructure *> current_realignment_cluster = realignment_clusters[var];

//	  std::ofstream myfile;
//	  myfile.open ("clusters "+ std::to_string(var) +".txt");

	  for (int i = 0; i < current_realignment_cluster.size(); ++i) {
		RealignmentStructure* realignment_structure = current_realignment_cluster[i];
		const SingleSequence *read = realignment_structure->sequence;

//		myfile << read->get_header() << std::endl;
	  }

//	  myfile.close();
	}

	#pragma omp parallel for num_threads(1) firstprivate(evalue_params) shared(parameters, sam_lines) schedule(dynamic, 1)
	for (int var = 0; var < realignment_clusters.size(); ++var) {
	  std::vector<RealignmentStructure *> current_realignment_cluster = realignment_clusters[var];

	  uint32_t thread_id = omp_get_thread_num();

	  int min_index = INT_MAX;
	  int max_index = 0;
	  int ref_number = 0;

	  int halvening = 0;
	  for(int i = 0; i < first_index->get_reference_lengths().size()/2; i++) {
		  int len = first_index->get_reference_lengths()[i];
		  halvening += first_index->get_reference_lengths()[i];
	  }

	  for (int j = 0; j < current_realignment_cluster.size(); ++j) {
		RealignmentStructure* rs = current_realignment_cluster[j];

		if(rs->ref_number < 0) {
		 continue;
		}

		ref_number = rs->ref_number;

		int64_t ref_data_len = first_index->get_reference_lengths()[rs->ref_number];


		int8_t *ref_data_int  = (int8_t *) &first_index->get_data()[0];
		const char *ref_data = (const char *) ref_data_int;

		int len_total = 0;

		int64_t *coverage_array = coverages_array[rs->ref_number];
		int counter = std::max(7, rs->start);

		if (rs->start < min_index) {
			min_index = std::max(0, rs->start);
		}

		if (rs->stop > max_index) {
			max_index = std::min((int64_t) rs->stop, ref_data_len);
		}

		if(ref_data_len < 50) {
			continue;
		}

		int gapCounter = 0;

		for (int i = 0; i < rs->raw_alignment.size(); ++i) {
		  char align_op = 255;
		  align_op = rs->raw_alignment[i];

		  if (rs->raw_alignment[i] != 1 && rs->raw_alignment[i] != 4) {
			if (rs->raw_alignment[i] != 7) {
			  if (gapCounter > 20) {
				coverage_array[counter-1] += 1;
				coverage_array[counter-2] += 1;
				coverage_array[counter-3] += 1;
				coverage_array[counter-4] += 1;
				coverage_array[counter-5] += 1;
				coverage_array[counter-6] += 1;
				coverage_array[counter-7] += 1;
			  }
			  gapCounter = 0;
			  coverage_array[counter] += 1;
			} else {
				  gapCounter += 1;
				  if (gapCounter < 8) {
					coverage_array[counter] += 1;
				  }
			}
			counter += 1;
		  } else {
			if (gapCounter > 20) {
			  coverage_array[counter-1] += 1;
			  coverage_array[counter-2] += 1;
			  coverage_array[counter-3] += 1;
			  coverage_array[counter-4] += 1;
			  coverage_array[counter-5] += 1;
			  coverage_array[counter-6] += 1;
			  coverage_array[counter-7] += 1;
			}
			gapCounter = 0;
		  }
		}
	  }

	  int64_t offset_ref = 0;

	  for (int i = 0; i < first_index->get_reference_lengths().size(); ++i) {
		if (ref_number == i) {
		  break;
		}
		int64_t ref_data_len = first_index->get_reference_lengths()[i];
		offset_ref += ref_data_len;
	  }

	  int current_start = min_index;
	  int current_stop = min_index;

	  std::vector<ExonMatch> covered_regions;

	  int sum = 0;
	  for (int i = min_index; i < max_index; ++i) {
		if (coverages_array[ref_number][i] > 4) {
		  sum += coverages_array[ref_number][i];
		  if (current_stop+10 < i) {
		    ExonMatch em = ExonMatch();
			em.start = current_start + offset_ref;
			em.stop = current_stop + offset_ref;
			em.coverage = sum;
			sum = 0;
			current_start = i;
			current_stop = i;
			if ((em.stop - em.start) > 4) {
			  covered_regions.push_back(em);
			}
		  } else {
			current_stop = i;
		  }
		}
	  }

	  ExonMatch em = ExonMatch();
	  em.start = current_start + offset_ref;
	  em.stop = current_stop + offset_ref;
	  em.coverage = sum;
	  covered_regions.push_back(em);

	  ExonsCluster exon_cluster = ExonsCluster(covered_regions);

	  if(covered_regions.size() <= 0) {
		continue;
	  }

	  int8_t *ref_data_int  = (int8_t *) &first_index->get_data()[0];
	  const char *ref_data = (const char *) ref_data_int;
	  int64_t ref_data_lenInner = first_index->get_reference_lengths()[0];

//	  if(!exon_cluster.isValid() || exon_cluster.exons.size() <= 0) {
//		  continue;
//	  }

	  int start = exon_cluster.exons[0].start;
	  int stop = exon_cluster.exons[exon_cluster.exons.size()-1].stop;

	  std::string cutted_reference;

	  for (int i = start; i < stop+1; ++i) {
		  cutted_reference += ref_data[i];
	  }

//	  if (cutted_reference.size() > 10000) {
//		continue;
//	  }

	  exon_cluster.cuttedRefStart = 0;
	  exon_cluster.cuttedRefEnd = cutted_reference.size();

	  std::vector<ExonInfo> exons;

//	  std::cout << "current_realignment_cluster.size() " << current_realignment_cluster.size() << std::endl;

	  int number_of_reads_processed = 0;

	  for (int i = 0; i < current_realignment_cluster.size(); ++i) {

		RealignmentStructure* realignment_structure = current_realignment_cluster[i];
		const SingleSequence *read = realignment_structure->sequence;

		double old_score = realignment_structure->score;

		SeqOrientation orientation = realignment_structure->orientation;
		std::string sam_line_realigned = "";
		auto mapping_data_realing = std::unique_ptr<MappingData>(new MappingData);

		double score = -10000;
		SingleSequence *ss = new SingleSequence();
		ss->CopyFrom(*read);

		std::vector<CigarExon> new_cigar_exons;

//		if (orientation == kReverse) {
//		  ss->ReverseComplement();
//		  score = RealignRead(ss, indexes_[0], &(*mapping_data_realing), parameters, cutted_reference, exon_cluster, orientation, ref_number, &new_cigar_exons);
//		} else {
//		  score = RealignRead(read, indexes_[0], &(*mapping_data_realing), parameters, cutted_reference, exon_cluster, orientation, ref_number, &new_cigar_exons);
//		}

		if (true) {
			number_of_reads_processed += 1;

			std::string read_string = (const char*) read->get_data();

			int sumaRef = 0;
			int suma = 0;

			int order_number = 0;

			std::vector<CigarExon> previousCigarExons;

			if (orientation == kReverse) {
				for(int i = realignment_structure->previousCigarExons.size()-1; i >= 0; i--) {
					previousCigarExons.push_back(realignment_structure->previousCigarExons[i]);
				}

				ss->ReverseComplement();
				read_string = (const char*) ss->get_data();

			} else {
				for(CigarExon &ce: realignment_structure->previousCigarExons) {
					previousCigarExons.push_back(ce);
				}
			}

			bool isEdgeExonSet = false;

			for(CigarExon &ce: previousCigarExons) {
				if(!ce.isGap) {
					int length = 0;
					int lengthRef = 0;
					for (is::CigarOp &op: ce.cigar) {
						if (op.op != 'N' && op.op != 'D') {
							length += op.count;
						}
						if (op.op != 'I' && op.op != 'S') {
							lengthRef += op.count;
						}
					}

					ExonInfo ei = ExonInfo(ce, order_number, false, 0, 0);

					if(!isEdgeExonSet) {
						ei.isStartExon = true;
						isEdgeExonSet = true;
					}

					order_number += 1;

					long location = realignment_structure->raw_start;

					if(location > halvening) {

						location = realignment_structure->raw_stop;

						long buffer_offset = 0;
						long desired_index = 0;
						bool found_index = false;

						for(int i = 0; i < first_index->get_reference_lengths().size()/2; i++) {
							long len = first_index->get_reference_lengths()[i];
							if(!found_index) {
								if(buffer_offset + len < (location-halvening)) {
									buffer_offset += len;
								} else {
									found_index = true;
									desired_index = i;
								}
							}
						}

						location = buffer_offset + (first_index->get_reference_lengths()[desired_index] - ((location - halvening) - buffer_offset));
					}

					std::string cut_reff;
					for (int iks = (sumaRef + location); iks < (sumaRef + location+lengthRef); ++iks) {
					  cut_reff += ref_data[iks];
					}

					ei.start = sumaRef + location;
					ei.stop = sumaRef + location + lengthRef;

					ei.content = read_string.substr(suma, length);
					ei.reference = cut_reff;
					ei.read_id = read->get_header();
					exons.push_back(ei);

					suma += length;
					sumaRef += lengthRef;
				} else {
					for (is::CigarOp &op: ce.cigar) {
						sumaRef += op.count;
					}
				}
			}

			exons[exons.size()-1].isEndExon = true;
		} else {
		}

		if (mapping_data_realing->final_mapping_ptrs.size() <= 0) {
			continue;
		}

//		std::vector<CigarExon> previous_cigar_exons = realignment_structure->previousCigarExons;

//		PathGraphEntry* entry = mapping_data_realing->final_mapping_ptrs.back();
//		std::vector<AlignmentResults> alignments = entry->get_alignments();
//		if (false) {
////			if (alignments.size() > 0) {
//		  AlignmentResults ar = alignments.back();
//		  if (orientation == kReverse) {
//			CollectAlignments(ss, parameters, &(*mapping_data_realing), sam_line_realigned);
//		  } else {
//			CollectAlignments(realignment_structure->sequence, parameters, &(*mapping_data_realing), sam_line_realigned);
//		  }
//		  mapping_data_realing->clear();
//
//		  double difference = abs(score-old_score);
//		  if (old_score > score+0.05) {
//		  } else if(difference < 0.05) {
//			if (previous_cigar_exons.size() != new_cigar_exons.size()) {
//			} else {
//			  std::vector<int> visitingGaps;
//
//			  int numberOfGaps = 0;
//			  for (int varGaps = 0; varGaps < previous_cigar_exons.size(); varGaps++) {
//				  if (previous_cigar_exons[varGaps].isGap) {
//					  numberOfGaps += 1;
//				  }
//			  }
//
//			  for (int varz = 0; varz < numberOfGaps; ++varz) {
//				  visitingGaps.push_back(0);
//			  }
//
//			  int currentGap = -1;
//			  int firstCigarSum = 0;
//			  int secondCigarSum = 0;
//
//			  for (int vary = 0; vary < previous_cigar_exons.size()-1; ++vary) {
//				  CigarExon c1 = previous_cigar_exons[vary];
//				  CigarExon c2 = new_cigar_exons[vary];
//
//				  if (abs(firstCigarSum-secondCigarSum) > 5) {
//					  visitingGaps[currentGap] = 1;
//				  }
//
//				  if (!c1.isGap) {
//					  currentGap += 1;
//				  }
//
//				  firstCigarSum += c1.length;
//				  secondCigarSum += c2.length;
//
//				  if (abs(firstCigarSum-secondCigarSum) > 5) {
//					  visitingGaps[currentGap] = 1;
//				  }
//			  }
//
//			  currentGap = 0;
//			  CigarExon previousExon1 = previous_cigar_exons[0];
//			  CigarExon previousExon2 = new_cigar_exons[0];
//
//			  bool isVisiting = false;
//			  for (int varvisiting = 0; varvisiting < visitingGaps.size(); ++varvisiting) {
//				  if (visitingGaps[varvisiting]) {
//					  isVisiting = true;
//					  break;
//				  }
//			  }
//
//			  int total1Sum = 0;
//			  int total2Sum = 0;
//
//			  for (int vary = 1; vary < previous_cigar_exons.size()-1; ++vary) {
//				  CigarExon currentExon1 = previous_cigar_exons[vary];
//				  CigarExon currentExon2 = new_cigar_exons[vary];
//
//				  if (!currentExon1.isGap) {
//					  previousExon1 = currentExon1;
//					  previousExon2 = currentExon2;
//					  continue;
//				  }
//
//				  CigarExon nextExon1 = previous_cigar_exons[vary+1];
//				  CigarExon nextExon2 = new_cigar_exons[vary+1];
//
//				  if (visitingGaps[currentGap]) {
//					  int sum1 = 0;
//					  int count1 = 0;
//					  for (auto& c: nextExon1.cigar) {
//							int numberOfBases = (int) c.count;
//							if (numberOfBases + count1 > 10) {
//								if (c.op == '=') {
//									sum1 += (10-count1) * 5;
//								} else {
//									sum1 += (10-count1) * -4;
//								}
//								break;
//							} else {
//								if (c.op == '=') {
//									sum1 += numberOfBases * 5;
//								} else {
//									sum1 += numberOfBases * -4;
//								}
//								count1 += numberOfBases;
//							}
//					  }
//
//					  int sum1Back = 0;
//					  int count1Back = 0;
//					  for (int unutar1 = previousExon1.cigar.size()-1; unutar1 >= 0; unutar1--) {
//
//							int numberOfBases = (int) previousExon1.cigar[unutar1].count;
//							if (numberOfBases + count1Back > 10) {
//								if (previousExon1.cigar[unutar1].op == '=') {
//									sum1Back += (10-count1Back) * 5;
//								} else {
//									sum1Back += (10-count1Back) * -4;
//								}
//								break;
//							} else {
//								if (previousExon1.cigar[unutar1].op == '=') {
//									sum1Back += numberOfBases * 5;
//								} else {
//									sum1Back += numberOfBases * -4;
//								}
//								count1Back += numberOfBases;
//							}
//					  }
//
//					  int sum2 = 0;
//					  int count2 = 0;
//					  for (auto& c: nextExon2.cigar) {
//							int numberOfBases = (int) c.count;
//							if (numberOfBases + count2 > 10) {
//								if (c.op == '=') {
//									sum2 += (10-count2) * 5;
//								} else {
//									sum2 += (10-count2) * -4;
//								}
//								break;
//							} else {
//								if (c.op == '=') {
//									sum2 += numberOfBases * 5;
//								} else {
//									sum2 += numberOfBases * -4;
//								}
//								count2 += numberOfBases;
//							}
//					  }
//
//					  int sum2Back = 0;
//					  int count2Back = 0;
//					  for (int unutar1 = previousExon2.cigar.size()-1; unutar1 >= 0; unutar1--) {
//
//							int numberOfBases = (int) previousExon2.cigar[unutar1].count;
//							if (numberOfBases + count2Back > 10) {
//								if (previousExon2.cigar[unutar1].op == '=') {
//									sum2Back += (10-count2Back) * 5;
//								} else {
//									sum2Back += (10-count2Back) * -4;
//								}
//								break;
//							} else {
//								if (previousExon2.cigar[unutar1].op == '=') {
//									sum2Back += numberOfBases * 5;
//								} else {
//									sum2Back += numberOfBases * -4;
//								}
//								count2Back += numberOfBases;
//							}
//					  }
//
//					  total1Sum += (sum1 + sum1Back);
//					  total2Sum += (sum2 + sum2Back);
//				  }
//
//				  currentGap += 1;
//			  }
//
//			  if (!isVisiting) {
//				  if (score >= old_score) {
//					  #pragma omp critical
//					  if (orientation == kReverse) {
////						  (*sam_lines)[realignment_structure->order_number] = "";
//					  } else {
////						  (*sam_lines)[realignment_structure->order_number] = sam_line_realigned;
//					  }
//				  }
//			  } else if (total1Sum <= total2Sum) {
//				  #pragma omp critical
//				  if (orientation == kReverse) {
////					  (*sam_lines)[realignment_structure->order_number] = "";
//				  } else {
////					  (*sam_lines)[realignment_structure->order_number] = sam_line_realigned;
//				  }
//			  }
//			}
//		  }
//		}
		delete ss;
	  }

//	  std::cout << "Number of exons " << exons.size() << std::endl;
	  std::sort(std::begin(exons), std::end(exons), [](ExonInfo a, ExonInfo b) {
		  if(a.start == b.start) {
			  return a.stop < b.stop;
		  } else {
			  return a.start < b.start;
		  }
	  });
	  std::vector<std::vector<ExonInfo>> clusters_of_exons_mid;

	  if(exons.size() > 0) {
		  std::vector<ExonInfo> current_exons;
		  long current_stop = exons[0].stop;

		  for(ExonInfo &exon: exons) {
			 if(current_stop < exon.start) {
				 clusters_of_exons_mid.push_back(current_exons);
				 current_exons.clear();
				 current_exons.push_back(exon);
				 current_stop = exon.stop;
			 } else {
				 current_exons.push_back(exon);
				 current_stop = std::max(exon.stop, current_stop);
			 }
		  }

		  clusters_of_exons_mid.push_back(current_exons);
	  }

	  std::vector<std::vector<ExonInfo>> clusters_of_exons;

	  for(std::vector<ExonInfo> &cluster_of_exons: clusters_of_exons_mid) {

		  if(cluster_of_exons.size() > 4000) {
			  continue;
		  }

		  std::vector<std::vector<int>> equivalency_matrix;

		  for(int i = 0; i < cluster_of_exons.size(); i++) {
			  std::vector<int> tmp_vector;
			  for (int j = 0; j < cluster_of_exons.size(); j++) {
				  tmp_vector.push_back(0);
			  }
			  equivalency_matrix.push_back(tmp_vector);
		  }

		  for(int i = 0; i < cluster_of_exons.size()-1; i++) {
			  ExonInfo ei1 = cluster_of_exons[i];
			  for(int j = i+1; j < cluster_of_exons.size(); j++) {
				  ExonInfo ei2 = cluster_of_exons[j];
				  if(abs(ei1.start - ei2.start) < 40 && (ei1.stop - ei2.stop) < 40) {
					  equivalency_matrix[i][j] = 1;
				  }
			  }
		  }

		  std::vector<std::vector<int>> classes;
		  std::vector<int> equivalency_class;
		  equivalency_class.push_back(0);
		  classes.push_back(equivalency_class);

		  std::vector<int> visited_nodes;

		  for(auto node: cluster_of_exons) {
			  visited_nodes.push_back(0);
		  }

		  FindEquivavlencyClasses(&classes, &equivalency_matrix, &visited_nodes, 0, cluster_of_exons.size());

		  for(int i = 1; i < cluster_of_exons.size(); i++) {
			  if(visited_nodes[i] == 0) {
				  std::vector<int> equivalency_class;
				  equivalency_class.push_back(i);
				  classes.push_back(equivalency_class);
				  FindEquivavlencyClasses(&classes, &equivalency_matrix, &visited_nodes, i, cluster_of_exons.size());
			  }
		  }

		  for(auto &eq_class: classes) {
			  std::vector<ExonInfo> cluster;
			  for(int index: eq_class) {
				  cluster.push_back(cluster_of_exons[index]);
			  }
			  clusters_of_exons.push_back(cluster);
		  }
	  }

	  for(std::vector<ExonInfo> &exon_cluster: clusters_of_exons) {

//		for(ExonInfo ei: exon_cluster) {
//			std::cout << ei.start << "-" << ei.stop << "  " << ei.isStartExon << " " << ei.isEndExon << " (" << ei.invalid << ")" << std::endl;
//		}
//		std::cout << std::endl;

	    	 if (exon_cluster.size() / (double) number_of_reads_processed > 0.1) {
	    		 std::vector<std::string> sequences;
	    		 std::sort(std::begin(exon_cluster), std::end(exon_cluster), [](ExonInfo a, ExonInfo b) { return a.reference.size() > b.reference.size(); });

	    		 long minLocation = LONG_MAX;
			 long maxLocation = 0;
			 for(ExonInfo &einfo: exon_cluster) {
				 if(einfo.stop > maxLocation) {
					 maxLocation = einfo.stop;
				 }
				 if(einfo.start < minLocation) {
					 minLocation = einfo.start;
				 }
			 }

			 int coverageSize = (int) maxLocation - minLocation;
			 int coverageArray[coverageSize];

			 int coverageArray2[coverageSize];
			 int coverageArray3[coverageSize];

			 for (int location = 0; location < coverageSize; ++location) {
				 coverageArray[location] = 0;
				 coverageArray2[location] = 0;
				 coverageArray3[location] = 0;
			 }

			 for(ExonInfo &einfo: exon_cluster) {
				 for (long location = einfo.start; location < einfo.stop; ++location) {
					 long value = location - minLocation;
					 coverageArray[(int)value] += 1;
				 }
			 }

			 int max_coverage = exon_cluster.size();
			 std::vector<int> start_pivots;
			 std::vector<int> end_pivots;

	    	 	 for (int location = 0; location < coverageSize; ++location) {
	    	 		 if((coverageArray2[location] / (double) max_coverage) > 0.1) {
	    	 			 start_pivots.push_back(location);
	    	 		 }
	    	 		 if((coverageArray3[location] / (double) max_coverage) > 0.1) {
    	 				 end_pivots.push_back(location);
	    	 		 }
	    	 	 }

			 bool isEnd = false;
		    	 bool isStart = true;
		    	 int startOffset = 0;
		    	 int endOffset = 0;

		    	 // checkpoint podesit parametre
		    	 for (int location = 0; location < coverageSize; ++location) {
		    		 if (isStart) {
		    			 if ( (coverageArray[location] / (double) max_coverage) < 0.5) {
		    				 startOffset += 1;
		    			 } else {
		    				 isStart = false;
		    			 }
		    		 } else if (isEnd) {
		    			 endOffset += 1;
		    		 } else {
		    			 if ( (coverageArray[location] / (double) max_coverage) < 0.5) {
		    				 endOffset += 1;
		    				 isEnd = true;
		    			 }
		    		 }
		    	 }

		    	 for(ExonInfo &einfo: exon_cluster) {

//			    	 int startOffset = 0;
//			    	 int endOffset = 0;
//
//			    	 int distance_start = max_coverage + 1;
//
//			    	 for(int pivot: start_pivots) {
//			    		 int current_location = (int)(einfo.start-minLocation);
//			    		 int new_distance = abs(pivot - current_location);
//			    		 if(new_distance < distance_start) {
//			    			 distance_start = new_distance;
//			    			 startOffset = pivot;
//			    		 }
//			    	 }
//
//			    	 int distance_end = max_coverage + 1;
//
//			    	 for(int pivot: end_pivots) {
//			    		 int current_location = (int)((einfo.stop-minLocation)-1);
////			    		 std::cout << current_location << std::endl;
//			    		 int new_distance = abs(pivot - current_location);
////			    		 std::cout << new_distance << std::endl;
//			    		 if(new_distance < distance_end) {
//			    			 distance_end = new_distance;
//			    			 endOffset = (coverageSize-1) - pivot;
//			    		 }
//			    	 }
//
//			    	 std::cout << "startOffset " << startOffset << " , endOffset " << endOffset << std::endl;
////			    	 std::cout << einfo.start << " " << einfo.stop << std::endl;
//
//			    	 if(!einfo.isStartExon) {
//			    		 einfo.leftOffset = (einfo.start - minLocation) - startOffset;
//			    	 } else {
//			    		 einfo.leftOffset = 0;
//			    	 }
//			    	 if(!einfo.isEndExon) {
//			    		 einfo.rightOffset = (maxLocation - einfo.stop) - endOffset;
//			    	 } else {
//			    		 einfo.rightOffset = 0;
//			    	 }

		    		 einfo.leftOffset = (einfo.start - minLocation) - startOffset;
		    		 einfo.rightOffset = (maxLocation - einfo.stop) - endOffset;
		    	 }
	    	 } else {
			 for(ExonInfo &einfo: exon_cluster) {
				if(einfo.stop - einfo.start < 15) {
					einfo.invalid = true;
				}
			 }
	    	 }
	  }

	  std::map<std::string, std::vector<ExonInfo>> map_of_exons;

	  for(std::vector<ExonInfo> &exon_cluster: clusters_of_exons) {
		for(ExonInfo &einfo: exon_cluster) {
			  if (map_of_exons.find(einfo.read_id) != map_of_exons.end()) {
				  std::map<std::string, std::vector<ExonInfo>>::iterator i = map_of_exons.find(einfo.read_id);
				  i->second.push_back(einfo);
			  } else {
				  std::vector<ExonInfo> read_exons;
				  read_exons.push_back(einfo);
				  map_of_exons[einfo.read_id] = read_exons;
			  }
		}
	  }

	  for (auto x : map_of_exons) {
		 std::vector<is::CigarOp> rez = ProcessReadExons(x.second, ref_data);
		 if (rez.size() > 0 && x.second.size() > 0) {

			  for (int i = 0; i < current_realignment_cluster.size(); ++i) {
				  RealignmentStructure* realignment_structure = current_realignment_cluster[i];

				  if (realignment_structure->sequence->get_header() == x.second[0].read_id) {
					  std::string sam_line_test = "";

					  auto mapping_data_test = std::unique_ptr<MappingData>(new MappingData);

					  SingleSequence *ss = new SingleSequence();
					  ss->CopyFrom(*realignment_structure->sequence);
					  ss->ReverseComplement();

					  bool is_aligned = GetMappingData(realignment_structure, indexes_[0], &(*mapping_data_test), parameters, rez, ss);

					  CollectAlignments(realignment_structure->sequence, parameters, &(*mapping_data_test), sam_line_test);

					  if (is_aligned) {
						  (*sam_lines)	[realignment_structure->order_number] = sam_line_test;
					  } else {
					  }

					  delete ss;

					  continue;
				  }
			  }
		}
	  }
	}

	for (int var = 0; var < number_of_refs; ++var) {
	  delete [] coverages_array[var];
	}
}

int GraphMap::ProcessSequenceFileInParallel(ProgramParameters *parameters, const SequenceFile *reads, clock_t *last_time, FILE *fp_out, int64_t *ret_num_mapped, int64_t *ret_num_unmapped) {
  int64_t num_reads = reads->get_sequences().size();
  std::vector<std::string> sam_lines;

  if (parameters->output_in_original_order == true) {
    sam_lines.resize(num_reads, std::string(""));
  }

  // Division by to to avoid hyperthreading cores, and limit on 24 to avoid clogging a shared SMP.
//  int64_t num_threads = std::min(24, ((int) omp_get_num_procs()) / 2);
//
//  if (parameters->num_threads > 0)
    int64_t num_threads = 12;

//  int64_t num_threads = (int64_t) parameters->num_threads;
  LogSystem::GetInstance().Log(VERBOSE_LEVEL_HIGH | VERBOSE_LEVEL_MED, true, FormatString("Using %ld threads.", num_threads), "ProcessReads");

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

  std::vector<RealignmentStructure *> realignment_structures;

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
//        ss << "\n";
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

    int mapped_state = STATE_UNMAPPED;
    std::string sam_line = "";
    // The actual interesting part.
    auto mapping_data = std::unique_ptr<MappingData>(new MappingData);
    ProcessRead(i, &(*mapping_data), reads->get_sequences()[i], parameters, evalue_params, &realignment_structures);

    // Generate the output.
    mapped_state = CollectAlignments(reads->get_sequences()[i], parameters, &(*mapping_data), sam_line);
    mapping_data->clear();

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
//    	  if(realignment_structures[realignment_structures.size()-1]->orientation == kReverse) {
//    	  }
    }
  }

  if (parameters->composite_parameters == "rnaseq") {
	  PostprocessRNAData(realignment_structures, &sam_lines, num_threads, parameters, evalue_params);
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
