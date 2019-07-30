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
		return temp;
	}

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


	for(ExonInfo exon: merged_exons) {
		std::string str(exon.read_id);
	}

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
				if (current_exon.cigar.size() > 0 && current_exon.cigar[current_exon.cigar.size()].op == 'S') {
					back_clipping = current_exon.cigar[current_exon.cigar.size()].count;
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

				// splitting cigar into two vectors
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
						}
						else if(reference_covered + count > first_ref_len) {
							is::CigarOp left_op = is::CigarOp(c.op, first_ref_len-reference_covered);
							first_vector.push_back(left_op);
							is::CigarOp right_op = is::CigarOp(c.op, (reference_covered + count)-first_ref_len);
							second_vector.push_back(right_op);

							read_covered += (first_ref_len-reference_covered);
							reference_covered += (first_ref_len-reference_covered);
							found_first_vector = true;
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
							found_first_vector = true;
						}
						else if(reference_covered + count > first_ref_len) {
							is::CigarOp left_op = is::CigarOp(c.op, first_ref_len-reference_covered);
							first_vector.push_back(left_op);
							is::CigarOp right_op = is::CigarOp(c.op, (reference_covered + count)-first_ref_len);
							second_vector.push_back(right_op);

							reference_covered += (first_ref_len-reference_covered);

							found_first_vector = true;
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
					front_clipping = 0;
				}

				for (auto& c: first_vector) {
					complete_cigar.push_back(c);
				}

				int gap_length = ((((current_exon.start - previous_exon.stop) + first_offset_left) - first_offset_right) - second_offset_left) + second_offset_right;
				complete_cigar.push_back(is::CigarOp('N', gap_length));


				container_vector = second_vector;

				if (index == merged_exons.size() - 1) {
					for (auto& c: second_vector) {
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
					for (auto& c: container_vector) {
						complete_cigar.push_back(c);
					}
				} else {
					for (auto& c: previous_exon.cigar) {
						complete_cigar.push_back(c);
					}
				}

				int gap_length = current_exon.start - previous_exon.stop;
				complete_cigar.push_back(is::CigarOp('N', gap_length));

				container_vector.clear();

				if (index == merged_exons.size() - 1) {
					for (auto& c: current_exon.cigar) {
						complete_cigar.push_back(c);
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
		if((*equivalency_matrix)[k][i] == 1) {
			if((*visited_nodes)[i] != 1) {
				(*classes)[classes->size()-1].push_back(i);
				(*visited_nodes)[i] = 1;
				FindEquivavlencyClasses(classes, equivalency_matrix, visited_nodes, i, n);
			}
		}
	}
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

//	for (int var = 0; var < realignment_clusters.size(); ++var) {
//	  std::vector<RealignmentStructure *> current_realignment_cluster = realignment_clusters[var];
//
//	  std::ofstream myfile;
//	  myfile.open ("clusters "+ std::to_string(var) +".txt");
//
//	  for (int i = 0; i < current_realignment_cluster.size(); ++i) {
//		RealignmentStructure* realignment_structure = current_realignment_cluster[i];
//		const SingleSequence *read = realignment_structure->sequence;
//
//		myfile << read->get_header() << std::endl;
//	  }
//
//	  myfile.close();
//	}

	#pragma omp parallel for num_threads(12) firstprivate(evalue_params) shared(parameters, sam_lines) schedule(dynamic, 1)
	for (int var = 0; var < realignment_clusters.size(); ++var) {
	  std::vector<RealignmentStructure *> current_realignment_cluster = realignment_clusters[var];


	  uint32_t thread_id = omp_get_thread_num();

	  int64_t min_index = INT_MAX;
	  int64_t max_index = 0;
	  int64_t ref_number = 0;

	  int64_t halvening = 0;
	  for(int64_t i = 0; i < first_index->get_reference_lengths().size()/2; i++) {
		  int64_t len = first_index->get_reference_lengths()[i];
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

		if (rs->start < min_index) {
			min_index = std::max((int64_t) 0, rs->start);
		}

		if (rs->stop > max_index) {
			max_index = std::min(rs->stop, ref_data_len);
		}
	  }

	  int current_start = min_index;
	  int current_stop = min_index;

	  int8_t *ref_data_int  = (int8_t *) &first_index->get_data()[0];
	  const char *ref_data = (const char *) ref_data_int;

	  std::vector<ExonInfo> exons;
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
		number_of_reads_processed += 1;

		std::string read_string = (const char*) read->get_data();

		int64_t suma_ref = 0;
		int64_t suma_read = 0;
		int64_t order_number = 0;

		std::vector<CigarExon> previousCigarExons;

		if (orientation == kReverse) {
			for(int i = realignment_structure->previousCigarExons.size()-1; i >= 0; i--) {
				CigarExon ce = realignment_structure->previousCigarExons[i];

				std::reverse(ce.cigar.begin(), ce.cigar.end());

				previousCigarExons.push_back(ce);
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
			  std::string qname = ((std::string) (read->get_header()));

			if(!ce.isGap) {
				int length_read = 0;
				int length_ref = 0;
				for (is::CigarOp &op: ce.cigar) {
					if (op.op != 'N' && op.op != 'D') {
						length_read += op.count;
					}
					if (op.op != 'I' && op.op != 'S') {
						length_ref += op.count;
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
					int64_t buffer_offset = 0;
					int64_t desired_index = 0;
					bool found_index = false;

					for(int j = 0; j < first_index->get_reference_lengths().size()/2; j++) {
						int64_t len = first_index->get_reference_lengths()[j];
						if(!found_index) {
							if(buffer_offset + len < (location-halvening)) {
								buffer_offset += len;
							} else {
								found_index = true;
								desired_index = j;
							}
						}
					}
					location = buffer_offset + (first_index->get_reference_lengths()[desired_index] - ((location - halvening) - buffer_offset));
				}

				std::string cut_reff;
				for (int j = (suma_ref + location); j < (suma_ref + location+length_ref); ++j) {
				  cut_reff += ref_data[j];
				}

				ei.start = suma_ref + location;
				ei.stop = suma_ref + location + length_ref;
				ei.content = read_string.substr(suma_read, length_read);
				ei.reference = cut_reff;
				ei.read_id = read->get_header();
				exons.push_back(ei);

				suma_read += length_read;
				suma_ref += length_ref;
			} else {
				for (is::CigarOp &op: ce.cigar) {
					suma_ref += op.count;
				}
			}
		}
		exons[exons.size()-1].isEndExon = true;
		delete ss;
	  }

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
		  int64_t current_stop = exons[0].stop;

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

		  for(int64_t i = 0; i < cluster_of_exons.size(); i++) {
			  std::vector<int> tmp_vector;
			  for (int64_t j = 0; j < cluster_of_exons.size(); j++) {
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
	    	 if (exon_cluster.size() / (double) number_of_reads_processed > 0.1) {
	    		 std::vector<std::string> sequences;
	    		 std::sort(std::begin(exon_cluster), std::end(exon_cluster), [](ExonInfo a, ExonInfo b) { return a.reference.size() > b.reference.size(); });

	    		 int64_t minLocation = LONG_MAX;
	    		 int64_t maxLocation = 0;
			 for(ExonInfo &einfo: exon_cluster) {
				 if(einfo.stop > maxLocation) {
					 maxLocation = einfo.stop;
				 }
				 if(einfo.start < minLocation) {
					 minLocation = einfo.start;
				 }
			 }

			 int64_t coverageSize = maxLocation - minLocation;

			 int64_t coverageArrayStart[coverageSize];
			 int64_t coverageArrayEnd[coverageSize];

			 for (int location = 0; location < coverageSize; ++location) {
				 coverageArrayStart[location] = 0;
				 coverageArrayEnd[location] = 0;
			 }

			 for(ExonInfo &einfo: exon_cluster) {
				 coverageArrayStart[einfo.start-minLocation] += 1;
				 coverageArrayEnd[einfo.stop-(1+minLocation)] += 1;
			 }

			 int max_coverage = exon_cluster.size();
			 std::vector<int> start_pivots;
			 std::vector<int> end_pivots;

	    	 	 for (int location = 0; location < coverageSize; ++location) {
	    	 		 if((coverageArrayStart[location] / (double) max_coverage) > 0.2) {
	    	 			 start_pivots.push_back(location);
	    	 		 }
	    	 		 if((coverageArrayEnd[location] / (double) max_coverage) > 0.2) {
    	 				 end_pivots.push_back(location);
	    	 		 }
	    	 	 }

	    	 	 std::vector<std::vector<int>> set_of_start_pivots;
	    	 	 std::vector<int> curent_start_pivots;

	    	 	 if(start_pivots.size() > 0) {
	    	 		int previous_pivot = start_pivots[0];

		    	 	 curent_start_pivots.push_back(previous_pivot);

		    	 	 for(int i = 1; i < start_pivots.size(); i++) {
		    	 		 int pivot = start_pivots[i];

		    	 		 if(pivot - previous_pivot < 6) {
		    	 			 curent_start_pivots.push_back(pivot);
		    	 			 previous_pivot = pivot;
		    	 		 } else {
		    	 			 set_of_start_pivots.push_back(curent_start_pivots);
		    	 			 curent_start_pivots.clear();
		    	 			 previous_pivot = pivot;
		    	 			curent_start_pivots.push_back(pivot);
		    	 		 }
		    	 	 }

		    	 	 set_of_start_pivots.push_back(curent_start_pivots);

		    	 	 start_pivots.clear();

		    	 	 for(std::vector<int> pivots_array: set_of_start_pivots) {
		    	 		 int max_pivot = pivots_array[0];
		    	 		 int max_value = coverageArrayStart[max_pivot];

		    	 		 for(int i = 1; i < pivots_array.size(); i++) {
		    	 			 int curr_pivot = pivots_array[i];
		    	 			 int max_value_curr = coverageArrayStart[curr_pivot];
		    	 			 if(max_value_curr > max_value) {
		    	 				 max_pivot = curr_pivot;
		    	 				 max_value = max_value_curr;
		    	 			 }
		    	 		 }

		    	 		 start_pivots.push_back(max_pivot);
		    	 	 }
	    	 	 }

	    	 	 std::vector<std::vector<int>> set_of_end_pivots;
	    	 	 std::vector<int> curent_end_pivots;

	    	 	 if(end_pivots.size() > 0) {
	    	 		 int previous_pivot = end_pivots[0];

	 	    	 	curent_end_pivots.push_back(previous_pivot);

	 	    	 	 for(int i = 1; i < end_pivots.size(); i++) {
	 	    	 		 int pivot = end_pivots[i];

	 	    	 		 if(pivot - previous_pivot < 6) {
	 	    	 			curent_end_pivots.push_back(pivot);
	 	    	 			previous_pivot = pivot;
	 	    	 		 } else {
	 	    	 			set_of_end_pivots.push_back(curent_end_pivots);
	 	    	 			curent_end_pivots.clear();
	 	    	 			previous_pivot = pivot;
	 	    	 			curent_end_pivots.push_back(pivot);
	 	    	 		 }
	 	    	 	 }

	 	    	 	set_of_end_pivots.push_back(curent_end_pivots);

	 			 end_pivots.clear();

	 			 for(std::vector<int> pivots_array: set_of_end_pivots) {
	 				 int max_pivot = pivots_array[0];
	 				 int max_value = coverageArrayEnd[max_pivot];

	 				 for(int i = 1; i < pivots_array.size(); i++) {
	 					 int curr_pivot = pivots_array[i];
	 					 int max_value_curr = coverageArrayEnd[curr_pivot];
	 					 if(max_value_curr > max_value) {
	 						 max_pivot = curr_pivot;
	 						 max_value = max_value_curr;
	 					 }
	 				 }

	 				 end_pivots.push_back(max_pivot);
	 			 }
	    	 	 }

		    	 for(ExonInfo &einfo: exon_cluster) {
			    	 int startOffset = 0;
			    	 int endOffset = 0;

			    	 int distance_start = max_coverage + 1;

			    	 for(int pivot: start_pivots) {
			    		 int current_location = (int) (einfo.start-minLocation);
			    		 int new_distance = abs(pivot - current_location);
			    		 if(new_distance < distance_start) {
			    			 distance_start = new_distance;
			    			 startOffset = pivot;
			    		 }
			    	 }

			    	 int distance_end = max_coverage + 1;

			    	 for(int pivot: end_pivots) {
			    		 int64_t current_location = (einfo.stop-minLocation)-1;
			    		 int64_t new_distance = abs(pivot - current_location);
			    		 if(new_distance < distance_end) {
			    			 distance_end = new_distance;
			    			 endOffset = (coverageSize-1) - pivot;
			    		 }
			    	 }

			    	 if(!einfo.isStartExon) {
			    		 einfo.leftOffset = (einfo.start - minLocation) - startOffset;
			    	 } else {
			    		 einfo.leftOffset = 0;
			    	 }
			    	 if(!einfo.isEndExon) {
			    		 einfo.rightOffset = (maxLocation - einfo.stop) - endOffset;
			    	 } else {
			    		 einfo.rightOffset = 0;
			    	 }

			    	 if(std::abs(einfo.leftOffset) > 25) {
			    		 einfo.leftOffset = 0;
			    	 }

			    	 if(std::abs(einfo.rightOffset) > 25) {
			    		 einfo.rightOffset = 0;
			    	 }
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
