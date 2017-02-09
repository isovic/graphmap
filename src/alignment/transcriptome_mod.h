/*
 * transcriptome_mod.h
 *
 *  Created on: Jan 5, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNMENT_TRANSCRIPTOME_MOD_H_
#define SRC_ALIGNMENT_TRANSCRIPTOME_MOD_H_

//#include "index/index.h"
//#include "index/index_spaced_hash_fast.h"
#include "minimizer_index/minimizer_index.h"
#include "containers/results.h"
#include "program_parameters.h"
#include "graphmap/transcriptome.h"
#include <tuple>

int ConvertFromTranscriptomeToGenomeAln(const ProgramParameters *parameters, std::shared_ptr<is::MinimizerIndex> index, std::shared_ptr<is::Transcriptome> transcriptome, AlignmentResults *aln);

#endif /* SRC_ALIGNMENT_TRANSCRIPTOME_MOD_H_ */
