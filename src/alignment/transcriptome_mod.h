/*
 * transcriptome_mod.h
 *
 *  Created on: Jan 5, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNMENT_TRANSCRIPTOME_MOD_H_
#define SRC_ALIGNMENT_TRANSCRIPTOME_MOD_H_

#include "index/index.h"
#include "index/index_spaced_hash_fast.h"
#include "containers/results.h"

int ConvertFromTranscriptomeToGenomeAln(const IndexSpacedHashFast *index, AlignmentResults *aln);

#endif /* SRC_ALIGNMENT_TRANSCRIPTOME_MOD_H_ */
