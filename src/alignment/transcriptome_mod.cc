/*
 * transcriptome_mod.cc
 *
 *  Created on: Jan 5, 2017
 *      Author: isovic
 */

#include "alignment/transcriptome_mod.h"

int ConvertFromTranscriptomeToGenomeAln(const IndexSpacedHashFast *index, AlignmentResults *aln) {
  if (index->is_transcriptome() == false) { return 1; }



  return 0;
}
