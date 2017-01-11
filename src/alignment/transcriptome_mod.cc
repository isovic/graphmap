/*
 * transcriptome_mod.cc
 *
 *  Created on: Jan 5, 2017
 *      Author: isovic
 */

#include "alignment/transcriptome_mod.h"
#include "alignment/cigargen.h"

int ConvertFromTranscriptomeToGenomeAln(const IndexSpacedHashFast *index, AlignmentResults *aln) {
  if (index->is_transcriptome() == false) { return 1; }

  std::vector<uint8_t> alignment;
  alignment.reserve(aln->alignment.size());

  int64_t pos_on_read = aln->query_start;
  int64_t pos_on_ref = aln->ref_start;

  for (int64_t i=0; i<aln->alignment.size(); i++) {
    int8_t op = aln->alignment[i];
    if (op == EDLIB_M || op == EDLIB_EQUAL || op == EDLIB_X) {
      pos_on_ref += 1;
      pos_on_read += 1;
    } else if (op == EDLIB_D) {
      pos_on_ref += 1;
    } else if (op == EDLIB_I || op == EDLIB_S) {
      pos_on_read += 1;
    }
  }

  return 0;
}
