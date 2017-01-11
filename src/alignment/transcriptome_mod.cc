/*
 * transcriptome_mod.cc
 *
 *  Created on: Jan 5, 2017
 *      Author: isovic
 */

#include "alignment/transcriptome_mod.h"
#include "alignment/cigargen.h"
#include "index/index_spaced_hash_fast.h"

int ConvertFromTranscriptomeToGenomeAln(const IndexSpacedHashFast *index, AlignmentResults *aln) {
  if (index->is_transcriptome() == false) { return 1; }

  auto& trans_id_to_regions = index->get_trans_id_to_regions();
  if (trans_id_to_regions.size() == 0) { return 2; }

  std::vector<uint8_t> alignment;
  alignment.reserve(aln->alignment.size());

  auto regions_it = trans_id_to_regions.find(aln->ref_header);
  if (regions_it == trans_id_to_regions.end()) { return 3; }
  auto& regions = regions_it->second;
  if (regions.size() == 0) { return 4; }

  int64_t pos_on_read = aln->query_start;
  int64_t pos_on_ref = aln->ref_start;

  int64_t reg_id = 0;
  int64_t region_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
  int64_t region_end = regions[reg_id].second - 1;      // They should be inclusive as well.

  for (int64_t i=0; i<aln->alignment.size(); i++) {
    int8_t op = aln->alignment[i];

    // The op needs to be copied regardless of whether there is a gap.
    alignment.push_back(op);

    // TODO: Or >=, depends whether the coordinates are 0 or 1 based.
    // Check if there should be a gap following the current base.
    // If so, add a bunch of 'N's.
    if (pos_on_ref >= region_end) {
      int64_t next_region_start = regions[reg_id+1].first;
      std::vector<uint8_t> jump(next_region_start - region_end + 1);
      alignment.insert(alignment.end(), (next_region_start - region_end - 1), EDLIB_N);

      reg_id += 1;
      region_start = regions[reg_id].first - 1;     // Coordinates are 1-based.
      region_end = regions[reg_id].second - 1;
    }

    if (op == EDLIB_M || op == EDLIB_EQUAL || op == EDLIB_X) {
      pos_on_ref += 1;
      pos_on_read += 1;
    } else if (op == EDLIB_D) {
      pos_on_ref += 1;
    } else if (op == EDLIB_I || op == EDLIB_S) {
      pos_on_read += 1;
    }
  }

  // Copy the new alignment.
  aln->alignment = alignment;
//  aln->ref_header = index->
//  aln->ref_start =
      // Trebam implementirati trans_id_to_genome_id, pa onda genome_headers za headere iz originalnog FASTA file-a, i onda izracunati poziciju za alignment.

  return 0;
}
