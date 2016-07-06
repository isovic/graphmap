/*
 * process_read.cc
 *
 *  Created on: Jun 6, 2016
 *      Author: isovic
 */

#include "owler2/owler2.h"

int Owler2::ProcessRead(const ProgramParameters& parameters, const IndexGappedMinimizer* index, const SingleSequence* read, const std::vector<CompiledShape> &lookup_shapes, OwlerResult *owler_result) {
  const auto &ls = index->get_lookup_shapes();
  std::vector<uint128_t> minimizers;
  minimizers.resize(ls.size() * read->get_sequence_length() * 2);     // Allocate maximum space for the seeds. Factor 2 is for the reverse complement.
  int64_t num_minimizers = 0;
  IndexGappedMinimizer::CollectMinimizers(read->get_data(), read->get_quality(), read->get_sequence_length(), -1, 0, ls, 5, minimizers, &num_minimizers);

  for (int64_t i=0; i<num_minimizers; i++) {
  }

  return 0;
}
