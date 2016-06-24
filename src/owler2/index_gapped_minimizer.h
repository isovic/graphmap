/*
 * index_gapped_minimizer.h
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#ifndef SRC_OWLER2_INDEX_GAPPED_MINIMIZER_H_
#define SRC_OWLER2_INDEX_GAPPED_MINIMIZER_H_

#include "sparsehash/dense_hash_map"

using google::dense_hash_map;      // namespace where class lives by default
//typedef dense_hash_map<int64_t, DenseType2, std::hash<int64_t> > DenseType;
//std::vector<DenseType> bins_map1;
//bins_map1[i].set_empty_key(-1);
//struct DenseType2 {
//  int32_t timestamp = 0;
//  float count = 0.0;
//};
//DenseType2 &hit = temp_map[position_bin];
//if (hit.timestamp == (i + 1)) { continue; }
//hit.count += 1.0f;
//hit.timestamp = (i + 1);

#include <vector>
#include <map>
#include <stdint.h>
#include <stdint.h>
#include "sequences/sequence_file.h"
#include "utility/utility_general.h"
#include "compiled_shape.h"

struct IndexPos {
  int32_t id;     // ID of the originating sequence.
  uint32_t pos;   // Position within the sequence. If >= (1 << 31), the sequence is reverse complemented.
};

class IndexGappedMinimizer {
 public:
  IndexGappedMinimizer();
  ~IndexGappedMinimizer();

  int AddSequence(const SingleSequence* seq, const std::vector<CompiledShape> &compiled_shapes, float min_avg_seed_qv, bool index_reverse_strand, int32_t num_threads);
  int CreateFromSequenceFile(const SequenceFile &seqs, const std::vector<CompiledShape> &compiled_shapes, float min_avg_seed_qv, bool index_reverse_strand, int32_t num_threads);
};

#endif /* SRC_OWLER2_INDEX_GAPPED_MINIMIZER_H_ */
