/*
 * aligner_containers.h
 *
 *  Created on: Jan 7, 2017
 *      Author: isovic
 */

#ifndef SRC_CONTAINERS_H_
#define SRC_CONTAINERS_H_

#include <stdint.h>
#include <limits>
#include <string>
#include <sstream>

#include "sam_parser.h"

namespace is {

enum class AlignmentReturnValue { // Strongly typed enum, C++11 feature.
  OK,                   // Everything went ok.
  Suboptimal,            // Alignment stepped out of defined band. Result is not optimal.
  InvalidOptions,       // In case parameters of values are invalid.
  QlenIsZero,
  TlenIsZero,
  WrongEditDist,
  AlignmentNotPerformed,  // A default value for an alignment which wasn't performed.
  NotImplementedYet     // For features in development.
};

enum class AlignmentType {  // Strongly typed enum, C++11 feature.
  Global,
  Local
};

class AlignmentPosition {
 public:
  AlignmentPosition() : qstart(0), qend(0), tstart(0), tend(0) { }
  AlignmentPosition(int32_t _qstart, int32_t _qend, int32_t _tstart, int32_t _tend) :
                    qstart(_qstart), qend(_qend), tstart(_tstart), tend(_tend) { }
  AlignmentPosition(const AlignmentPosition& op) :
                    AlignmentPosition(op.qstart, op.qend, op.tstart, op.tend) { }
  AlignmentPosition& operator=(const AlignmentPosition& op) {
    qstart = op.qstart;
    qend = op.qend;
    tstart = op.tstart;
    tend = op.tend;
    return *this;
  }

  int32_t qstart, qend;              // Query and target alignment start and end positions. End position
  int32_t tstart, tend;              // is inclusive (the position of the last base).
};

class AlignmentResult {
 public:
  AlignmentResult() : score(0), edit_dist(0), position(), k(-1), rv(AlignmentReturnValue::AlignmentNotPerformed) {
  }

  AlignmentResult(const AlignmentResult& op) :
    score(op.score), edit_dist(op.edit_dist),
    position(op.position), cigar(op.cigar), k(op.k), rv(op.rv) { // Copy constructor.
  }

  ~AlignmentResult() { };

  AlignmentResult& operator=(const AlignmentResult& op) {
    score = op.score;
    edit_dist = op.edit_dist;
    position = op.position;
    cigar = op.cigar;
    k = op.k;
    rv = op.rv;
    return *this;
  }

  // Alignment results.
  int64_t score;
  int64_t edit_dist;
  is::AlignmentPosition position;                   // There can be multiple alignments with the same score.
                                                    // Only the first position and the corresponding alignment
  std::vector<is::CigarOp> cigar;                   // are reported
  int64_t k;                                        // Value of band k used in the final alignment.
  AlignmentReturnValue rv;                          // Return value of the aligner.
};

// If any global margin is true, then the corresponding will be penalized.
// Concretely, if top/left are true, then the first row/column will be initialized
// to the multiple of the gap extend penalty in global alignment.
// If bottom is false, the maximum of last row will be found instead of taking
// the bottom right corner for global alignment.
// If right is false, the maximum of last column will be found instead of taking
// the bottom right corner for global alignment.
class GlobalMargins {
 public:
  GlobalMargins()
      : top(true),
        left(true),
        bottom(true),
        right(true) {
  }
  GlobalMargins(bool _top, bool _left, bool _bottom, bool _right)
      : top(_top),
        left(_left),
        bottom(_bottom),
        right(_right) {
  }
  bool top, left, bottom, right;
};

class AlignmentOptions {
 public:
  AlignmentOptions() :  k(-1),
                        do_traceback(true) {
  }

  int32_t k;                // Band for banded alignment. If < 0, banded alignment is turned off.
  bool do_traceback;        // If traceback is not needed, then there is no need to alocate a large
                            // matrix to store directions.
  GlobalMargins gm;
};

} /* namespace is */



#endif /* SRC_CONTAINERS_H_ */
