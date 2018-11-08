/*
 * aligner_base.h
 *
 *  Created on: Jan 7, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNER_BASE_H_
#define SRC_ALIGNER_ALIGNER_BASE_H_

#include <memory>
#include <vector>
#include "aligner_containers.h"
#include "pairwise_penalties.h"

namespace is {

class AlignerBase {
 public:
  virtual ~AlignerBase() { }

  // virtual AlignmentReturnValue Align(const char* q, int64_t qlen, const char* t, int64_t tlen, AlignmentType type) = 0;    // Selects the alignment mode based on a parameter.

  virtual AlignmentReturnValue Global(const char* q, int64_t qlen, const char* t, int64_t tlen, bool type) = 0;      // Global alignment mode.

  virtual AlignmentReturnValue Local(const char* q, int64_t qlen, const char* t, int64_t tlen) = 0;       // Local alignment mode.

  virtual AlignmentReturnValue Semiglobal(const char* q, int64_t qlen, const char* t, int64_t tlen) = 0;  // Semiglobal alignment mode.

  virtual AlignmentReturnValue Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen,     // Extend alignment mode. Does not necessarily
                                      int32_t bandwidth, int32_t zdrop) = 0;                              //  produce CIGAR,but generate max alignment coords

  virtual std::shared_ptr<AlignmentResult> getResults() = 0;

};

} /* namespace is */

#endif /* SRC_ALIGNER_ALIGNER_BASE_H_ */
