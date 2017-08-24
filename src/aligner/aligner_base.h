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

  virtual AlignmentReturnValue Global(const char* q, int64_t qlen, const char* t, int64_t tlen) = 0;   // Global alignment mode.

  virtual AlignmentReturnValue Local(const char* q, int64_t qlen, const char* t, int64_t tlen) = 0;    // Local alignment mode.

  virtual AlignmentReturnValue Semiglobal(const char* q, int64_t qlen, const char* t, int64_t tlen) = 0;   // Semiglobal alignment mode.

  virtual std::shared_ptr<AlignmentResult> getResults() = 0;

//  protected:
//   virtual AlignerBase(const is::PiecewisePenalties &p, const is::AlignmentOptions &opt) = 0;   // We don't want users attempting to instantiate manually, even though the class is virtual.

//  private:
//   AlignerBase(const AlignerBase&) = delete;                       // No copying.
//   AlignerBase& operator=(const AlignerBase&) = delete;            // No copying.
//   AlignerBase(AlignerBase&&) = delete;                            // No move constructor.
//   AlignerBase& operator=(const AlignerBase&&) = delete;           // No copying.

};

} /* namespace is */

#endif /* SRC_ALIGNER_ALIGNER_BASE_H_ */
