/*
 * aligner_base.h
 *
 *  Created on: Jan 7, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNER_KSW2_H_
#define SRC_ALIGNER_ALIGNER_KSW2_H_

#include <memory>
#include <vector>
#include "aligner_base.h"
#include "aligner_containers.h"
#include "pairwise_penalties.h"
#include "aligner_util.hpp"
#include "ksw2/ksw2.h"

namespace is {

class AlignerKSW2;

std::shared_ptr<AlignerBase> createAlignerKSW2(const is::PiecewisePenalties &p, const is::AlignmentOptions &opt);

class AlignerKSW2 : public AlignerBase {
 public:
  friend std::shared_ptr<AlignerBase> createAlignerKSW2(const is::PiecewisePenalties &p, const is::AlignmentOptions &opt);

  ~AlignerKSW2();

  AlignmentReturnValue Global(const char* q, int64_t qlen, const char* t, int64_t tlen);   // Global alignment mode.

  AlignmentReturnValue Local(const char* q, int64_t qlen, const char* t, int64_t tlen);    // Local alignment mode.

  AlignmentReturnValue Semiglobal(const char* q, int64_t qlen, const char* t, int64_t tlen);   // Semiglobal alignment mode.

  std::shared_ptr<AlignmentResult> getResults();

 protected:
  AlignerKSW2(const is::PiecewisePenalties &p, const is::AlignmentOptions &opt);   // We don't want users attempting to instantiate manually, even though the class is virtual.

 private:
  AlignerKSW2(const AlignerKSW2&) = delete;                       // No copying.
  AlignerKSW2& operator=(const AlignerKSW2&) = delete;            // No copying.
  AlignerKSW2(AlignerKSW2&&) = delete;                            // No move constructor.
  AlignerKSW2& operator=(const AlignerKSW2&&) = delete;           // No copying.

  void KSW2GlobalAlnWrapper_(void *km,
                        const int8_t *qseq_, int qlen, const int8_t *tseq_, int tlen,
                        int8_t m, const int8_t *mat,
                        int8_t q, int8_t e, int8_t q2, int8_t e2,
                        int w, int zdrop, int flag, ksw_extz_t *ez);

  const is::PiecewisePenalties& p_;
  const is::AlignmentOptions& opt_;
  std::shared_ptr<is::AlignmentResult> result_;
};

} /* namespace is */

#endif /* SRC_ALIGNER_ALIGNER_BASE_H_ */