#include "aligner_ksw2.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
// #include "ksw2/kseq.h"
#include "aligner_util.hpp"

// KSEQ_INIT(gzFile, gzread)

namespace is {

uint8_t seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

std::shared_ptr<AlignerBase> createAlignerKSW2(const is::PiecewisePenalties &p, const is::AlignmentOptions &opt) {
  return std::shared_ptr<AlignerBase>(new AlignerKSW2(p, opt));
}

static void print_aln(const char *tname, const char *qname, ksw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez->score);
	printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
	if (ez->n_cigar > 0) {
		int i;
		putchar('\t');
		for (i = 0; i < ez->n_cigar; ++i)
			printf("%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
	}
	putchar('\n');
}


AlignerKSW2::AlignerKSW2(const is::PiecewisePenalties &p, const is::AlignmentOptions &opt) : p_(p), opt_(opt), result_(nullptr) {

}

AlignerKSW2::~AlignerKSW2() {

}

is::AlignmentReturnValue AlignerKSW2::Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) {
	void *km = 0;
	ksw_extz_t ez;  // Alignment result.
	int w = -1, flag = 0, zdrop = -1;

  #ifdef HAVE_KALLOC
    km = km_init();
  #endif

	memset(&ez, 0, sizeof(ksw_extz_t));

  auto mat = GenerateSimpleMatchMatrix<int8_t>((int8_t) p_.match, (int8_t) p_.mismatch, 5);
  // In GraphMap definition, penalties are negative. KSW2 expects positive values.
  int8_t q = -p_.w[0].v;  // Gap open. The intercept component of the affine function.
  int8_t e = -p_.w[0].u;  // Gap extend. The slope of the affine function.
  int8_t q2 = -p_.w[1].v;
  int8_t e2 = -p_.w[1].u;

  KSW2GlobalAlnWrapper_(km, (const int8_t*) qseq, qlen, (const int8_t*) tseq, tlen, 5, &mat[0], q, e, q2, e2, w, zdrop, flag, &ez);

  // print_aln("Query", "Target", &ez);

  result_ = std::shared_ptr<is::AlignmentResult>(new is::AlignmentResult);
  result_->score = ez.score;
  result_->position = is::AlignmentPosition(0, qlen, 0, tlen);
  result_->k = -1;
  result_->rv = is::AlignmentReturnValue::OK;

  result_->cigar.clear();
  std::vector<is::CigarOp> basic_cigar;
  for (size_t i=0; i<ez.n_cigar; i++) {
    basic_cigar.push_back(is::CigarOp("MID"[ez.cigar[i]&0xf], ez.cigar[i]>>4));
  }
  result_->cigar = is::ConvertBasicToExtCIGAR(qseq, qlen, tseq, tlen, basic_cigar);

  result_->edit_dist = EditDistFromExtCIGAR(result_->cigar);

  // printf ("Converted CIGAR:\n");
  // for (size_t i=0; i<result_->cigar.size(); i++) {
  //   printf ("%d%c", result_->cigar[i].count, result_->cigar[i].op);
  // }
  // printf ("\n");
  // printf ("Edit distance: %ld\n", result_->edit_dist);

  kfree(km, ez.cigar);
  #ifdef HAVE_KALLOC
    km_destroy(km);
  #endif

  return is::AlignmentReturnValue::OK;
}

is::AlignmentReturnValue AlignerKSW2::Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen, int32_t bandwidth, int32_t zdrop) {
  result_ = std::shared_ptr<is::AlignmentResult>(new is::AlignmentResult);

  if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
    return is::AlignmentReturnValue::InvalidOptions;
  }

	void *km = 0;
	ksw_extz_t ez;  // Alignment result.
	int flag = KSW_EZ_SCORE_ONLY | KSW_EZ_EXTZ_ONLY;

  #ifdef HAVE_KALLOC
    km = km_init();
  #endif

	memset(&ez, 0, sizeof(ksw_extz_t));

  auto mat = GenerateSimpleMatchMatrix<int8_t>((int8_t) p_.match, (int8_t) p_.mismatch, 5);
  // In GraphMap definition, penalties are negative. KSW2 expects positive values for affine pieces.
  int8_t q = -p_.w[0].v;  // Gap open. The intercept component of the affine function.
  int8_t e = -p_.w[0].u;  // Gap extend. The slope of the affine function.
  int8_t q2 = -p_.w[1].v;
  int8_t e2 = -p_.w[1].u;

  KSW2GlobalAlnWrapper_(km, (const int8_t*) qseq, qlen, (const int8_t*) tseq, tlen, 5, &mat[0], q, e, q2, e2, bandwidth, zdrop, flag, &ez);

  // print_aln("Query", "Target", &ez);

  result_->score = ez.score;
  result_->position = is::AlignmentPosition(0, qlen, 0, tlen);
  result_->k = -1;
  result_->rv = is::AlignmentReturnValue::OK;
  result_->max_score = ez.max;
  result_->max_q_pos = ez.max_q;
  result_->max_t_pos = ez.max_t;

  result_->cigar.clear();
  std::vector<is::CigarOp> basic_cigar;
  for (size_t i=0; i<ez.n_cigar; i++) {
    basic_cigar.push_back(is::CigarOp("MID"[ez.cigar[i]&0xf], ez.cigar[i]>>4));
  }
  result_->cigar = is::ConvertBasicToExtCIGAR(qseq, qlen, tseq, tlen, basic_cigar);

  result_->edit_dist = EditDistFromExtCIGAR(result_->cigar);

  // printf ("Converted CIGAR:\n");
  // for (size_t i=0; i<result_->cigar.size(); i++) {
  //   printf ("%d%c", result_->cigar[i].count, result_->cigar[i].op);
  // }
  // printf ("\n");
  // printf ("Edit distance: %ld\n", result_->edit_dist);

  kfree(km, ez.cigar);
  #ifdef HAVE_KALLOC
    km_destroy(km);
  #endif

  return is::AlignmentReturnValue::OK;
}

is::AlignmentReturnValue AlignerKSW2::Local(const char* q, int64_t qlen, const char* t, int64_t tlen) {
  return is::AlignmentReturnValue::NotImplementedYet;
}

is::AlignmentReturnValue AlignerKSW2::Semiglobal(const char* q, int64_t qlen, const char* t, int64_t tlen) {
  return is::AlignmentReturnValue::NotImplementedYet;
}

std::shared_ptr<is::AlignmentResult> AlignerKSW2::getResults() {
  return result_;
}

void AlignerKSW2::KSW2GlobalAlnWrapper_(void *km,
                       const int8_t *qseq_, int qlen, const int8_t *tseq_, int tlen,
                       int8_t m, const int8_t *mat,
                       int8_t q, int8_t e, int8_t q2, int8_t e2,
                       int w, int zdrop, int flag, ksw_extz_t *ez) {
	int i;
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
	ez->n_cigar = 0;

  auto qseq = ConvertSeqAlphabet(qseq_, qlen, &seq_nt4_table[0]);
  auto tseq = ConvertSeqAlphabet(tseq_, tlen, &seq_nt4_table[0]);

	ksw_extd2_sse(km, qlen, (const uint8_t*) &qseq[0],
                tlen, (const uint8_t*) &tseq[0],
                m, mat, q, e, q2, e2, w, zdrop, flag, ez);

  // const char *algo = "extd2_sse";
	// if (strcmp(algo, "extz2_sse") == 0)   ksw_extz2_sse(km, qlen, (const uint8_t*)&qseq[0], tlen, (const uint8_t*)&tseq[0], m, mat, q, e, w, zdrop, flag, ez);
	// else if (strcmp(algo, "extd2_sse") == 0)   ksw_extd2_sse(km, qlen, (const uint8_t*)&qseq[0], tlen, (const uint8_t*)&tseq[0], m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	// // else if (strcmp(algo, "extf2_sse") == 0)   ksw_extf2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, mat[0], mat[1], e, w, zdrop, ez);
	// else {
	// 	fprintf(stderr, "ERROR: can't find algorithm '%s'\n", algo);
	// 	exit(1);
	// }
}



}
