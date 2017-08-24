#include "aligner_ksw2.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "ksw2/ksw2.h"
#include "ksw2/kseq.h"

KSEQ_INIT(gzFile, gzread)

namespace is {

std::shared_ptr<AlignerBase> createAlignerKSW2(const is::PiecewisePenalties &p, const is::AlignmentOptions &opt) {
  return std::shared_ptr<AlignerBase>(new AlignerKSW2(p, opt));
}

unsigned char seq_nt4_table[256] = {
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

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

// static void global_aln(const char *algo, void *km,
//                        const char *qseq_, int qlen, const char *tseq_, int tlen,
//                        int8_t m, const int8_t *mat,
//                        int8_t q, int8_t e, int8_t q2, int8_t e2,
//                        int w, int zdrop, int flag, ksw_extz_t *ez)
// {
// 	int i;
// 	uint8_t *qseq, *tseq;
// 	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
// 	ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
// 	ez->n_cigar = 0;
// 	qseq = (uint8_t*)calloc(qlen + 33, 1); // 32 for gaba
// 	tseq = (uint8_t*)calloc(tlen + 33, 1);
// 	for (i = 0; i < qlen; ++i)
// 		qseq[i] = seq_nt4_table[(uint8_t)qseq_[i]];
// 	for (i = 0; i < tlen; ++i)
// 		tseq[i] = seq_nt4_table[(uint8_t)tseq_[i]];
// 	if (strcmp(algo, "gg") == 0) {
// 		if (flag & KSW_EZ_SCORE_ONLY) ez->score = ksw_gg(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, 0, 0, 0);
// 		else ez->score = ksw_gg(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
// 	} else if (strcmp(algo, "gg2") == 0) {
// 		if (flag & KSW_EZ_SCORE_ONLY) ez->score = ksw_gg2(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, 0, 0, 0);
// 		else ez->score = ksw_gg2(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
// 	}
// 	else if (strcmp(algo, "gg2_sse") == 0)     ez->score = ksw_gg2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
// 	else if (strcmp(algo, "extz") == 0)        ksw_extz(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, zdrop, flag, ez);
// 	else if (strcmp(algo, "extz2_sse") == 0)   ksw_extz2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, zdrop, flag, ez);
// 	else if (strcmp(algo, "extd") == 0)        ksw_extd(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
// 	else if (strcmp(algo, "extd2_sse") == 0)   ksw_extd2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
// 	else if (strcmp(algo, "extf2_sse") == 0)   ksw_extf2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, mat[0], mat[1], e, w, zdrop, ez);
// 	else if (strcmp(algo, "test") == 0) ksw_extd2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, 4, 2, 24, 1, 751, 400, 8, ez);
// 	else {
// 		fprintf(stderr, "ERROR: can't find algorithm '%s'\n", algo);
// 		exit(1);
// 	}
// 	free(qseq); free(tseq);
// }

static void global_aln(const char *algo, void *km,
                       const char *qseq_, int qlen, const char *tseq_, int tlen,
                       int8_t m, const int8_t *mat,
                       int8_t q, int8_t e, int8_t q2, int8_t e2,
                       int w, int zdrop, int flag, ksw_extz_t *ez)
{
	int i;
	uint8_t *qseq, *tseq;
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
	ez->n_cigar = 0;
	qseq = (uint8_t*)calloc(qlen + 33, 1); // 32 for gaba
	tseq = (uint8_t*)calloc(tlen + 33, 1);

	for (i = 0; i < qlen; ++i)
		qseq[i] = seq_nt4_table[(uint8_t)qseq_[i]];
	for (i = 0; i < tlen; ++i)
		tseq[i] = seq_nt4_table[(uint8_t)tseq_[i]];

	if (strcmp(algo, "extz2_sse") == 0)   ksw_extz2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, w, zdrop, flag, ez);
	else if (strcmp(algo, "extd2_sse") == 0)   ksw_extd2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, m, mat, q, e, q2, e2, w, zdrop, flag, ez);
	// else if (strcmp(algo, "extf2_sse") == 0)   ksw_extf2_sse(km, qlen, (uint8_t*)qseq, tlen, (uint8_t*)tseq, mat[0], mat[1], e, w, zdrop, ez);
	else {
		fprintf(stderr, "ERROR: can't find algorithm '%s'\n", algo);
		exit(1);
	}

	free(qseq); free(tseq);
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
	int8_t a = 2; // Match.
  int8_t b = 4; // Mismatch.
  int8_t q = 4, e = 2, q2 = 13, e2 = 1; // Piecewise gap penalties.
	std::string algo = "extd2_sse";
  char *s;
	int8_t mat[25]; // Match/mismatch matrix.
	ksw_extz_t ez;  // Alignment result.
	int w = -1, flag = 0, zdrop = -1;

  #ifdef HAVE_KALLOC
    km = km_init();
  #endif

	memset(&ez, 0, sizeof(ksw_extz_t));
	ksw_gen_simple_mat(5, mat, a, -b);

  // printf ("Query:\n");
  // for (int64_t i = 0; i < qlen; i++) {
  //   printf ("%c", qseq[i]);
  // }
  // printf ("\n");
  // printf ("\n");
  // printf ("Target:\n");
  // for (int64_t i = 0; i < tlen; i++) {
  //   printf ("%c", tseq[i]);
  // }
  // printf ("\n");
  // printf ("\n");

  global_aln(algo.c_str(), km, qseq, qlen, tseq, tlen, 5, mat, q, e, q2, e2, w, zdrop, flag, &ez);
  print_aln("Query", "Target", &ez);

  kfree(km, ez.cigar);
  #ifdef HAVE_KALLOC
    km_destroy(km);
  #endif

  return is::AlignmentReturnValue::NotImplementedYet;
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

}
