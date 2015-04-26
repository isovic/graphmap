/*
 * evalue.h
 *
 *  Created on: Feb 23, 2015
 *      Author: isovic
 */

#ifndef EVALUE_H_
#define EVALUE_H_

#include <stdio.h>
#include <stdlib.h>
#include "utility/evalue_constants.h"



struct EValueParams {
    double lambda;
    double K;
    double logK;
    double H;
    double a;
    double C;
    double alpha;
    double sigma;
    double b;
    double beta;
    double tau;
    double G;
    double aUn;
    double alphaUn;
    long long length;
    int isDna;
};

struct Scorer {

    char* name;
    int nameLen;

    int gapOpen;
    int gapExtend;

    int* table;

    int maxCode;
    int maxScore;
    int scalar;
};

typedef struct ScorerConstants {
    const char* matrix;
    int gapOpen;
    int gapExtend;
    double lambda;
    double K;
    double H;
    double a;
    double C;
    double alpha;
    double sigma;
    int isDna;
} ScorerConstants;

static ScorerConstants scorerConstants[] = {
//     name         gapO   gapE lambda  K      H       a       C         alpha     sigma     isDna
//    { "BLOSUM_62",  -1,     -1, 0.3176, 0.134, 0.4012, 0.7916, 0.623757, 4.964660, 4.964660, 0},
//    { "BLOSUM_62",  11,     2, 0.297, 0.082, 0.27, 1.1, 0.641766, 12.673800, 12.757600, 0 },
//    { "BLOSUM_62",  10,     2, 0.291, 0.075, 0.23, 1.3, 0.649362, 16.474000, 16.602600, 0 },
//    { "BLOSUM_62",  9,      2, 0.279, 0.058, 0.19, 1.5, 0.659245, 22.751900, 22.950000, 0 },
//    { "BLOSUM_62",  8,      2, 0.264, 0.045, 0.15, 1.8, 0.672692, 35.483800, 35.821300, 0 },
//    { "BLOSUM_62",  7, 2, 0.239, 0.027, 0.10, 2.5, 0.702056, 61.238300, 61.886000, 0 },
//    { "BLOSUM_62",  6, 2, 0.201, 0.012, 0.061, 3.3, 0.740802, 140.417000, 141.882000, 0 },
//    { "BLOSUM_62",  13, 1, 0.292, 0.071, 0.23, 1.2, 0.647715, 19.506300, 19.893100, 0 },
//    { "BLOSUM_62",  12, 1, 0.283, 0.059, 0.19, 1.5, 0.656391, 27.856200, 28.469900, 0 },
//    { "BLOSUM_62",  11, 1, 0.267, 0.041, 0.14, 1.9, 0.669720, 42.602800, 43.636200, 0 },
//    { "BLOSUM_62",  10, 1, 0.243, 0.024, 0.10, 2.5, 0.693267, 83.178700, 85.065600, 0 },
//    { "BLOSUM_62",  9, 1, 0.206, 0.010, 0.052, 4.0, 0.731887, 210.333000, 214.842000, 0 },
    { "EDNA_FULL_1_1", 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99 },
    { "EDNA_FULL_1_1", 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97 },
    { "EDNA_FULL_1_1", 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92 },
    { "EDNA_FULL_1_1", 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72 },
    { "EDNA_FULL_1_1", 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98 },
    { "EDNA_FULL_1_1", 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96 },
    { "EDNA_FULL_1_1", 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90 },
    { "EDNA_FULL_5_4", 10, 6, 0.163, 0.068, 0.16, 0, 0, 0, 0, 1 },
    { "EDNA_FULL_5_4", 8, 6, 0.146, 0.039, 0.11, 0, 0, 0, 0, 1 },
    { "EDNA_FULL_4_5", 6, 5, 0.28,  0.21, 0.47, 0.6 , -7, 93 },
    { "EDNA_FULL_4_5", 0, 0, 0.22, 0.061, 0.22, 1.0, -15, 74 },
    { "EDNA_FULL_4_5", 3, 5, 0.23, 0.065, 0.25, 0.9, -11, 76 }
};

///** Number of statistical parameters in each row of the precomputed tables. */
//#define BLAST_NUM_STAT_VALUES 11  /**< originally 8, now 11 to support Spouge's FSC. see notes below */
//
///** Holds values (gap-opening, extension, etc.) for a matrix. */
//typedef double array_of_8[BLAST_NUM_STAT_VALUES];
//
///** Karlin-Altschul parameter values for substitution scores 1 and -1. */
//static const array_of_8 blastn_values_1_1[] = {
//    { 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99 },
//    { 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97 },
//    { 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92 },
//    { 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72 },
//    { 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98 },
//    { 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96 },
//    { 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90 }
//};

typedef struct ScorerEntry {
    const char* name;
    int (*table)[26 * 26];
} ScorerEntry;

// to register a scorer just add his name and corresponding table to this array
static ScorerEntry scorers[] = {
//    { "BLOSUM_62", &BLOSUM_62_TABLE }, // default one
//    { "BLOSUM_45", &BLOSUM_45_TABLE },
//    { "BLOSUM_50", &BLOSUM_50_TABLE },
//    { "BLOSUM_80", &BLOSUM_80_TABLE },
//    { "BLOSUM_90", &BLOSUM_90_TABLE },
//    { "PAM_30", &PAM_30_TABLE },
//    { "PAM_70", &PAM_70_TABLE },
//    { "PAM_250", &PAM_250_TABLE },
    { "EDNA_FULL_5_4", &EDNA_FULL_TABLE_5_4 },
    { "EDNA_FULL_1_1", &EDNA_FULL_TABLE_1_1 },
    { "EDNA_FULL_4_5", &EDNA_FULL_TABLE_4_5 }
};

typedef struct EValueParams EValueParams;
typedef struct Scorer Scorer;

#define SCORER_CONSTANTS_LEN (sizeof(scorerConstants) / sizeof(ScorerConstants))

EValueParams* CreateEValueParams(long long length, Scorer* scorer);
void DeleteEValueParams(EValueParams* eValueParams);

/*!
@brief Scorer constructor.

Input scores table must be an array which length is equal to maxCode * maxCode.
Input scores table must be organized so that columns and rows correspond to the
codes shown in scorerEncode(char).

@param name scorer name, copy is made
@param scores similarity table, copy is made
@param maxCode maximum code that scorer should work with
@param gapOpen gap open penalty given as a positive integer
@param gapExtend gap extend penalty given as a positive integer

@return scorer object
*/
Scorer *ScorerCreate(const char* name, int* scores, char maxCode,
    int gapOpen, int gapExtend);
void ScorerDelete(Scorer* scorer);

int IsScalar(Scorer* scorer);

int MaxScore(Scorer* scorer);

int CalculateEValueDNA(int64_t alignment_score, int64_t query_length, int64_t target_length, const EValueParams* eValueParams, double *ret_evalue);

int SetupScorer(char* matrix_name, int64_t reference_length, int gap_open, int gap_extend, EValueParams **ret_eValueParams);





#endif /* EVALUE_H_ */
