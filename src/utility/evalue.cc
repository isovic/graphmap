/*
 * evalue.cc
 *
 *  Created on: Feb 23, 2015
 *      Author: isovic
 */

#include <string.h>
#include <math.h>
#include "utility/evalue.h"
#include "log_system/log_system.h"

/* Short usage instructions:
 * Scorer *scorer;
 * ScorerCreateMatrix(&scorer, "EDNA_FULL", gap_open, gap_extend);
 *
 */

#define SCORERS_LEN (sizeof(scorers) / sizeof(ScorerEntry))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void ScorerCreateMatrix(Scorer** scorer, char* name, int gapOpen, int gapExtend) {

  int index = -1;

  // scorers is an array consisting of the scoring matrix name (i.e. EDNA_FULL, BLOSUM62, ...) and the actual values of the
  // scoring matrix. Matrices are stored as 26x26 arrays of ints.
  // This for loop searches for the matrix of a correct name, and returns its index.
  int i;
  for (i = 0; i < SCORERS_LEN; ++i) {
    if (strcmp(name, scorers[i].name) == 0) {
      index = i;
      break;
    }
  }

//    ASSERT(index != -1, "unknown table %s", name);

  // Check if the scoring matrix name was found.
  if (index >= 0) {
    // Initialize a new set of values. ScorerEntry basically only consists of the matrix name, match/mismatch scores
    // given by the scorers[] data (that is, the scoring matrix), the gap open and gap extend penalties, the size of the
    // scoring matrix maxCode (the scoring matrix is of size maxCode * maxCode), and two additional values: maxScore which
    // is of value of the maximum match/mismatch score in the scoring matrix; and 'scalar' which is 1 if the scoring
    // matrix contains the same values everywhere (diagonal has one value, and the rest of the matrix is exactly the same),
    // and 0 if the matrix is not of equal values.

    ScorerEntry* entry = &(scorers[index]);
    *scorer = ScorerCreate(entry->name, *(entry->table), 26, gapOpen, gapExtend);

//    ErrorReporting::GetInstance().VerboseLog(VERBOSE_LEVEL_ALL, true, FormatString(), "ScorerCreateMatrix");

  } else {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Could not find the correct scoring matrix!"));

  }
}

Scorer* ScorerCreate(const char* name, int* scores, char maxCode, int gapOpen, int gapExtend) {

//    ASSERT(maxCode > 0, "scorer table must have at least one element");
//    ASSERT(gapOpen > 0, "gap open is defined as positive integer");
//    ASSERT(gapExtend > 0, "gap extend is defined as positive integer");
//    ASSERT(gapOpen >= gapExtend, "gap extend must be equal or less to gap open");

//    int i;
//    int j;
//    for (i = 0; i < maxCode; ++i) {
//        for (j = i + 1; j < maxCode; ++j) {
//            int a = scores[i * maxCode + j];
//            int b = scores[j * maxCode + i];
////            ASSERT(a == b, "scorer table must be symmetrical");
//        }
//    }

  Scorer* scorer = (Scorer*) malloc(sizeof(struct Scorer));

  scorer->nameLen = strlen(name) + 1;
  scorer->name = (char*) malloc(scorer->nameLen * sizeof(char));
  scorer->name[scorer->nameLen - 1] = '\0';
  memcpy(scorer->name, name, (scorer->nameLen - 1) * sizeof(char));

  scorer->gapOpen = gapOpen;
  scorer->gapExtend = gapExtend;

  size_t tableSize = maxCode * maxCode * sizeof(int);
  scorer->table = (int*) malloc(tableSize);
  memcpy(scorer->table, scores, tableSize);

  scorer->maxCode = maxCode;
  scorer->maxScore = MaxScore(scorer);
  scorer->scalar = IsScalar(scorer);

  return scorer;
}

void ScorerDelete(Scorer* scorer) {
  if (scorer == NULL)
    return;

  if (scorer->name)
    free(scorer->name);
  if (scorer->table)
    free(scorer->table);
  free(scorer);
}

int IsScalar(Scorer* scorer) {

  int x, i, j;

  int scorerLen = scorer->maxCode;

  x = scorer->table[0];
  for (i = 1; i < scorerLen; ++i) {
    if (scorer->table[i * scorerLen + i] != x) {
      return 0;
    }
  }

  x = scorer->table[1];
  for (i = 0; i < scorerLen; ++i) {
    for (j = 0; j < scorerLen; ++j) {
      if (i != j && scorer->table[i * scorerLen + j] != x) {
        return 0;
      }
    }
  }

  return 1;
}

int MaxScore(Scorer* scorer) {

  int i, j;

  int scorerLen = scorer->maxCode;

  int max = scorer->table[0];
  for (i = 0; i < scorerLen; ++i) {
    for (j = 0; j < scorerLen; ++j) {
      max = MAX(max, scorer->table[i * scorerLen + j]);
    }
  }

  return max;
}

int SetupScorer(char* matrix_name, int64_t reference_length, int gap_open, int gap_extend, EValueParams **ret_eValueParams) {
//  char* matrix_name = "EDNA_FULL";
//  int gapOpen = 10;
//  int gapExtend = 1;

  Scorer* scorer = NULL;
  ScorerCreateMatrix(&scorer, matrix_name, gap_open, gap_extend);

//  long long cells;
//  statFastaChains(&chains, &cells, databasePath);
//  long long reference_length = 1;
  if (scorer == NULL) {
    LogSystem::GetInstance().Log(SEVERITY_INT_FATAL, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Scorer not initialized!"));
    return 1;
  }

  *ret_eValueParams = CreateEValueParams(reference_length, scorer);
  ScorerDelete(scorer);

  return 0;
}

EValueParams* CreateEValueParams(long long length, Scorer* scorer) {

  const char* matrix = scorer->name;
  int gapOpen = scorer->gapOpen;
  int gapExtend = scorer->gapExtend;

  int index = -1;
  int indexUn = -1;

  for (int i = 0; i < (int) SCORER_CONSTANTS_LEN; ++i) {

    ScorerConstants* entry = &(scorerConstants[i]);

    int sameName = strcmp(entry->matrix, matrix) == 0;

    if (sameName && indexUn == -1) {
      indexUn = i;
    }

    if (sameName && entry->gapOpen == gapOpen && entry->gapExtend == gapExtend) {
      index = i;
      break;
    }
  }

  // ignore miss of default matrix parameters
  if (indexUn == -1 && index != -1) {
    indexUn = index;
  }

  if (index == -1 || indexUn == -1) {

    if (indexUn == -1) {
      index = 0;
      indexUn = 0;
    } else {
      index = indexUn;
    }

    ScorerConstants* entry = &(scorerConstants[indexUn]);

//        WARNING(1, "no e-value params found, using %s %d %d", entry->matrix,
//            entry->gapOpen, entry->gapExtend);
  }

  double alphaUn = scorerConstants[indexUn].alpha;
  double aUn = scorerConstants[indexUn].a;
  double G = gapOpen + gapExtend;

  EValueParams* params = (EValueParams*) malloc(sizeof(struct EValueParams));

  params->G = G;
  params->aUn = aUn;
  params->alphaUn = alphaUn;
  params->lambda = scorerConstants[index].lambda;
  params->K = scorerConstants[index].K;
  params->logK = log(params->K);
  params->H = scorerConstants[index].H;
  params->a = scorerConstants[index].a;
  params->C = scorerConstants[index].C;
  params->alpha = scorerConstants[index].alpha;
  params->sigma = scorerConstants[index].sigma;
  params->b = 2.0 * G * (params->aUn - params->a);
  params->beta = 2.0 * G * (params->alphaUn - params->alpha);
  params->tau = 2.0 * G * (params->alphaUn - params->sigma);
  params->length = length;
  params->isDna = scorerConstants[index].isDna;

//    printf("Using: lambda = %.3lf, K = %.3lf, H = %.3lf\n",
//        params->lambda, params->K, params->H);

  return params;
}

void DeleteEValueParams(EValueParams* eValueParams) {
  free(eValueParams);
  eValueParams = NULL;
}

//static void eValuesCpu(double* values, int* scores, Chain* query,
//    Chain** database, int databaseLen, EValueParams* eValueParams) {
//
//    double (*function) (int, int, int, EValueParams*);
//
//    if (eValueParams->isDna) {
//        function = calculateEValueDna;
//    } else {
//        function = calculateEValueProt;
//    }
//
//    int queryLen = chainGetLength(query);
//
//    for (int i = 0; i < databaseLen; ++i) {
//
//        int score = scores[i];
//        int targetLen = chainGetLength(database[i]);
//
//        if (score == NO_SCORE) {
//            values[i] = INFINITY;
//            continue;
//        }
//
//        values[i] = function(score, queryLen, targetLen, eValueParams);
//    }
//}

int CalculateEValueDNA(int64_t alignment_score, int64_t query_length, int64_t target_length, const EValueParams* eValueParams, double *ret_evalue) {

  double lambda = eValueParams->lambda;
  double logK = eValueParams->logK;

  *ret_evalue = (double) query_length * target_length * exp(-lambda * alignment_score + logK);

  return 0;
}

