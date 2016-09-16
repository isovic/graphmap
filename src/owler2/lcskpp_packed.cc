/*
 * lcskpp.cc
 *
 *  Created on: Nov 24, 2015
 *      Author: isovic
 */

#include "lcskpp_packed.h"

#include "qsort.h"
#include "omp_sort.hpp"

inline double dist(double x1, double y1, double x2, double y2) {
  return (sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)));
}

int LCSkPacked(const uint128_t *packed_hits, int64_t n_hits, int32_t k, int64_t *ret_lcskpp_length, std::vector<int32_t> *ret_lcskpp_indices, bool print_debug) {
  if (packed_hits == NULL || n_hits <= 0)
    return 1;

  unsigned __int128 *events = NULL;
  uint64_t *matches_starts = NULL;

  uint32_t n = 0;
  int64_t num_matches = 0;

  int64_t num_events = 0;
  int64_t lcskpp_length = 0;

  int32_t min_ref = get_lcsk128_rpos(packed_hits[0]), min_query = get_lcsk128_qpos(packed_hits[0]);

  for (uint32_t i=0; i<n_hits; i++) {
    min_ref = std::min(min_ref, get_lcsk128_rpos(packed_hits[i]));
    min_query = std::min(min_query, get_lcsk128_qpos(packed_hits[i]));
  }

  events = (unsigned __int128 *) malloc(sizeof(unsigned __int128) * n_hits * 2);
  matches_starts = (uint64_t *) malloc(sizeof(uint64_t) * n_hits);

  clock_t time_preprocessing = clock();
  for (uint32_t i=0; i<n_hits; i++) {
    uint32_t ref_start = get_lcsk128_rpos(packed_hits[i]) - min_ref;
    uint32_t ref_end = ref_start + k;
    uint32_t query_start = get_lcsk128_qpos(packed_hits[i]) - min_query;
    uint32_t query_end = query_start + k;

    unsigned __int128 event1 = (((unsigned __int128) ref_start) << (8 * 8)) | (((unsigned __int128) query_start) << (4 * 8)) | (((unsigned __int128) (i + n_hits)));
    events[num_events] = event1;
    num_events += 1;

    unsigned __int128 event2 = (((unsigned __int128) ref_end) << (8 * 8)) | (((unsigned __int128) query_end) << (4 * 8)) | ((((unsigned __int128) i)));
    events[num_events] = event2;
    num_events += 1;

    matches_starts[num_matches] = (uint64_t) (event1 >> (4 * 8));

    num_matches += 1;

    n = std::max(n, ref_end);
    n = std::max(n, query_end);
  }
  if (print_debug == true) {
    printf ("Time for preprocessing: %f\n", (((float) (clock() - time_preprocessing))/CLOCKS_PER_SEC));
    fflush(stdout);
  }

  clock_t time_sorting = clock();
  std::sort(events, (events + num_events));
//  quickSort(events, num_events);
//  #define islt_event(a,b) (a < b)
//  QSORT(unsigned __int128, events, num_events, islt_event);

  if (print_debug == true) {
    printf ("Time for sorting: %f\n", (((float) (clock() - time_sorting))/CLOCKS_PER_SEC));
    fflush(stdout);
  }

  clock_t time_pp = clock();
  // Indexed by column, first:dp value, second:index in matches.
  FenwickMax<std::pair<int, int> > dp_col_max(n);
  std::vector<int> dp(num_matches);
  std::vector<int> recon(num_matches);
  std::vector<int> continues(num_matches, -1);
  for (int64_t curr = 0; curr < num_matches; curr++) {
    uint64_t G = 0;
    uint64_t G1 = (((matches_starts[curr] >> (4 * 8)) & (0x00000000FFFFFFFF)) - 1);
    uint64_t G2 = (((matches_starts[curr]) & (0x00000000FFFFFFFF)) - 1);
    G = (G1 << (4 * 8)) | G2;

    auto prev = std::lower_bound(matches_starts, (matches_starts + num_matches), G);

    if (prev != (matches_starts + num_matches) && *prev == G) {
      continues[curr] = prev - matches_starts;
    }
  }
  if (print_debug == true) {
    printf ("Time for finding extensions %f\n", (((float) (clock() - time_pp))/CLOCKS_PER_SEC));
    fflush(stdout);
  }

  clock_t time_processing = clock();
  int best_idx = 0;
  lcskpp_length = 0;
  for (int64_t current_event=0; current_event<num_events; current_event++) {
    int64_t raw_idx = (int64_t) ((events[current_event] & (0x00000000FFFFFFFF)));

    int64_t idx = (raw_idx >= ((int64_t) num_matches)) ? (raw_idx - ((int64_t) num_matches)) : (raw_idx);
    bool is_beginning = (raw_idx >= num_matches);
    uint64_t i = (uint64_t) ((events[current_event] >> (8 * 8)) & (0x00000000FFFFFFFF));
    uint64_t j = (uint64_t) ((events[current_event] >> (4 * 8)) & (0x00000000FFFFFFFF));
    int primary_diagonal = n - 1 + i - j;

    if (is_beginning) { // begin
      std::pair<int, int> prev_dp = dp_col_max.get(j);
      dp[idx] = k;
      recon[idx] = -1;

      if (prev_dp.first > 0) {
//        uint64_t prev_i = (uint64_t) ((events[prev_dp.second] >> (8 * 8)) & (0x00000000FFFFFFFF));
//        uint64_t prev_j = (uint64_t) ((events[prev_dp.second] >> (4 * 8)) & (0x00000000FFFFFFFF));
//        int penalty = PenaltyAngular(0, 0, i, j, prev_i, prev_j);
        dp[idx] = prev_dp.first + k; // - penalty;
        recon[idx] = prev_dp.second;
      }
    } else {
      if (continues[idx] != -1) {
        if (dp[continues[idx]] + 1 > dp[idx]) {
          dp[idx] = dp[continues[idx]] + 1;
          recon[idx] = continues[idx];
        }
      }

      dp_col_max.update(j, std::make_pair(dp[idx], idx));

      if (dp[idx] > lcskpp_length) {
        lcskpp_length = dp[idx];
        best_idx = idx;
      }
    }
  }
  if (print_debug == true) {
    printf ("Time for processing: %f\n", (((float) (clock() - time_processing))/CLOCKS_PER_SEC));
    fflush(stdout);
  }

  clock_t time_reconstruction = clock();
  ret_lcskpp_indices->clear();
  if (best_idx != -1 && recon.size() > 0) {
    int64_t num_recon = 1;
    for (int i1 = best_idx; i1 != -1; i1 = recon[i1]) { if (recon[i1] != -1) { num_recon += 1; } }
    ret_lcskpp_indices->resize(num_recon);

    int64_t curr_recon = 0;
    (*ret_lcskpp_indices)[curr_recon++] = best_idx;
    for (int i1 = best_idx; i1 != -1; i1 = recon[i1]) {
      if (recon[i1] != -1) {
        (*ret_lcskpp_indices)[curr_recon++] = recon[i1];
      }
    }
  }
  if (print_debug == true) {
    printf ("Time for reconstruction: %f\n", (((float) (clock() - time_reconstruction))/CLOCKS_PER_SEC));
    fflush(stdout);
  }

  *ret_lcskpp_length = lcskpp_length;

  if (events)
    free(events);
  events = NULL;
  if (matches_starts)
    free(matches_starts);
  matches_starts = NULL;

  return 0;
}

int WriteLCSkDebug(std::string out_file, std::string qname, int64_t qlen, std::string rname, int64_t rlen, int32_t k, const uint128_t *packed_hits, int64_t n_hits, std::vector<int32_t> *lcskpp_indices, std::vector<int32_t> *cluster_ids) {
  if (n_hits == 0 || packed_hits == NULL) {
    return 1;
  }

  FILE *fp_out = fopen(out_file.c_str(), "w");
  if (fp_out == NULL) {
//    fprintf (stderr, "ERROR: Could not open file '%s' for writing!\n", out_file.c_str());
    return 2;
  }

  int32_t qid = get_lcsk128_qid(packed_hits[0]);
  int32_t rid = get_lcsk128_rid(packed_hits[0]);
  fprintf (fp_out, "%s\t%ld\t%ld\t%s\t%ld\t%ld\n", qname.c_str(), qid, qlen, rname.c_str(), rid, rlen);

  if (lcskpp_indices != NULL && cluster_ids == NULL) {
    /// In this case, output only the filtered indices.
    for (int64_t i=0; i<lcskpp_indices->size(); i++) {
      int32_t rpos = get_lcsk128_rpos(packed_hits[(*lcskpp_indices)[i]]) & 0x07FFFFFFF;
      int32_t qpos = get_lcsk128_qpos(packed_hits[(*lcskpp_indices)[i]]) & 0x07FFFFFFF;
      fprintf (fp_out, "%ld\t%ld\t0\n", qpos, rpos);
      fprintf (fp_out, "%ld\t%ld\t0\n", qpos+k, rpos+k);
    }
  } else if(lcskpp_indices != NULL && cluster_ids != NULL) {
    /// Output filtered indices and include cluster colors as well.
    for (int64_t i=0; i<lcskpp_indices->size(); i++) {
      int32_t rpos = get_lcsk128_rpos(packed_hits[(*lcskpp_indices)[i]]) & 0x07FFFFFFF;
      int32_t qpos = get_lcsk128_qpos(packed_hits[(*lcskpp_indices)[i]]) & 0x07FFFFFFF;
      fprintf (fp_out, "%ld\t%ld\t%d\n", qpos, rpos, (*cluster_ids)[i]);
      fprintf (fp_out, "%ld\t%ld\t%d\n", qpos+k, rpos+k, (*cluster_ids)[i]);
    }

  } else {
    /// Else, output everything.
    for (int64_t i=0; i<n_hits; i++) {
      int32_t rpos = get_lcsk128_rpos(packed_hits[i]); // & 0x07FFFFFFF;
      int32_t qpos = get_lcsk128_qpos(packed_hits[i]); // & 0x07FFFFFFF;
      fprintf (fp_out, "%ld\t%ld\t0\n", qpos, rpos);
      fprintf (fp_out, "%ld\t%ld\t0\n", qpos+k, rpos+k);
    }
  }

  fclose(fp_out);

  return 0;
}

void TestLCSk() {
  int32_t k = 12;
  int64_t qlen = 1000;
  int64_t rlen = 1000;

  std::vector<uint128_t> hits;

  /// Generate a truth line.
  for (int64_t i=200; i<(qlen - 200 - k); i+=(k+1)) {
    /// Linear.
    hits.push_back(pack_lcsk128(i, i, 0, 0));

    /// Experimenting with other functions.
//    int32_t x = i;
//    float i_flt = i;
//    float scaling = 3.14f / ((float) qlen);
//    int32_t y = (int32_t) (atan((i_flt / (qlen/2.0f) - 1.0f)*10) / scaling);
//    hits.push_back(pack128(x, y, 0, 0));
//    y = [math.atan((val / (1000/2.0) - 1.0)*10) for val in xflt]
  }

  /// Generate random seed hits.
//  srand (time(NULL));
  for (int64_t i=0; i<250; i++) {
    int64_t qpos = rand() % qlen;
    int64_t rpos = rand() % rlen;
    hits.push_back(pack_lcsk128(qpos, rpos, 0, 0));
  }

  /// Randomize the ordering.
  std::random_shuffle(hits.begin(), hits.end());

  std::sort(hits.begin(), hits.end());

  /// Write everything down into a file.
  WriteLCSkDebug("temp/test-raw.csv", "test1", qlen, "test2", rlen, k, &hits[0], hits.size(), NULL, NULL);

  /// Perform LCSk on the seed hits to obtain the original line.
  int64_t lcskpp_len = 0;
  std::vector<int32_t> lcskpp_indices;
  LCSkPacked(&hits[0], hits.size(), k, &lcskpp_len, &lcskpp_indices);

  /// Write the LCSk results to a file.
  WriteLCSkDebug("temp/test-lcsk.csv", "test1", qlen, "test2", rlen, k, &hits[0], hits.size(), &lcskpp_indices, NULL);
}

void TestLCSk2() {
  int32_t k = 12;
  int64_t qlen = 100;
  int64_t rlen = 100;

  std::vector<uint128_t> hits;

  /// Generate a truth line.
//  for (int64_t i=0; i<(qlen - k); i+=(k/2)) {
////    if ((i%4)==0 || (i%6)==0) {
//    if ((i%12) == 0) {
////      i += 7;
//      continue;
//    }
//    /// Linear.
//    hits.push_back(pack128(i, i, 0, 0));
//    printf ("i = %ld\n", i);
//    fflush(stdout);
//  }

  int32_t x = 0;
  x = 0; hits.push_back(pack_lcsk128(x, x, 0, 0));
  x = 6; hits.push_back(pack_lcsk128(x, x, 0, 0));
  x = 8; hits.push_back(pack_lcsk128(x, x, 0, 0));

  x = 30; hits.push_back(pack_lcsk128(x, x, 0, 0));
  x = 32; hits.push_back(pack_lcsk128(x, x, 0, 0));
  x = 42; hits.push_back(pack_lcsk128(x, x, 0, 0));

  x = 60; hits.push_back(pack_lcsk128(x, x, 0, 0));

  x = 80; hits.push_back(pack_lcsk128(x, x, 0, 0));

  /// Randomize the ordering.
  std::random_shuffle(hits.begin(), hits.end());

  std::sort(hits.begin(), hits.end());

  /// Write everything down into a file.
  WriteLCSkDebug("temp/test-raw.csv", "test1", qlen, "test2", rlen, k, &hits[0], hits.size(), NULL, NULL);

  /// Perform LCSk on the seed hits to obtain the original line.
  int64_t lcskpp_len = 0;
  std::vector<int32_t> lcskpp_indices;
  LCSkPacked(&hits[0], hits.size(), k, &lcskpp_len, &lcskpp_indices);

  /// Write the LCSk results to a file.
  WriteLCSkDebug("temp/test-lcsk.csv", "test1", qlen, "test2", rlen, k, &hits[0], hits.size(), &lcskpp_indices, NULL);
}

void TestLCSk3() {
  int32_t k = 12;

  for (int64_t N=100; N<=100000; N*=10) {
    int64_t qlen = N;
    int64_t rlen = N;
    std::vector<uint128_t> hits;
    /// Generate a truth line.
    int64_t step = (qlen - 400 - k) / (0.65f*N) + 1;
    for (int64_t i=200; i<(qlen - 200 - k); i+=step) {
      /// Linear.
      hits.push_back(pack_lcsk128(i, i, 0, 0));
    }

    /// Generate random seed hits.
  //  srand (time(NULL));
    for (int64_t i=0; i<(0.35f*N); i++) {
      int64_t qpos = rand() % qlen;
      int64_t rpos = rand() % rlen;
      hits.push_back(pack_lcsk128(qpos, rpos, 0, 0));
    }

    /// Randomize the ordering.
    std::random_shuffle(hits.begin(), hits.end());

    std::sort(hits.begin(), hits.end());

    /// Perform LCSk on the seed hits to obtain the original line.
    int64_t lcskpp_len = 0;
    std::vector<int32_t> lcskpp_indices;

    clock_t time1 = clock();
    printf ("Testing a batch of %ld (hits.size() = %ld) input points.\n", N, hits.size());
    fflush(stdout);
    LCSkPacked(&hits[0], hits.size(), k, &lcskpp_len, &lcskpp_indices, true);
    printf ("Time for batch of %ld (hits.size() = %ld) input points: %f\n", N, hits.size(), (((float) (clock() - time1))/CLOCKS_PER_SEC));
    printf ("\n");
    fflush(stdout);
  }
}
