#include <math.h>
#include "dpfilter.h"

//#define DPFILTER_DEBUG_I 24

inline double dist(double x1, double y1, double x2, double y2) {
  return (sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)));
}

//#define penalty penalty_indel
#define penalty penalty_angular1

inline float penalty_angular1(int32_t len_ref, int32_t len_query, const Khit *curr, const Khit *prev) {
  if (curr == NULL && prev == NULL) { /// Something wrong was specified.
    return 0;

  } else if (curr != NULL && prev == NULL) {  /// Start gap penalty.
    return -45;

  } else if (curr == NULL && prev != NULL) {  /// End gap penalty.
    return -45;

  } else {  /// A normal case where a maximum number of mismatches + remaining indels are taken (max diag. + indels).
    float angle = 0.0f;
    if (curr->x != prev->x) {
      angle = atan(((double) (curr->y - prev->y)) / ((double) (curr->x - prev->x))) * 180.0f / 3.14159f - 45.0f;
    } else {
      angle = 45;
      printf ("Tu sam 1!\n");
    }
    return -abs(angle);

  }
  return 0.0;
}

/// Parameter description:
///  curr must always be larger or equal to prev (sorting is expected). These should be derived from a LCSk routine.
inline float penalty_angular(int32_t len_ref, int32_t len_query, const Khit *curr, const Khit *prev) {
  if (curr == NULL && prev == NULL) { /// Something wrong was specified.
    return 0;
  } else if (curr != NULL && prev == NULL) {  /// Start gap penalty.
    int32_t l = curr->y - curr->x;
    if (l >= 0) {
      return -(dist(0, l, curr->x, curr->y) + l);
    } else {
      return -(dist(-l, 0, curr->x, curr->y) - l);
    }

  } else if (curr == NULL && prev != NULL) {  /// End gap penalty.
    return -(len_query - prev->y);

  } else {  /// A normal case where a maximum number of mismatches + remaining indels are taken (max diag. + indels).
    float angle = 0.0f;
    if (curr->x != prev->x) {
      angle = atan((curr->y - prev->y) / (curr->x - prev->x)) * 180.0f / 3.14159f - 45.0f;
    } else {
      angle = 45;
    }
    return -abs(angle);

  }
  return 0.0;
}

/// Parameter description:
///  curr must always be larger or equal to prev (sorting is expected). These should be derived from a LCSk routine.
inline float penalty_indel(int32_t len_ref, int32_t len_query, const Khit *curr, const Khit *prev) {
  int32_t ret = 0;
  if (curr == NULL && prev == NULL) { /// Something wrong was specified.
    return 0;

  } else if (curr != NULL && prev == NULL) {  /// Start gap penalty.
    return -std::min(curr->x, curr->y);

  } else if (curr == NULL && prev != NULL) {  /// End gap penalty.
    return -abs((float) len_query - prev->x);

  } else {  /// A normal case where a maximum number of mismatches + remaining indels are taken (max diag. + indels).
    return -(std::max(abs(curr->x - prev->x), abs(curr->y - prev->y)) - std::min(abs(curr->x - prev->x), abs(curr->y - prev->y)));

  }

  return ret;
}

void dp_filter(int32_t len_ref, int32_t len_query, const std::vector<Khit> &lcsk_hits, std::vector<int32_t> &filtered_indices) {
	std::vector<double> diags(lcsk_hits.size(), 0);
	std::vector<double> diffs(lcsk_hits.size(), 0);

	int32_t minx = 0, miny=0;
	if (lcsk_hits.size() > 0) {
		minx = lcsk_hits[0].x;
		miny = lcsk_hits[0].y;
	}
	for (int32_t i=0; i<lcsk_hits.size(); i++) {
		minx = std::min(minx, lcsk_hits[i].x);
		miny = std::min(miny, lcsk_hits[i].y);
	}
	
	int32_t k = 12;
	int32_t match = k;
/*	std::vector<Khit> */
	
	int32_t n = lcsk_hits.size();
//	std::vector<std::vector<int32_t> > dp;
//	dp.resize(n+1);
//	for (int32_t i=0; i<(n+1); i++) {
//	  dp[i].resize(n+1, 0);
//	}

	std::vector<int32_t> dp_max_id(n+2, -1);
	std::vector<float> dp_max_val(n+2, 0);
	dp_max_id[0] = 0;
	
	for (int32_t i=1; i<(n+1); i++) {
#ifdef DPFILTER_DEBUG_I
    if (i == DPFILTER_DEBUG_I) {
      printf ("i = %d\n", i);
      if (i > 1) {
        printf ("Prev: x = %d, y = %d\n", lcsk_hits[i-2].x, lcsk_hits[i-2].y);
      }
      printf ("Curr: x = %d, y = %d\n", lcsk_hits[i-1].x, lcsk_hits[i-1].y);
    }
#endif

		for (int32_t j=0; j<i; j++) {
			float dp = 0;
			if (j > 0) {
			  dp = match + penalty(len_ref, len_query, &lcsk_hits[i-1], &lcsk_hits[j-1]) + dp_max_val[j];
			} else {
        dp = penalty(len_ref, len_query, &lcsk_hits[i-1], NULL);
			}

			if (j == 0 || dp >= dp_max_val[i]) {
			  dp_max_val[i] = dp;
			  dp_max_id[i] = j;
			}

#ifdef DPFILTER_DEBUG_I
      if (i == DPFILTER_DEBUG_I) {
        float penalty_val = (j > 0) ? (dp - match - dp_max_val[j]) : dp;
        printf ("\t[j = %d] dp = %.2f, penalty = %.2f, dp_max_val[i] = %.2f, dp_max_val[j] = %.2f, dp_max_id[i] = %d, x = %d, y = %d\n", j, dp, penalty_val, dp_max_val[i], dp_max_val[j], dp_max_id[i], lcsk_hits[j].x, lcsk_hits[j].y);
      }
#endif
		}
	}
  for (int32_t j=1; j<(n+1); j++) {
    float dp = 0 + penalty(len_ref, len_query, NULL, &lcsk_hits[j-1]) + dp_max_val[j];
    if (j == 1 || dp >= dp_max_val[n+1]) {
      dp_max_val[n+1] = dp;
      dp_max_id[n+1] = j;
    }
#ifdef DPFILTER_DEBUG_I
    if ((n + 1) == DPFILTER_DEBUG_I) {
      float penalty_val = (j > 0) ? (dp - match - dp_max_val[j]) : dp;
      printf ("\t[j = %d] dp = %.2f, penalty = %.2f, dp_max_val[i] = %.2f, dp_max_id[i] = %d\n", j, dp, penalty_val, dp_max_val[n+1], dp_max_id[n+1]);
    }
#endif
  }

  filtered_indices.clear();
  int32_t recon_id = dp_max_id[n+1];
//#ifdef DEBUG_I
//  printf ("dp_max_id[n+1] = %d\n", dp_max_id[n+1]);
//#endif
  while (recon_id > 0) {
//#ifdef DEBUG_I
//    printf ("dp_max_id[recon_id] = %d\n", dp_max_id[recon_id]);
//    printf ("recon_id = %d", recon_id);
//#endif
//    if (dp_max_id[recon_id] == 0) { break; }

    filtered_indices.push_back(recon_id-1);
    recon_id = dp_max_id[recon_id];
//#ifdef DEBUG_I
//    printf ("-> %d\n", recon_id);
//#endif
  }

#ifdef DPFILTER_DEBUG_I
  for (int32_t i=0; i<(n+2); i++) {
    if (i > 0) {
      printf ("[%d] %.2f\t%d\tx = %d\ty = %d\n", i, dp_max_val[i], dp_max_id[i], lcsk_hits[i-1].x, lcsk_hits[i-1].y);
    } else {
      printf ("[%d] %.2f\t%d\n", i, dp_max_val[i], dp_max_id[i]);
    }
    fflush(stdout);
  }
#endif

}




inline float OwlerPenaltyAngular(int32_t len_ref, int32_t len_query, SeedHit2 *curr, SeedHit2 *prev) {
  if (curr == NULL && prev == NULL) { /// Something wrong was specified.
    return 0;

  } else if (curr != NULL && prev == NULL) {  /// Start gap penalty.
    return -45;

  } else if (curr == NULL && prev != NULL) {  /// End gap penalty.
    return -45;

  } else {  /// A normal case where a maximum number of mismatches + remaining indels are taken (max diag. + indels).
    float angle = 0.0f;
    if (curr->query_pos != prev->query_pos) {
      angle = atan(((double) (curr->ref_pos - prev->ref_pos)) / ((double) (curr->query_pos - prev->query_pos))) * 180.0f / 3.14159f - 45.0f;
    } else {
      angle = 45;
    }
    return -abs(angle);

  }
  return 0.0;
}

void OwlerDPFilter(std::vector<int> &lcskpp_indices, std::vector<SeedHit2> &seed_hits, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length, int32_t len_ref, int32_t len_query, std::vector<int32_t>& filtered_indices) {
  std::vector<double> diags(lcskpp_indices.size(), 0);
  std::vector<double> diffs(lcskpp_indices.size(), 0);

  int32_t minx = 0, miny=0;
  if (lcskpp_indices.size() > 0) {
    int64_t lcskpp_index = lcskpp_indices.at(0) + ref_hits_start;
    minx = (int32_t) seed_hits[lcskpp_index].query_pos;
    miny = (int32_t) seed_hits[lcskpp_index].ref_pos;
  }
  for (int32_t i=0; i<lcskpp_indices.size(); i++) {
    int64_t lcskpp_index = lcskpp_indices.at(i) + ref_hits_start;
    minx = std::min(minx, (int32_t) seed_hits[lcskpp_index].query_pos);
    miny = std::min(miny, (int32_t) seed_hits[lcskpp_index].ref_pos);
  }

  int32_t match = seed_length;

  int32_t n = lcskpp_indices.size();

  std::vector<int32_t> dp_max_id(n+2, -1);
  std::vector<float> dp_max_val(n+2, 0);
  dp_max_id[0] = 0;

  for (int32_t i=1; i<(n+1); i++) {
    int64_t lcskpp_index_i = lcskpp_indices.at(i-1) + ref_hits_start;

    for (int32_t j=0; j<i; j++) {
      int64_t lcskpp_index_j = lcskpp_indices.at(j-1) + ref_hits_start;

      float dp = 0;
      if (j > 0) {
        dp = match + OwlerPenaltyAngular(len_ref, len_query, &seed_hits[lcskpp_index_i], &seed_hits[lcskpp_index_j]) + dp_max_val[j];
      } else {
        dp = OwlerPenaltyAngular(len_ref, len_query, &seed_hits[lcskpp_index_i], NULL);
      }

      if (j == 0 || dp >= dp_max_val[i]) {
        dp_max_val[i] = dp;
        dp_max_id[i] = j;
      }
    }
  }
  for (int32_t j=1; j<(n+1); j++) {
    int64_t lcskpp_index_j = lcskpp_indices.at(j-1) + ref_hits_start;
    float dp = 0 + OwlerPenaltyAngular(len_ref, len_query, NULL, &seed_hits[lcskpp_index_j]) + dp_max_val[j];
    if (j == 1 || dp >= dp_max_val[n+1]) {
      dp_max_val[n+1] = dp;
      dp_max_id[n+1] = j;
    }
  }

  filtered_indices.clear();
  int32_t recon_id = dp_max_id[n+1];
  while (recon_id > 0) {
    filtered_indices.push_back(lcskpp_indices.at(recon_id-1));
    recon_id = dp_max_id[recon_id];
  }
}
