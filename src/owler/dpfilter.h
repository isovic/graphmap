#ifndef DP_FILTER_
#define DP_FILTER_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <stdint.h>

#include "owler/owler_data.h"

struct Khit {
	int32_t x = 0, y = 0;
	float c = 0.0f;
//	int32_t x, y;
//	float c;
};

inline float penalty(int32_t len_ref, int32_t len_query, const Khit *curr, const Khit *prev);
void dp_filter(int32_t len_ref, int32_t len_query, const std::vector<Khit> &lcsk_hits, std::vector<int32_t> &filtered_indices);
void OwlerDPFilter(std::vector<int> &lcskpp_indices, std::vector<SeedHit2> &seed_hits, int64_t ref_hits_start, int64_t ref_hits_end, int64_t seed_length, int32_t len_ref, int32_t len_query, std::vector<int32_t> &filtered_indices);

#endif
