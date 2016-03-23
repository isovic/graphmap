/*
 * mapping_data.h
 *
 *  Created on: Mar 19, 2015
 *      Author: isovic
 */

#ifndef MAPPING_DATA_H_
#define MAPPING_DATA_H_

#include "log_system/log_system.h"
#include "alignment/local_realignment.h"
#include "utility/program_parameters.h"
#include "utility/utility_general.h"
#include "containers/region.h"
#include "index/index.h"
#include "index/index_hash.h"
#include "index/index_sa.h"
#include "index/index_spaced_hash.h"
#include "containers/vertices.h"
#include "utility/evalue.h"

//#define UNMAPPED_CODE_NO_VALID_GRAPH_PATHS  (1 << 0)

#define MAPPED_CODE_READ_UNPROCESSED_YET      (0)
#define MAPPED_CODE_UNIQUE_MAPPING            (1 << 0)
#define MAPPED_CODE_MULTIPLE_EQ_MAPPINGS      (1 << 1)

#define ITERATION_RESET_LIMIT ((int64_t) 0x1000000000000000)



struct ChromosomeBin {
  int64_t reference_id = 0;
  int64_t bin_id = 0;
  float bin_value = 0.0f;
};

struct bins_greater_than_key
{
    inline bool operator() (const ChromosomeBin& op1, const ChromosomeBin& op2) {
      if (op1.bin_value > op2.bin_value)
        return true;
      return false;
    }
};

class MappingData {
 public:
  MappingData();
  ~MappingData();

  Vertices vertices;
  std::vector<ChromosomeBin> bins;
  std::vector<PathGraphEntry *> intermediate_mappings;
  std::vector<PathGraphEntry *> final_mapping_ptrs;

  int64_t bin_size;
  int64_t num_seeds_over_limit;
  int64_t num_seeds_with_no_hits;
  int64_t num_seeds_errors;
  int64_t iteration;

  int64_t num_similar_mappings;                  // Number of found mapping positions with very similar (estimated) scores. E.g. to within some difference from the top mapping.
  int64_t num_same_mappings;
  int64_t avg_covered_bases_of_all_mappings;
  int64_t std_covered_bases_of_all_mappings;
  int64_t median_covered_bases_of_all_mappings;

  std::string unmapped_reason;

  int64_t num_region_iterations;
  int8_t mapping_quality;
  int64_t metagen_alignment_score;

  double time_region_selection;
  double time_mapping;
  double time_alignment;
  double time_region_seed_lookup;
  double time_region_hitsort;
  double time_region_conversion;
  double time_region_alloc;
  double time_region_counting;

  bool IsMapped();
  bool IsAligned();

  std::string VerboseFinalMappingsToString(const Index *index, const SingleSequence *read) const;
  std::string VerboseIntermediateMappingsToString(const Index *index, const SingleSequence *read) const;

 private:
  std::string VerboseMappingDataToString_(const std::vector<PathGraphEntry *> *mapping_data, const Index *index, const SingleSequence *read) const;

};

#endif /* MAPPING_DATA_H_ */
