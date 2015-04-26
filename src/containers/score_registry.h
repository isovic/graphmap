/*
 * score_registry.h
 *
 *  Created on: Jul 14, 2014
 *      Author: ivan
 */

#ifndef SCORE_REGISTRY_H_
#define SCORE_REGISTRY_H_

#include <string>
#include <vector>
#include <list>
#include <sstream>
#include "containers/region.h"
#include "sequences/single_sequence.h"
#include "sequences/sequence_file.h"
#include "containers/vertices.h"

class ScoreRegistry {
 public:
  ScoreRegistry();
  ScoreRegistry(const Region& region, int64_t scores_id);
  ~ScoreRegistry();

  /// Empties the registry and sets all values to zero.
  void Clear();

  /// Simply appends the data to the end of the registry and updates the top score.
  /// No additional checks are performed.
  void Add(Vertices &src_vertices, int64_t vertex_idx);

  /// If the data has a registry number >= 0, then the entry with that index will be updated.
  /// Otherwise, if registry number < 0 or if the suppress is smaller than the existing one,
  /// the new data will only be appended to the end of the registry, and its registry number
  /// will be updated.
  void Register(Vertices &src_vertices, int64_t vertex_idx);

  // Allocates space for vertices.
  void Reserve(int64_t size);

  /// Formats the debug verbose to a std::string.
  std::string VerboseToString();
  const Region& get_region() const;
  void set_region(const Region& region);
  int64_t get_scores_id() const;
  void set_scores_id(int64_t scoresId);
  const Vertices& get_registry_entries() const;
  void set_registry_entries(Vertices& registryEntries);

 private:
  Vertices registry_entries_;
  Region region_;
  int64_t scores_id_;
};

#endif /* SCORE_REGISTRY_H_ */
