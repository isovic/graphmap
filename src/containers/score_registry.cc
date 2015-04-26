/*
 * score_registry.cc
 *
 *  Created on: Jul 14, 2014
 *      Author: ivan
 */

#include "containers/score_registry.h"

ScoreRegistry::ScoreRegistry() {
  scores_id_ = 0;
}

ScoreRegistry::ScoreRegistry(const Region& region, int64_t scores_id) {
  set_region(region);
  set_scores_id(scores_id);
}

ScoreRegistry::~ScoreRegistry() {
  Clear();
}

void ScoreRegistry::Clear() {
//  registry_.clear();
  registry_entries_.Clear();
  scores_id_ = 0;
}

void ScoreRegistry::Add(Vertices &src_vertices, int64_t vertex_idx) {
//  registry_.push_back(vertex_data);
  registry_entries_.Add(src_vertices, vertex_idx);
}

void ScoreRegistry::Register(Vertices &src_vertices, int64_t vertex_idx) {
  if (src_vertices.registry_numbers[vertex_idx] < 0) { // || vertex_data.registry_number >= registry_.size()) {
    src_vertices.registry_numbers[vertex_idx] = registry_entries_.num_vertices;
    registry_entries_.Add(src_vertices, vertex_idx);

  }
  else {
    // Handle the case where a repeating kmer causes a 'jump' in the middle of an existing long path.
    // Edit 07.11.2014.: Because of the condition that a kmer needs to be within l iterations from the
    // vertex's path that it want's to extend, the kmer cannot hit it somewhere in the middle of the path.
    // It can only occur near the end of the path, and can only cause the path to have a more-or-less
    // even/uneven length in the reference and the query. For this reason, I think that forking a path
    // is perhaps not a good option, but instead to check for its ratio in query and in reference, and
    // choose to extend the path with the new kmer only if the ratio is closer to 1.0f.
    // For precaution sake, I'll keep the previous version here in comments.
//    if (vertex_data.covered_bases < registry_[vertex_data.registry_number].covered_bases) {
//      vertex_data.registry_number = registry_.size();
//      registry_.push_back(vertex_data);
//    } else {
//      registry_[vertex_data.registry_number] = vertex_data;
//    }

    int64_t registry_number = src_vertices.registry_numbers[vertex_idx];

    if ((src_vertices.num_kmers[vertex_idx] > registry_entries_.num_kmers[registry_number]) ||
        (src_vertices.num_kmers[vertex_idx] <= registry_entries_.num_kmers[registry_number] &&
            src_vertices.CalculateSuppress(vertex_idx) < registry_entries_.CalculateSuppress(registry_number))) {

        registry_entries_.CopyValuesFromOut(src_vertices, vertex_idx, registry_number);
    }
  }
}

std::string ScoreRegistry::VerboseToString() {
  std::stringstream ss;

  ss << "Num scores: " << registry_entries_.num_vertices << std::endl;

  for (int64_t i=0; i<registry_entries_.num_vertices; i++) {
    ss << "[" << i << "] (" << registry_entries_.VerboseToString(i) << ")" << std::endl;
  }

//  ss << std::endl;

  return ss.str();
}

const Region& ScoreRegistry::get_region() const {
  return region_;
}

void ScoreRegistry::set_region(const Region& region) {
  region_ = region;
}

int64_t ScoreRegistry::get_scores_id() const {
  return scores_id_;
}

void ScoreRegistry::set_scores_id(int64_t scoresId) {
  scores_id_ = scoresId;
}

const Vertices& ScoreRegistry::get_registry_entries() const {
  return registry_entries_;
}

void ScoreRegistry::Reserve(int64_t size) {
  registry_entries_.Reserve(size);
}

void ScoreRegistry::set_registry_entries(Vertices& registryEntries) {
  registry_entries_ = registryEntries;
}
