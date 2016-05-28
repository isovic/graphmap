/*
 * rna.cc
 *
 *  Created on: May 27, 2016
 *      Author: isovic
 */

#include "graphmap/graphmap.h"
#include "graphmap/filter_anchors.h"

void ClearFile(std::string &out_file);
int VerbosePathGraphEntryCluster(std::string &out_file, const Region &region, const MappingResults& mapping_results, const std::vector<Index *> &indexes,
                                  const SingleSequence* read, const ProgramParameters* parameters);


// This function filters clusters in different mapping regions, to preserve only the ones which might be construed as valid RNA-seq mappings.
int GraphMap::RNAFilterClusters_(MappingData* mapping_data, const std::vector<Index *> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {

  // Commented out debug info. This is used to create the datasets containing clusters, meant for RNA-seq support development (knapsack).
//  #ifndef RELEASE_VERSION
//    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//      std::string temp_cluster_file = FormatString("temp/clusters/new-clusters-read-%ld.csv", read->get_sequence_id());
//      ClearFile(temp_cluster_file);
//    }
//  #endif

  int64_t read_len = read->get_sequence_length();

  // Clusters are stored in the intermediate mappings. Intermediate mapping corresponds to one processed region on the reference.
  for (int32_t i=0; i<mapping_data->intermediate_mappings.size(); i++) {
    // All info about the region is given here (such as: reference_id, start and end coordinates of the region, etc.).
    Region& region = mapping_data->intermediate_mappings[i]->region_data();
    int64_t ref_id = region.reference_id % indexes[0]->get_num_sequences_forward();     // If there are N indexed sequences, then the index contains 2*N sequences: first N are the forward strand, followed by the same N sequences reverse-complemented. The reference_id is then the absolute reference ID in the index, which means that if it refers to the reverse complement of the sequence, reference_id will be > N. Modulo needs to be taken.
    int64_t ref_len = indexes[0]->get_reference_lengths()[region.reference_id];

    // Each intermediate mapping contains a vector of clusters.
    MappingResults& mapping_results = mapping_data->intermediate_mappings[i]->mapping_data();
    auto& cluster_vector = mapping_results.clusters;

    // Commented out debug info. This is used to create the datasets containing clusters, meant for RNA-seq support development (knapsack).
//    #ifndef RELEASE_VERSION
//      if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//        std::string temp_cluster_file = FormatString("temp/clusters/new-clusters-read-%ld.csv", read->get_sequence_id());
//        VerbosePathGraphEntryCluster(temp_cluster_file, region, mapping_results, indexes, read, parameters);
//      }
//    #endif

    // Clusters can be accessed like so.
    for (int32_t j=0; j<cluster_vector.size(); j++) {
      Cluster& cluster = cluster_vector[j];
      // Members of Cluster contain these values:
      // cluster.query.start
      // cluster.query.end
      // cluster.ref.start
      // cluster.ref.end
      // cluster.num_anchors
      // cluster.coverage
      // If the cluster is supposed to be used, set cluster.valid to true, otherwise set it to false (mandatory; it will be true by default).
    }
  }

  return 0;
}



// Just a debug helper function.
void ClearFile(std::string &out_file) {
  FILE *fp_cluster_path = fopen(out_file.c_str(), "w");
  if (fp_cluster_path) {
    fclose(fp_cluster_path);
  }
}

// Just a debug helper function.
int VerbosePathGraphEntryCluster(std::string &out_file, const Region &region, const MappingResults& mapping_results, const std::vector<Index *> &indexes,
                                  const SingleSequence* read, const ProgramParameters* parameters) {
  FILE *fp = fopen(out_file.c_str(), "a");
  if (fp == NULL) { return 1; }

  // 1. Number of clusters, 2. Read ID, 3. Read len, 4. Target ID, 4. Target len, 5. Target reversed
  fprintf (fp, "#region\tID\tnum_clusters\tread_id\tread_len\tref_id\tref_len\tis_rev\n");
  fprintf (fp, "region\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%d\n", mapping_results.local_score_id, mapping_results.clusters.size(), read->get_sequence_id(), read->get_sequence_length(),
           (region.reference_id % indexes[0]->get_num_sequences_forward()), indexes[0]->get_reference_lengths()[region.reference_id], (int32_t) region.reference_id >= indexes[0]->get_num_sequences_forward());
  fprintf (fp, "#cluster\tID\tread_start\tread_end\tref_start\tref_end\tnum_anchors\tcoverage\n");

  int64_t current_cluster = 0;
  for (int64_t i=0; i<mapping_results.clusters.size(); i++) {
    if (mapping_results.clusters[i].valid == true && mapping_results.clusters[i].num_anchors > 0) {
      current_cluster += 1;
      // Cluster line:
      //  cluster_id qstart qend rstart rend num_anchors num_covered_bases
      fprintf (fp, "cluster\t%ld\t%ld\t%ld\t%ld\t%ld\t%d\t%d\n",
               current_cluster, mapping_results.clusters[i].query.start, mapping_results.clusters[i].query.end,
               mapping_results.clusters[i].ref.start, mapping_results.clusters[i].ref.end,
               mapping_results.clusters[i].num_anchors, mapping_results.clusters[i].coverage);
    }
  }
  fprintf (fp, "#\n");
  fclose(fp);
  return 0;
}

