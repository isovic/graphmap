/*
 * rna.cc
 *
 *  Created on: May 27, 2016
 *      Author: isovic
 */

#include "graphmap/graphmap.h"
#include "graphmap/filter_anchors.h"

#include <set>

#define VERBOSE_RNASEQ_DEBUG

void ClearFile(std::string &out_file);
int VerbosePathGraphEntryCluster(std::string &out_file, const Region &region, const MappingResults& mapping_results, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes,
                                  const SingleSequence* read, const ProgramParameters* parameters);

struct ClusterGroup {
	int start;
	int stop;
};

struct MyCluster {
	int64_t readStart;
	int64_t readEnd;
	int64_t refStart;
	int64_t refEnd;
	int coveredBases;
	int regionI;
	int clusterJ;
	int x1, y1;
	int x2, y2;

	MyCluster(int64_t _readStart, int64_t _readEnd, int64_t _refStart, int64_t _refEnd, int _coveredBases,
			int _regionI, int _clusterJ) :
		readStart(_readStart), readEnd(_readEnd), refStart(_refStart), refEnd(_refEnd),
		coveredBases(_coveredBases), regionI(_regionI), clusterJ(_clusterJ) {}

	friend bool operator <(const MyCluster& a, const MyCluster& b) {
		if (a.readStart != b.readStart) {
			return a.readStart < b.readStart;
		}
		if (a.refStart != b.refStart) {
			return a.refStart < b.refStart;
		}
		if (a.readEnd != b.readEnd) {
			return a.readEnd < b.readEnd;
		}
		if (a.refEnd != b.refEnd) {
			return a.refEnd < b.refEnd;
		}
		return a.coveredBases < b.coveredBases;
	}
};

#define MAXN 3001
#define STRANDS 2

//Margin for overlap of neighboring clusters.
#define FUZZ 10

//int backtrack[STRANDS][MAXN];
//int dp[STRANDS][MAXN];

void calculateDP(int **dp, int **backtrack, std::vector<MyCluster> *clusters) {
	for (int strand = 0; strand < STRANDS; strand++) {
		if (clusters[strand].size() == 0) {
			continue;
		}
		for (int n = clusters[strand].size(), x = n; x >= 0; x--) {
			int index = x - 1;
			int64_t xLeft = (index == -1) ? 0 : clusters[strand][index].readEnd, yLeft = (index == -1) ? 0 : clusters[strand][index].refEnd;
			int &ref = dp[strand][x];
			for (int i = x; i < n; i++) {
				if (clusters[strand][i].readStart < xLeft || clusters[strand][i].refStart < yLeft) {
					continue;
				}
				int tmp = clusters[strand][i].coveredBases + dp[strand][i + 1];
				if (tmp > ref) {
					backtrack[strand][x] = i + 1;
					ref = tmp;
				}
			}
		}
	}
}

void empty_clusters(std::vector<MyCluster> *clusters) {
  for (int strand = 0; strand < STRANDS; strand++) {
  	clusters[strand].clear();
  }
}

void GraphMap::RNAFilterClustersWithHugeGaps_(MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {
    std::vector<Cluster> clustersTmp;
    std::vector<ClusterGroup> clustersGroups;

    for (int32_t i=0; i<mapping_data->intermediate_mappings.size(); i++) {
  	  MappingResults& mapping_results = mapping_data->intermediate_mappings[i]->mapping_data();
  	  auto& cluster_vector = mapping_results.clusters;
      if (mapping_data->intermediate_mappings[i]->mapping_data().is_mapped == true) {
        for (int32_t j=0; j<cluster_vector.size(); j++) {
          Cluster& cluster = cluster_vector[j];
          clustersTmp.push_back(cluster);
        }
      }
    }

    if (clustersTmp.size() > 1) {
	    int start = 0;
	    int stop = 0;
	    Cluster currentCluster = clustersTmp[0];

	    for (int i = 1; i < clustersTmp.size(); i++) {
	        Cluster cluster = clustersTmp[i];
	        // Currently set 100000 as gap size, this should be discussed
	        if (abs(cluster.ref.start - currentCluster.ref.end) < 100000) {
	        		stop = i;
	        } else {
	        		ClusterGroup group;
	        		group.start = start;
	        		group.stop = stop;
	        		clustersGroups.push_back(group);
	        		start = i;
	        		stop = i;
	        }
			currentCluster = cluster;
	    }

		ClusterGroup group;
		group.start = start;
		group.stop = stop;
		clustersGroups.push_back(group);

		int maxCoverage = 0;
		ClusterGroup bestGroup = clustersGroups[0];

	    for (int i = 0; i < clustersGroups.size(); i++) {
	    		int coverage = 0;
	    		ClusterGroup tmpGroup;
	    		for(int j = tmpGroup.start; j < tmpGroup.stop; j++) {
	    			coverage += clustersTmp[j].coverage;
	    		}
	    		if(coverage > maxCoverage) {
	    			maxCoverage = coverage;
	    			bestGroup = tmpGroup;
	    		}
	    }

		int index = 0;

	    for (int32_t i=0; i<mapping_data->intermediate_mappings.size(); i++) {
	  	  MappingResults& mapping_results = mapping_data->intermediate_mappings[i]->mapping_data();
	  	  auto& cluster_vector = mapping_results.clusters;
	      if (mapping_data->intermediate_mappings[i]->mapping_data().is_mapped == true) {
	        for (int32_t j=0; j<cluster_vector.size(); j++) {
	          Cluster& cluster = cluster_vector[j];
	          if(index < bestGroup.start || index > bestGroup.stop) {
	        	  	  cluster.valid = false;
	          }
	          index += 1;
	        }
	      }
	    }

	    int ret_value_rna_fin2 = CleanupIntermediateMappings_(mapping_data, indexes, read, parameters);
    }
}

// This function filters clusters in different mapping regions, to preserve only the ones which might be construed as valid RNA-seq mappings.
int GraphMap::RNAFilterClusters_(MappingData* mapping_data, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes, const SingleSequence* read, const ProgramParameters* parameters) {

  // Commented out debug info. This is used to create the datasets containing clusters, meant for RNA-seq support development (knapsack).
//  #ifndef RELEASE_VERSION
//    if (parameters->verbose_level > 5 && read->get_sequence_id() == parameters->debug_read) {
//      std::string temp_cluster_file = FormatString("temp/clusters/new-clusters-read-%ld.csv", read->get_sequence_id());
//      ClearFile(temp_cluster_file);
//    }
//  #endif

  int64_t read_len = read->get_sequence_length();

  std::set<MyCluster> clusterSet[STRANDS];
	std::vector<MyCluster> clusters[STRANDS];

	if (mapping_data->intermediate_mappings.size() == 0) {
		//printf ("aloooo!!! polje vectora je 0!\n");
		//empty_clusters();
		return 0;
	}
  std::vector<Cluster*> cluster_data[mapping_data->intermediate_mappings.size()];
  //printf ("\npopis clustera po regijama: \n");
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

    int reverseStrand = (int32_t) region.reference_id >= indexes[0]->get_num_sequences_forward();

		#ifdef VERBOSE_RNASEQ_DEBUG
    	LOG_DEBUG_SPEC_NO_HEADER ("\ni: %d, cluster_vector.size(): %d\n", i, cluster_vector.size());
		#endif

    // Clusters can be accessed like so.
    for (int32_t j=0; j<cluster_vector.size(); j++) {
			#ifdef VERBOSE_RNASEQ_DEBUG
				LOG_DEBUG_SPEC_NO_HEADER("\tj: %d\n", j);
      	  	#endif

			Cluster& cluster = cluster_vector[j];
			clusterSet[reverseStrand].insert(MyCluster(cluster.query.start, std::max(cluster.query.start, cluster.query.end - FUZZ),
																								 cluster.ref.start, std::max(cluster.ref.start, cluster.ref.end - FUZZ),
																								 cluster.coverage, i, j));
			cluster.valid = false;
			cluster_data[i].push_back(&cluster);

			#ifdef VERBOSE_RNASEQ_DEBUG
				LOG_DEBUG_SPEC_NO_HEADER ("readStart: %lld, readEnd: %lld, refStart: %lld, refEnd: %lld, coveredBases: %d, i: %d, j: %d\n",
																	cluster.query.start, cluster.query.end, cluster.ref.start, cluster.ref.end,
																	cluster.coverage, i, j);
			#endif
    }
  }

	#ifdef VERBOSE_RNASEQ_DEBUG
		LOG_DEBUG_SPEC("Initialized cluster data. Initializing sets.\n");
	#endif

	int emptySets = 0;
  for (int strand = 0; strand < STRANDS; strand++) {
		#ifdef VERBOSE_RNASEQ_DEBUG
  		LOG_DEBUG_SPEC_NO_HEADER ("\nclusterSet[%d].size(): %d\n", strand, clusterSet[strand].size());
		#endif

  	if (clusterSet[strand].size() == 0) {
  		emptySets++;
  		continue;
  	}
		for (const MyCluster &cluster : clusterSet[strand]) {
			clusters[strand].push_back(cluster);
		}
	}

	#ifdef VERBOSE_RNASEQ_DEBUG
		LOG_DEBUG_SPEC("Moving from set to vector done.\n");
	#endif

	if (emptySets == STRANDS) {
		LOG_DEBUG_SPEC("The cluster set is empty, nothing to do.\n");
		return 0;
	}

	#ifdef VERBOSE_RNASEQ_DEBUG
		LOG_DEBUG_SPEC("Initializing the dp and backtrack matrices.\n");
	#endif

  int **dp = (int **) calloc(STRANDS, sizeof(int));
  int **backtrack = (int **) calloc(STRANDS, sizeof(int));
  for (int strand = 0; strand < STRANDS; strand++) {
  	dp[strand] = (int *) calloc(MAXN, sizeof(int));
  	backtrack[strand] = (int *) calloc(MAXN, sizeof(int));
  	for (int i = 0; i < MAXN; i++) {
  		backtrack[strand][i] = -1;
  	}
  }

	#ifdef VERBOSE_RNASEQ_DEBUG
		LOG_DEBUG_SPEC("Running calculateDP.\n");
	#endif

  calculateDP(dp, backtrack, clusters);

	#ifdef VERBOSE_RNASEQ_DEBUG
		LOG_DEBUG_SPEC("calculateDP done.\n");
	#endif

	int strand = (dp[0][0] > dp[1][0]) ? 0 : 1;

	std::set<int>	done;
  int curr = backtrack[strand][0];

	#ifdef VERBOSE_RNASEQ_DEBUG
		LOG_DEBUG_SPEC("Backtracking:\n");
	#endif

  while (true) {
  	if (curr == -1) {
  		break;
  	}
  	if (done.count(curr) > 0) {
			#ifdef VERBOSE_RNASEQ_DEBUG
				LOG_DEBUG_SPEC("Already been here!\n");
  			LOG_DEBUG_SPEC ("number of regions: %d\n", mapping_data->intermediate_mappings.size());
			#endif
  		return 0;
  	}
  	done.insert(curr);
  	int currClusterIndex = curr - 1;
  	if (curr == backtrack[strand][curr]) {
			#ifdef VERBOSE_RNASEQ_DEBUG
				LOG_DEBUG_SPEC("Next one is the same as current one (circular backtrack?)!\n");
  			LOG_DEBUG_SPEC ("number of regions: %d\n", mapping_data->intermediate_mappings.size());
			#endif
  		return 0;
  	}
  	int i = clusters[strand][currClusterIndex].regionI;
  	int j = clusters[strand][currClusterIndex].clusterJ;
		#ifdef VERBOSE_RNASEQ_DEBUG
  	LOG_DEBUG_SPEC_NO_HEADER ("i: %d, j: %d\n", i, j);
  	LOG_DEBUG_SPEC_NO_HEADER ("curr: %d\n", curr);
  	LOG_DEBUG_SPEC_NO_HEADER ("readStart: %lld, readEnd: %lld, refStart: %lld, refEnd: %lld, coveredBases: %d, i: %d, j: %d\n",
  		clusters[strand][currClusterIndex].readStart,
  		clusters[strand][currClusterIndex].readEnd,
  		clusters[strand][currClusterIndex].refStart,
  		clusters[strand][currClusterIndex].refEnd,
  		clusters[strand][currClusterIndex].coveredBases,
  		clusters[strand][currClusterIndex].regionI,
  		clusters[strand][currClusterIndex].clusterJ);
		#endif
  	////#XY printf ("backtrack[strand][curr]: %d\n", backtrack[strand][curr]);

  	cluster_data[i][j]->valid = true;
  	if (backtrack[strand][curr] == -1) {
  		break;
  	}
  	curr = backtrack[strand][curr];
  }

	#ifdef VERBOSE_RNASEQ_DEBUG
		LOG_DEBUG_SPEC("RNAFilterClusters_ done.\n");
	#endif

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
int VerbosePathGraphEntryCluster(std::string &out_file, const Region &region, const MappingResults& mapping_results, std::vector<std::shared_ptr<is::MinimizerIndex>> &indexes,
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

