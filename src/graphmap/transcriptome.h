/*
 * transcriptome.h
 *
 *  Created on: Feb 6, 2017
 *      Author: isovic
 */

#ifndef SRC_GRAPHMAP_TRANSCRIPTOME_H_
#define SRC_GRAPHMAP_TRANSCRIPTOME_H_

#include <stdint.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "sequences/sequence_file.h"

namespace is {

class Transcriptome;

std::shared_ptr<Transcriptome> createTranscriptome();

class Transcriptome {
 public:
  friend std::shared_ptr<Transcriptome> createTranscriptome();
  ~Transcriptome();

  /* Parses exons and extracts regions from the given GTF file.
   */
  int LoadGTF(const std::string &gtf_path);

  /* Constructs transcriptome sequences from the preloaded GTF file.
   */
  std::shared_ptr<SequenceFile> GenerateTranscriptomeSeqs(const std::shared_ptr<SequenceFile> sequences);

  /* Generates a header for a SAM file. The header is composed of
   * genomic sequence names.
   */
  std::string GenerateSAMHeaders();

  const std::map<std::string, std::vector<std::pair<std::string, char> > >& get_genome_id_to_trans_id() const {
    return genome_id_to_trans_id_;
  }

  const std::map<std::string, std::vector<std::pair<int64_t, int64_t> > >& get_trans_id_to_exons() const {
    return trans_id_to_exons_;
  }

  const std::map<std::string, std::vector<std::pair<int64_t, int64_t> > >& get_trans_id_to_regions() const {
    return trans_id_to_regions_;
  }

  const std::map<std::string, std::pair<std::string, char>>& get_trans_id_to_genome_id() const {
    return trans_id_to_genome_id_;
  }

  const std::map<std::string, int64_t>& get_genome_id_to_len() const {
    return genome_id_to_len_;
  }

 private:
  std::string gtf_path_;

  // A map from genome (chromosome) name (e.g. header split to first space) to a vector containing all transcriptomes which can be generated from that chromosome.
  // Each pair is a (transcript_id, strand), where strand is either '+' or '-';
  std::map<std::string, std::vector<std::pair<std::string, char>>> genome_id_to_trans_id_;
  // Reverse map, to obtain the chromosome name when converting from transcriptome space back to genome space.
  // Second parameter of the pair is the orientation on the genome.
  std::map<std::string, std::pair<std::string, char>> trans_id_to_genome_id_;
  // A map from transcript_id to a vector containing pairs of coordinates. Each pair of coordinates presents one exon which makes the transcriptome.
  std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> trans_id_to_exons_;
  // A list of exons in such way that it combines overlapping exons into regions.
  std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> trans_id_to_regions_;
  // Length of each chromosome in genome space. Needed for reversing the mapping if transcriptome was reverse complemented.
  std::map<std::string, int64_t> genome_id_to_len_;

  Transcriptome();      // Private constructor, prevent memory leaks;
  Transcriptome(const Transcriptome&) = delete;
  const Transcriptome& operator=(const Transcriptome&) = delete;

  // Creates a transcriptome from a given reference sequence and a path to a file with gene annotations.
  // Parameters:
  // @param annotations_path Path to a GFF file (or another supported format) which contains the annotations of exonic regions.
  // @param references A SequenceFile object which contains reference sequences already loaded from disk.
  // @param transcripts A SequenceFile which will contain the generated transcriptomes.
  // @return 0 if everything went fine (C-style).
  int MakeTranscript_(const std::map<std::string, std::vector<std::pair<std::string, char>>> &genome_id_to_trans_id,
                      const std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> &trans_id_to_exons,
                      const std::shared_ptr<SequenceFile> references, std::shared_ptr<SequenceFile> transcripts) const;
  /** Resolves lists of exons in such way that it combines overlapping exons into regions.
   * Returns dict that maps transcript id to list of regions.
   * @param trans_id_to_exons A map from transcriptome ID (name) to a vector of exons which make this transcriptome.
   * @param trans_id_to_regions Generated return map from transcriptome ID (name) to a vector containing regions.
   * @return 0 if everything went fine (C-style).
   */
  int MakeRegions_(const std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> &trans_id_to_exons,
                   std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> &trans_id_to_regions) const;
  int ParseExons_(const std::string &annotations_path,
                  std::map<std::string, std::vector<std::pair<std::string, char>>> &genomeToTrans,
                  std::map<std::string, std::pair<std::string, char>> &transIdToGenomeId,
                  std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> &transToExons) const;
  void HashGenomeLengths_(const std::shared_ptr<SequenceFile> sequences, std::map<std::string, int64_t> &rlens) const;
  std::string trim_(std::string s) const;
  std::vector<std::string> split_(std::string s, char c) const;
  std::string getSequenceName_(const SingleSequence &seq) const;
  std::string getTID_(const std::string &chr_name, const std::string &attributes) const;
//  void outputSeq_(char *header, size_t headerLen, const int8_t *seq, size_t seqLen) const;

};

} /* namespace is */

#endif /* SRC_GRAPHMAP_TRANSCRIPTOME_H_ */
