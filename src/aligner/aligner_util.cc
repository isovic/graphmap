#include "aligner_util.hpp"
#include "assert.h"

#include <sstream>
#include <string>

namespace is {

std::vector<int8_t> ConvertSeqAlphabet(const int8_t* seq, size_t seqlen, const uint8_t* conv_table) {
  std::vector<int8_t> ret(seqlen + 33); // 32 for gaba
  for (size_t i=0; i<seqlen; i++) {
    ret[i] = (int8_t) conv_table[(uint8_t) seq[i]];
  }
  return ret;
}

std::vector<is::CigarOp> ConvertBasicToExtCIGAR(const char* qseq, int64_t qlen,
                                                const char* tseq, int64_t tlen,
                                                const std::vector<is::CigarOp>& basic_cigar) {
  std::vector<is::CigarOp> ret;

  int64_t qpos = 0, tpos = 0;
  for (size_t i=0; i<basic_cigar.size(); i++) {
    char op = basic_cigar[i].op;
    int64_t count = basic_cigar[i].count;

    if (op != 'M') {
      ret.push_back(basic_cigar[i]);

      if (op == 'I' || op == 'S') {
        qpos += count;
      }
      if (op == 'D' || op == 'N') {
        tpos += count;
      }
    } else {
      char prev_m = 0;
      int64_t curr_count = 0;
      for (int64_t j=0; j<count; j++) {
        char curr_m = (qseq[qpos] == tseq[tpos]) ? '=' : 'X';
        if (j == 0) { prev_m = curr_m; }
        if (curr_m == prev_m) {
          curr_count += 1;
        } else {
          ret.push_back(is::CigarOp(prev_m, curr_count));
          prev_m = curr_m;
          curr_count = 1;
        }
        qpos += 1;
        tpos += 1;
      }
      if (curr_count > 0) {
        ret.push_back(is::CigarOp(prev_m, curr_count));
      }
    }
  }

  return ret;
}

int64_t EditDistFromExtCIGAR(const std::vector<is::CigarOp>& extended_cigar) {
  int64_t edit_dist = 0;
  for (size_t i=0; i<extended_cigar.size(); i++) {
    char op = extended_cigar[i].op;
    assert(op != 'M');
    if (op == 'X' || op == 'I' || op == 'D') {
      edit_dist += extended_cigar[i].count;
    }
  }
  return edit_dist;
}

std::vector<is::CigarOp> ExtractCigarBetweenQueryCoords(const std::vector<is::CigarOp>& cigar, int64_t qstart, int64_t qend, int64_t *cigar_length) {
  std::vector<is::CigarOp> ret;

  int64_t qpos = 0;

  int lengthOfRef = 0;

  for (auto& c: cigar) {

    int64_t qpos_next = (c.op == 'M' || c.op == '=' || c.op == 'X' || c.op == 'I' || c.op == 'S') ? (qpos + c.count) : qpos;

    if (qpos > qend) { break; }

    if (qpos_next < qstart) {
      qpos = qpos_next;
      continue;
    }

    int64_t b = 0, e = c.count;

    if (qstart >= qpos && qstart < qpos_next) { b = qstart - qpos; }
    if (qend >= qpos && qend < qpos_next) { e = qend - qpos; }

    if ((e - b) > 0) {
      ret.emplace_back(is::CigarOp(c.op, (e - b)));

      if (c.op != 'I') {
    	  lengthOfRef += (e - b);
      }
    }

    qpos = qpos_next;
  }

  *cigar_length = lengthOfRef;

  return ret;
}

std::string CigarToString(const std::vector<is::CigarOp>& cigar) {
  std::stringstream ss;
  for (size_t i=0; i<cigar.size(); i++) {
    ss << cigar[i].count << cigar[i].op;
  }
  return ss.str();
}

}
