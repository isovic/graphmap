/*
 * mhap.h
 *
 *  Created on: May 29, 2016
 *      Author: isovic
 */

#ifndef SRC_MHAP_H_
#define SRC_MHAP_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <map>

#include "log_system/log_system.h"
#include "sequences/sequence_file.h"
#include "utility/utility_general.h"

class OverlapLine {
 public:
  OverlapLine() :
    Aid(0), Bid(0), Aname(""), Bname(""), perc_err(0.0), shared_minmers(0), Arev(0), Astart(0), Aend(0), Alen(0), Brev(0), Bstart(0), Bend(0), Blen(0) {
  }

  OverlapLine(int64_t _Aid, int64_t _Bid, std::string _Aname, std::string _Bname, double _perc_err, int64_t _shared_minmers,
              int64_t _Arev, int64_t _Astart, int64_t _Aend, int64_t _Alen,
              int64_t _Brev, int64_t _Bstart, int64_t _Bend, int64_t _Blen) :
    Aid(_Aid), Bid(_Bid), Aname(_Aname), Bname(_Bname), perc_err(_perc_err), shared_minmers(_shared_minmers),
    Arev(_Arev), Astart(_Astart), Aend(_Aend), Alen(_Alen), Brev(_Brev), Bstart(_Bstart), Bend(_Bend), Blen(_Blen) {
  }

  int ParseMHAP(std::string &line) {
    std::istringstream iss(line);
    if (!(iss >> Aid >> Bid >> perc_err >> shared_minmers >> Arev >> Astart >> Aend >> Alen >> Brev >> Bstart >> Bend >> Blen)) {
      ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Overlaps are not formatted in the MHAP format. Exiting.");
    }
    Aname = FormatString("%ld", Aid);
    Bname = FormatString("%ld", Bid);
    return 0;
  }

  std::string WriteMHAP() {
    std::stringstream ss;
    ss << Aid << " " << Bid << " " << perc_err << " " << shared_minmers << " " << Arev << " " << Astart << " " << Aend << " " << Alen << " " << Brev << " " << Bstart << " " << Bend << " " << Blen;
    return ss.str();
  }

  int ParsePAF(std::string &line, const std::map<std::string, int64_t> &qname_to_ids) {
    std::istringstream iss(line);
    std::string tempAstrand, cm;
    int64_t num_residue_matches=0, aln_block_len=0, mapq=0;
    if (!(iss >> Aname >> Alen >> Astart >> Aend >> tempAstrand >> Bname >> Blen >> Bstart >> Bend >> num_residue_matches >> aln_block_len >> mapq >> cm)) {
      ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Overlaps are not formatted in the PAF format. Exiting.");
    }

    if (tempAstrand == "+") {
      Arev = 0; Brev = 0;
    } else {
      Arev = 0; Brev = 1;
    }

    // Convert the scores.
    perc_err = ((double) num_residue_matches) / ((double) aln_block_len);
    std::size_t found = cm.find_last_of(":");
    sscanf (cm.substr(found+1).c_str(), "%ld", &shared_minmers);

    // Find the ID of A.
    auto it_a = qname_to_ids.find(Aname);
    if (it_a == qname_to_ids.end()) {
      ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Could not find qname '%s' in the input reads file! Exiting.", Aname.c_str());
    }
    Aid = it_a->second + 1;

    // Find the ID of B.
    auto it_b = qname_to_ids.find(Bname);
    if (it_b == qname_to_ids.end()) {
      ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Could not find qname '%s' in the input reads file! Exiting.", Bname.c_str());
    }
    Bid = it_b->second + 1;

    return 0;
  }

  void Switch() {
    OverlapLine T = *this;
    Aid = T.Bid; Bid = T.Aid;
    Arev = T.Brev; Astart = T.Bstart; Aend = T.Bend; Alen = T.Blen;
    Brev = T.Arev; Bstart = T.Astart; Bend = T.Aend; Blen = T.Alen;
  }

  void ReverseComplement() {
    OverlapLine T = *this;

    Arev= 1 - T.Arev;
    Astart = T.Alen - T.Aend - 1;
    Aend = T.Alen - T.Astart - 1;

    Brev = 1 - T.Brev;
    Bstart = T.Blen - T.Bend - 1;
    Bend = T.Blen - T.Bstart - 1;
  }

  std::string Verbose() const {
    std::stringstream ss;
    ss << Aid << " " << Bid << " " << perc_err << " " << shared_minmers << " " << Arev << " " << Astart << " " << Aend << " " << Alen << " " << Brev << " " << Bstart << " " << Bend << " " << Blen;
    return ss.str();
  }

  int CheckConstraints(double max_dist_ratio) const {
    double Adist = Aend - Astart;
    double Bdist = Bend - Bstart;
    double ratio = (Adist > Bdist) ? (1.0f - Bdist / Adist) : (1.0f - Adist/Bdist);
    if (ratio > max_dist_ratio) { return 1; }
    return 0;
  }

  int64_t Aid, Bid;
  std::string Aname, Bname;
  double perc_err;
  int64_t shared_minmers;
  int64_t Arev, Astart, Aend, Alen;     // start is zero-based, end points to a position right after the last inclusive base.
  int64_t Brev, Bstart, Bend, Blen;
};

#endif /* SRC_MHAP_H_ */
