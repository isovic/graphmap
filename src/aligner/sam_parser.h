#ifndef SRC_SAM_PARSER_H_
#define SRC_SAM_PARSER_H_

#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>

namespace is {

#define is_cigar_op(x)  (x == 'M' || x == '=' || x == 'X' || x == 'I' || x == 'D' || x == 'S' || x == 'H')
#define is_cigar_match(x)  (x == 'M' || x == '=' || x == 'X')
#define is_cigar_ins(x)  (x == 'I')
#define is_cigar_del(x)  (x == 'D')
#define is_cigar_soft(x)  (x == 'S')
#define is_cigar_hard(x)  (x == 'H')
#define is_cigar_ref(x)  (x == 'M' || x == '=' || x == 'X' || x == 'D')
#define is_cigar_read(x)  (x == 'M' || x == '=' || x == 'X' || x == 'I' || x == 'S')

// class CigarOp {
//  public:
//   char op = '-';
//   int32_t count = 0;
//   int64_t pos_ref = -1;     // Relative to the pos_ field of the corresponding SequenceAlignment object. pos_ref starts from zero, eventhough the actuall alignment starts at an arbitrary position on the reference.
//   int64_t pos_query = - 1;

//   CigarOp() { }
//   CigarOp(char _op, int32_t _count, int64_t _pos_ref, int64_t _pos_query) : op(_op), count(_count), pos_ref(_pos_ref), pos_query(_pos_query) { }

// };

/** @brief A container for a single CIGAR operation.
 *
 */
class CigarOp {
 public:
  CigarOp() : op(0), count(0) { }
  CigarOp(char _op, int32_t _count) : op(_op), count(_count) { }
  CigarOp(const CigarOp& t) : CigarOp(t.op, t.count) { }
  ~CigarOp() { }
  CigarOp& operator=(const CigarOp t) {
    op = t.op;
    count = t.count;
    return *this;
  }
  std::string get() { std::stringstream ss; ss << count << op; return ss.str(); }

  char op;
  int64_t count;
};


int SplitCigar(const std::string &cigar_str, std::vector<CigarOp>& ret);
int64_t CalcReferenceLengthFromCigar(const std::vector<CigarOp>& split_cigar);

class SamLine {
 public:
  SamLine();
  SamLine(const std::string& line);
//   SamLine();
// SequenceAlignment::SequenceAlignment(uint32_t _flag, std::string &rname, int64_t pos, int32_t mapq, std::string &cigar_string, std::string &rnext, int64_t pnext, int64_t tlen, std::vector<std::string> &optional)
// : flag_(flag), rname_(rname), pos_(pos), mapq_(mapq), rnext_(rnext), pnext_(pnext), tlen_(tlen), optional_(optional) {
//   SplitCigar(cigar_string, cigar_);
//   ProcessOptional();
// }
  ~SamLine();

  int ParseLine(const std::string& line);
  std::string YieldString();
  bool IsMapped();
  bool IsReverse();
  int FindAlignmentPosition(int64_t& q_start, int64_t& q_end,
                                    int64_t& r_start, int64_t& r_end);

  std::string qname;      // Field #1.
  uint32_t flag;          // Field #2.
  std::string rname;      // Field #3.
  int64_t pos;            // Field #4.
  int32_t mapq;           // Field #5.
  //  std::string cigar;  // Field #6.
  std::vector<CigarOp> cigar;
  std::string rnext;      // Field #7.
  int64_t pnext;          // Field #8.
  int64_t tlen;           // Field #9.
  std::string seq;        // Field #10.
  std::string qual;       // Field #11.

  // Optional fields in the SAM format:
  int64_t as;           // Alignment score.
  double evalue;        // E-value. There is no dedicated field in the SAM format, but GraphMap uses ZE to specify the E-value.
  std::vector<std::string> optional;  // Raw values (strings) of optional fields, not explicitly converted to expected values;



 private:
  void Tokenize_(const std::string& str, const char delimiter, std::vector<std::string>& words);

};

}

#endif
