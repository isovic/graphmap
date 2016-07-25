// #include "make_transcript.h"
#include "index_spaced_hash_fast.h"

int IndexSpacedHashFast::MakeTranscript_(std::string annotations_path, const SequenceFile &references, SequenceFile &transcripts) {
	std::map<std::string, std::vector<std::pair<std::string, char>>> seqToTrans;
	std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> transToExons;
	GenerateExons(annotations_path, seqToTrans, transToExons);
	transcripts.Clear();
	int64_t refN = references.get_sequences().size();
	int64_t id = 1;
	for(int64_t i = 0; i < refN; i++) {
		auto seq = references.get_sequences()[i];
		auto seqName = getSequenceName(*seq);
		for(auto &trans : seqToTrans[seqName]) {
			std::string transSeq = "";
			for(auto &coord : transToExons[trans.first]) {
				for(int64_t j = coord.first - 1; j < coord.second; j++) {
					transSeq += seq->get_data()[j];
				}
			}
			if(transSeq == "") continue;
			SingleSequence *s = new SingleSequence;

			auto seqLen = transSeq.size();
			int8_t *seq = new int8_t[seqLen];
			std::copy(transSeq.begin(), transSeq.end(), seq);

			auto headerLen = trans.first.size();
			char *header = new char[headerLen];
			std::copy(trans.first.begin(), trans.first.end(), header);

			s->InitHeaderAndDataFromAscii(header, headerLen, seq, seqLen, id++);
			if(trans.second == '-') {
				s->ReverseComplement();
			}
			transcripts.AddSequence(s, true);
			// outputSeq(header, headerLen, s->get_data(), seqLen);
		}
	}
	return 0;
}

int IndexSpacedHashFast::GenerateExons(
	std::string annotations_path,
	std::map<std::string, std::vector<std::pair<std::string, char>>> &seqToTrans,
	std::map<std::string, std::vector<std::pair<int64_t, int64_t>>> &transToExons
) {
	std::ifstream annotations(annotations_path);
	if(!annotations.is_open()) {
      FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", annotations_path.c_str());
      return 1;
    }
	std::string line;
	while(getline(annotations, line)) {
		auto fields = split(line, '\t');
		if(fields.size() < 9 || fields[2] != "exon") {
			continue;
		}
		std::string tid = getTID(fields[8]);
		if(transToExons[tid].empty()) {
			seqToTrans[split(fields[0], ' ')[0]].push_back(std::make_pair(tid, fields[6][0]));
		}
		int64_t left, right;
	    std::stringstream ss(fields[3] + " " + fields[4]);
		ss >> left >> right;
		transToExons[tid].push_back(std::make_pair(left, right));
	}
	return 0;
}

std::string IndexSpacedHashFast::getTID(std::string attributes) {
	for(auto s : split(attributes, ';')) {
		s = trim(s);
		auto keyValue = split(s, ' ');
		if(keyValue[0] == "transcript_id") {
			return split(keyValue[1], '"')[1];
		}
	}
	return "";
}

std::string IndexSpacedHashFast::trim(std::string s) {
	std::string ret = "";
	int i, j;
	for(i = 0; i < s.size(); i++) {
		if(!std::isspace(s[i])) {
			break;
		}
	}
	for(j = s.size() - 1; j >= 0; j--) {
		if(!std::isspace(s[j])) {
			break;
		}
	}
	for(; i <= j; i++) ret += s[i];
	return ret;
}

std::vector<std::string> IndexSpacedHashFast::split(std::string s, char c) {
	std::vector<std::string> v;
    std::string temp = "";
    for(int i = 0; i < s.size(); i++) {
        if(s[i] == c) {
            v.push_back(temp);
            temp = "";
        }
        else
            temp += s[i];
    }
    v.push_back(temp);
	return v;
}

std::string IndexSpacedHashFast::getSequenceName(SingleSequence &seq) {
	std::string name = "";
	for(int i = 0; i < seq.get_header_length() && seq.get_header()[i] != ' '; i++) {
		name += seq.get_header()[i];
	}
	return name;
}

void IndexSpacedHashFast::outputSeq(char *header, size_t headerLen, const int8_t *seq, size_t seqLen) {
	std::cout << ">";
	for(int i = 0; i< headerLen; i++)
		std::cout << header[i];
	for(int i = 0; i < seqLen; i++) {
		if(i % 80 == 0) std::cout << "\n";
		std::cout << seq[i];
	}
	std::cout << "\n";
}
