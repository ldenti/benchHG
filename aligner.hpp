#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <map>
#include <regex>
#include <string>
#include <vector>

#include <bdsg/hash_graph.hpp>

#include <align.hpp>
#include <graphLoad.hpp>

using namespace std;

struct Alignment {
  int id;
  int l;
  int il;
  string cigar;
  vector<int> path;
  float score;

  void set_score() {
    // FIXME: is this the best way to do this?
    map<char, float> OPs;
    OPs['='] = 0.0;
    OPs['X'] = 0.0;
    OPs['I'] = 0.0;
    OPs['D'] = 0.0;
    regex word_regex("([0-9]+[=XID])");
    auto words_begin = sregex_iterator(cigar.begin(), cigar.end(), word_regex);
    auto words_end = std::sregex_iterator();
    for (sregex_iterator i = words_begin; i != words_end; ++i) {
      smatch match = *i;
      string match_str = match.str();
      char op = match_str.back();
      match_str.pop_back();
      int opl = stoi(match_str);
      OPs[op] += opl;
    }
    il = OPs['X'] + OPs['='] + OPs['I'];
    // score = 1.0 - (OPs['='] + OPs['D'] + OPs['X'] + abs(l - il)) / l; //
    // CHECKME
    score =
        OPs['='] / (OPs['='] + OPs['D'] + OPs['X'] + OPs['I'] + abs(l - il));
  }
};

class Aligner {
public:
  psgl::graphLoader graph;
  vector<string> queries;
  vector<psgl::BestScoreInfo> bestScoreVector;
  psgl::Parameters parameters;
  vector<Alignment> alignments;

public:
  Aligner(const bdsg::HashGraph &, const vector<string> &, const int);
  void align();
};

#endif // ALIGNER_H_
