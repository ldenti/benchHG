#ifndef CSCORER_H_
#define CSCORER_H_

#include <string>
#include <vector>

#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include "aligner.hpp"
#include "graph.hpp"

using namespace std;

class CScorer {
public:
  bcf_srs_t *vcf;
  bcf_hdr_t *hdr;
  bcf1_t *rec;
  string seq_name;
  int start;
  int stop;
  map<string, float> results;
  CScorer(const string &, const string &, const int, const int);
  void compute(const vector<Alignment> &, const Graph &);

private:
  bool is_subpath(const vector<int> &, const vector<int> &);
  bool check_ins(const vector<int> &, const vector<int> &);
  bool check_del(const vector<int> &, const vector<int> &, const vector<int> &);
};

#endif // CSCORER_H_
