#ifndef TSCORER_H_
#define TSCORER_H_

#include <iostream>
#include <string>
#include <vector>

#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include <aligner.hpp>

using namespace std;

class TScorer {
public:
  bcf_srs_t *vcf;
  bcf_hdr_t *hdr;
  bcf1_t *rec;
  int stop;
  map<string, float> results;
  TScorer(const string &, const string &, const int, const int);
  void compute(const vector<Alignment> &);
};

#endif // TSCORER_H_
