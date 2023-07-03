#ifndef SCORER_H_
#define SCORER_H_

#include <iostream>
#include <string>
#include <vector>

#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include <aligner.hpp>

using namespace std;

class Scorer {
public:
  bcf_srs_t *vcf;
  bcf_hdr_t *hdr;
  bcf1_t *rec;
  string seq_name;
  int start;
  int stop;
  float score1;
  float score2;
  map<string, float> results;
  Scorer(const string &, const string &, const int, const int, const float,
         const float);
  void compute();
};

#endif // TSCORER_H_
