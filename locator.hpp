#ifndef LOCATOR_H_
#define LOCATOR_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <interval-tree/interval_tree.hpp>

using namespace std;
using namespace lib_interval_tree;

class Locator {
public:
  int k;
  int w;
  map<string, interval_tree_t<int>> conf_trees;
  map<string, interval_tree_t<int>> trf_trees;
  map<string, interval_tree_t<int>> true_trees;
  map<string, interval_tree_t<int>> call_trees;
  map<string, interval_tree_t<int>> regions;

  Locator(const int, const int);
  void add_conf(const string &);
  void add_trf(const string &);
  void parse_truth(faidx_t *fai, const string &, const string &);
  void parse_call(faidx_t *fai, const string &, const string &);
  vector<string> get_regions(bool) const;
  void intersect(bool);

private:
  map<string, interval_tree_t<int>> build_tree(const string &, const int);
  pair<int, int> get_unique_kmers(const string &seq);
  map<string, interval_tree_t<int>> parse(faidx_t *fai, const string &,
                                          const string &);
};

#endif // LOCATOR_H_
