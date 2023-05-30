#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

#include <htslib/faidx.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include <bdsg/hash_graph.hpp>
#include <constructor.hpp>
// #include "gfa.hpp"
#include <region.hpp> // vg

#include <spdlog/spdlog.h>

#include "aligner.hpp"
#include "argument_parser.hpp"
#include "consenser.hpp"
#include "cscorer.hpp"
#include "graph.hpp"
#include "tscorer.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  parse_arguments(argc, argv);

  vector<map<string, map<string, float>>> t_results(opt::threads);
  vector<map<string, map<string, float>>> c_results(opt::threads);

  // spdlog::info("Initializing..");

  vector<string> regions;
  string line;
  ifstream region_f(opt::region);
  if (region_f.is_open()) {
    while (getline(region_f, line)) {
      regions.push_back(line);
      string seq_name;
      int64_t start_pos = -1, stop_pos = -1;
      vg::parse_region(line, seq_name, start_pos, stop_pos);
      for (int i = 0; i < opt::threads; ++i) {
        t_results[i][seq_name] = map<string, float>();
        c_results[i][seq_name] = map<string, float>();
      }
    }
    region_f.close();
  }

  // TODO: split more intelligently?
  omp_set_dynamic(0);
  omp_set_num_threads(opt::threads);
#pragma omp parallel for
  for (uint i = 0; i < regions.size(); ++i) {
    char *region = (char *)malloc(regions[i].size() + 1);
    strcpy(region, regions[i].c_str());
    region[regions[i].size()] = '\0';

    string seq_name;
    int64_t start_pos = -1, stop_pos = -1;
    vg::parse_region(region, seq_name, start_pos, stop_pos);
    cerr << omp_get_thread_num() << " " << regions[i] << " ("
         << stop_pos - start_pos + 1 << ")" << endl;

    char *tvcf_path = (char *)malloc(opt::tvcf_path.size() + 1);
    strcpy(tvcf_path, opt::tvcf_path.c_str());
    tvcf_path[opt::tvcf_path.size()] = '\0';

    // Extract region from .fai
    faidx_t *fai = fai_load3_format(opt::fa_path.c_str(), NULL, NULL,
                                    FAI_CREATE, FAI_FASTA);
    hts_pos_t seq_len;
    char *region_seq = fai_fetch64(fai, region, &seq_len);
    fai_destroy(fai);

    // spdlog::info("Building consensus..");
    string hap1, hap2;
    Consenser c1(tvcf_path, 1);
    hap1 = c1.build(region, region_seq);
    c1.destroy_data();
    Consenser c2(tvcf_path, 2);
    hap2 = c2.build(region, region_seq);
    c2.destroy_data();

    // spdlog::info("Building graph..");
    Graph graph(opt::fa_path, opt::cvcf_path, string(region));
    graph.build();
    graph.analyze();

    // spdlog::info("Aligning to graph..");
    vector<string> haps;
    haps.push_back(hap1);
    haps.push_back(hap2);
    Aligner al(graph.hg, haps, 2);
    al.align();

    // Iterating over truth
    TScorer ts(opt::tvcf_path, seq_name, start_pos, stop_pos);
    ts.compute(al.alignments);

    CScorer cs(opt::cvcf_path, seq_name, start_pos, stop_pos);
    cs.compute(al.alignments, graph);

    for (const pair<string, float> res : ts.results)
      t_results[omp_get_thread_num()][seq_name][res.first] = res.second;
    for (const pair<string, float> res : cs.results)
      c_results[omp_get_thread_num()][seq_name][res.first] = res.second;

    free(region);
    free(region_seq);
    free(tvcf_path);
  }

  // OUTPUT
  // 1. Linearize
  map<string, map<string, float>> truth_results;
  map<string, map<string, float>> call_results;
  for (int i = 0; i < opt::threads; ++i) {
    for (auto it1 = t_results[i].begin(); it1 != t_results[i].end(); ++it1) {
      for (auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
        truth_results[it1->first][it2->first] = it2->second;
      }
    }
  }

  for (int i = 0; i < opt::threads; ++i) {
    for (auto it1 = c_results[i].begin(); it1 != c_results[i].end(); ++it1) {
      for (auto it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
        call_results[it1->first][it2->first] = it2->second;
      }
    }
  }

  // 2. Iterate over truth and assign scores
  htsFile *vcf = bcf_open(opt::tvcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  bcf1_t *vcf_record = bcf_init();
  string seq_name, idx;
  float score;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    seq_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
    idx = vcf_record->d.id;
    score = -1;
    if (truth_results[seq_name].find(idx) != truth_results[seq_name].end())
      score = truth_results[seq_name][idx];
    cout << "T " << idx << " " << score << endl;
  }
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);

  // 3. Iterate over call and assign scores
  vcf = bcf_open(opt::cvcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  vcf_record = bcf_init();
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    seq_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
    idx = vcf_record->d.id;
    score = -1;
    if (call_results[seq_name].find(idx) != call_results[seq_name].end())
      score = call_results[seq_name][idx];
    cout << "C " << idx << " " << score << endl;
  }
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);
  bcf_destroy(vcf_record);

  return 0;
}
