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

  // spdlog::info("Initializing..");

  char *region = (char *)malloc(opt::region.size() + 1);
  strcpy(region, opt::region.c_str());
  region[opt::region.size()] = '\0';

  string seq_name;
  int64_t start_pos = -1, stop_pos = -1;
  vg::parse_region(region, seq_name, start_pos, stop_pos);

  char *tvcf_path = (char *)malloc(opt::tvcf_path.size() + 1);
  strcpy(tvcf_path, opt::tvcf_path.c_str());
  tvcf_path[opt::tvcf_path.size()] = '\0';

  // Extract region from .fai
  faidx_t *fai =
      fai_load3_format(opt::fa_path.c_str(), NULL, NULL, FAI_CREATE, FAI_FASTA);
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
  Graph graph(opt::fa_path, opt::cvcf_path, opt::region);
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

  free(region);
  free(region_seq);
  free(tvcf_path);
  return 0;
}
