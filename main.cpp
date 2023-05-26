#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "bdsg/hash_graph.hpp"
#include "constructor.hpp"
#include "htslib/faidx.h"
// #include "io/save_handle_graph.hpp"
#include "gfa.hpp"
#include "kseq.h"
#include "region.hpp"

#include "align.hpp"
#include "graphLoad.hpp"

KSEQ_INIT(gzFile, gzread)

using namespace std;

extern "C" {
int main_consensus(int argc, char *argv[]);
}

int main(int argc, char *argv[]) {
  char *fa_path = argv[1];
  char *vcf_path = argv[2];
  char *region = argv[3];

  // Create or load the .fai
  faidx_t *fai = fai_load3_format(fa_path, NULL, NULL, FAI_CREATE, FAI_FASTA);

  // Extract region from .fa
  hts_pos_t seq_len;
  char *region_seq = fai_fetch64(fai, region, &seq_len);

  // Store region to file
  char *region_fa_path = (char *)malloc(strlen(region) + 4);
  strcpy(region_fa_path, region);
  strcat(region_fa_path, ".fa");
  ofstream outfa(region_fa_path);
  outfa << ">" << region << "\n" << region_seq << endl;
  outfa.close();
  free(region_seq);
  outfa.close();

  // Extract subhaplotypes and store to files
  char *hap1_path = (char *)malloc(strlen(region) + 6);
  strcpy(hap1_path, region);
  strcat(hap1_path, ".1.fa");
  char *cmd1[] = {"consensus",    "-H",     "1",  "-f",
                  region_fa_path, vcf_path, "-o", hap1_path};
  main_consensus(8, cmd1);

  char *hap2_path = (char *)malloc(strlen(region) + 6);
  strcpy(hap2_path, region);
  strcat(hap2_path, ".2.fa");
  char *cmd2[] = {"consensus",    "-H",     "2",  "-f",
                  region_fa_path, vcf_path, "-o", hap2_path};
  optind = 1;
  main_consensus(8, cmd2);

  // Read subhaplotypes
  gzFile fp = gzopen(hap1_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  char *hap1;
  while ((l = kseq_read(seq)) >= 0) {
    hap1 = (char *)malloc(l);
    strcpy(hap1, seq->seq.s);
  }
  gzclose(fp);

  fp = gzopen(hap2_path, "r");
  seq = kseq_init(fp);
  char *hap2;
  while ((l = kseq_read(seq)) >= 0) {
    hap2 = (char *)malloc(l);
    strcpy(hap2, seq->seq.s);
  }
  gzclose(fp);
  kseq_destroy(seq);

  free(region_fa_path);
  free(hap1_path);
  free(hap2_path);

  // Build the graph
  // vg construct -N -a -r {input.fa} -v {input.cvcf} -R {wildcards.region} >
  // {output.vg} # -S -i -f ?
  bdsg::HashGraph constructed;

  string seq_name;
  int64_t start_pos = -1, stop_pos = -1;
  vg::parse_region(region, seq_name, start_pos, stop_pos);

  vg::Constructor constructor;
  constructor.alt_paths = true;
  constructor.allowed_vcf_names.insert(seq_name);
  constructor.allowed_vcf_regions[seq_name] =
      make_pair(start_pos - 1, stop_pos);

  vector<string> fasta_filenames;
  fasta_filenames.push_back(fa_path);
  vector<string> vcf_filenames;
  vcf_filenames.push_back(vcf_path);
  vector<string> insertion_filenames;

  constructor.max_node_size = 32;
  constructor.construct_graph(fasta_filenames, vcf_filenames,
                              insertion_filenames, &constructed);
  stringstream ss;
  constructed.serialize(ss);
  // cout << ss.str();

  // vg convert --gfa-out {output.vg} > {output.gfa}
  // const vg::PathHandleGraph *graph_to_write =
  //     dynamic_cast<const vg::PathHandleGraph *>(&constructed);
  // set<string> rgfa_paths;
  // bool rgfa_pline = false;
  // bool wline = true;
  // vg::graph_to_gfa(graph_to_write, std::cout, rgfa_paths, rgfa_pline, wline);

  // Align to the graph
  psgl::graphLoader g;
  g.loadFromHG(constructed);

  vector<string> haps;
  haps.push_back(hap1);
  haps.push_back(hap2);
  vector<psgl::ContigInfo> qmetadata;
  qmetadata.push_back(psgl::ContigInfo{"hap1", strlen(hap1)});
  qmetadata.push_back(psgl::ContigInfo{"hap2", strlen(hap2)});

  vector<psgl::BestScoreInfo> bestScoreVector;
  psgl::Parameters parameters = {"", "", "", "", 4, 1, 1, 1, 1};
  alignToDAGLocal(haps, g.diCharGraph, parameters, bestScoreVector);

  for (auto &e : bestScoreVector) {
    cerr << e.qryId << " " << e.cigar << " ";
    std::vector<uint> path;
    path.push_back(g.diCharGraph.originalVertexId[e.refColumnStart].first);
    cerr << g.diCharGraph.originalVertexId[e.refColumnStart].first;
    for (const int32_t c : e.refColumns) {
      if (c >= e.refColumnStart && c <= e.refColumnEnd) {
        int32_t n = g.diCharGraph.originalVertexId[c].first;
        if (n != path.back()) {
          cerr << " " << n;
          path.push_back(n);
        }
      }
    }
    cerr << endl;
  }

  free(hap1);
  free(hap2);
  return 0;
}
