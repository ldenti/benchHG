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
#include "graph.hpp"

using namespace std;

bool is_subpath(const vector<int> &A, const vector<int> &B) {
  // FIXME: is this the best way to do this?
  string sA = "";
  for (const int &a : A)
    sA += to_string(a);
  string sB = "";
  for (const int &b : B)
    sB += to_string(b);
  return sA.find(sB) != sA.npos;
}

bool check_ins(const vector<int> &true_path, const vector<int> &alt_path) {
  return is_subpath(true_path, alt_path);
}

bool check_del(const vector<int> &true_path, const vector<int> &pre_nodes,
               const vector<int> &post_nodes) {
  for (const int &pre_n : pre_nodes)
    for (const int &post_n : post_nodes)
      if (is_subpath(true_path, {pre_n, post_n}))
        return true; // CHECKME: is this correct? One is enough?
  return false;
}

int main(int argc, char *argv[]) {
  spdlog::info("Welcome to spdlog!");
  parse_arguments(argc, argv);

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

  // Extract subhaplotypes
  string hap1, hap2;
  Consenser c1(tvcf_path, 1);
  hap1 = c1.build(region, region_seq);
  c1.destroy_data();
  Consenser c2(tvcf_path, 2);
  hap2 = c2.build(region, region_seq);
  c2.destroy_data();

  // Build the graph
  Graph graph(opt::fa_path, opt::cvcf_path, opt::region);
  graph.build();
  graph.analyze();

  // Align to the graph
  vector<string> haps;
  haps.push_back(hap1);
  haps.push_back(hap2);
  Aligner al(graph.hg, haps, 2);
  al.align();

  // Iterating over truth
  bcf_srs_t *vcf = bcf_sr_init();
  vcf->require_index = 1;
  bcf_hdr_t *hdr;
  if (!bcf_sr_add_reader(vcf, opt::tvcf_path.c_str()))
    cerr << "Failed to read from " << opt::tvcf_path << endl;
  hdr = vcf->readers[0].header;
  bcf_sr_seek(vcf, seq_name.c_str(), start_pos);
  bcf1_t *rec = bcf_init();
  int score = 0;
  while (bcf_sr_next_line(vcf)) {
    rec = vcf->readers[0].buffer[0];
    if (rec->pos > stop_pos)
      break;
    bcf_unpack(rec, BCF_UN_ALL);
    char *idx(rec->d.id);
    bcf_fmt_t *fmt = bcf_get_fmt(hdr, rec, "GT");
    if (!fmt)
      return 1;
    uint8_t *ptr = fmt->p; // First sample
    int a1 = bcf_gt_allele(ptr[0]);
    int a2 = bcf_gt_allele(ptr[1]);
    // FIXME: assuming we have always two alignments
    score = 0;
    if (a1 == 1 && a2 == 0)
      score = al.alignments[0].score;
    else if (a1 == 0 && a2 == 1)
      score = al.alignments[1].score;
    else if (a1 == 0 && a2 == 0)
      score = -1;
    else if (a1 == a2)
      score = max(al.alignments[0].score, al.alignments[1].score);
    else
      score = (al.alignments[0].score + al.alignments[1].score) / 2;
    cerr << "T " << idx << " " << score << endl;
    // TRUTHS[idx] = score
  }
  bcf_sr_destroy(vcf);

  // Assign scores to calls
  vcf = bcf_sr_init();
  vcf->require_index = 1;
  if (!bcf_sr_add_reader(vcf, opt::cvcf_path.c_str()))
    cerr << "Failed to read from " << opt::cvcf_path << endl;
  hdr = vcf->readers[0].header;
  bcf_sr_seek(vcf, seq_name.c_str(), start_pos);
  while (bcf_sr_next_line(vcf)) {
    rec = vcf->readers[0].buffer[0];
    if (rec->pos > stop_pos)
      break;
    bcf_unpack(rec, BCF_UN_ALL);

    // ID
    char *idx(rec->d.id);
    string seq_name = bcf_hdr_id2name(hdr, rec->rid);
    int pos = rec->pos;
    string refall = rec->d.allele[0];
    transform(refall.begin(), refall.end(), refall.begin(), ::toupper);
    vector<pair<vector<int>, bool>> alts;
    // retrieve path by matching the sequences
    assert(rec->n_allele <= 3); // FIXME: assuming diploid
    for (int i = 1; i < rec->n_allele; ++i) {
      string altall = rec->d.allele[i];
      transform(altall.begin(), altall.end(), altall.begin(), ::toupper);
      bool is_ins = refall.size() < altall.size();
      for (const pair<vector<int>, string> &path : graph.alt_paths) {
        if (is_ins && path.second.compare(altall.substr(1)) == 0) {
          alts.push_back(make_pair(path.first, is_ins));
          break;
        } else if (!is_ins && path.second.compare(refall.substr(1)) == 0) {
          alts.push_back(make_pair(path.first, is_ins));
          break;
        }
      }
    }
    // string quality = record->qual == "nan" ? "." : record->qual;

    // We assume that we must have found the paths since graph is built from
    // same VCF. If not, issue could be in graph construction
    score = 0;
    if (rec->n_allele == 2) {
      // Just one alternate allele
      assert(alts.size() == 1);
      // we assign the score if the allele is covered by at least one path
      if ((alts[0].second &&
           (check_ins(al.alignments[0].path, alts[0].first) ||
            check_ins(al.alignments[1].path, alts[0].first))) ||
          (!alts[0].second &&
           (check_del(al.alignments[0].path,
                      graph.in_edges[alts[0].first.front()],
                      graph.out_edges[alts[0].first.back()]) ||
            check_del(al.alignments[1].path,
                      graph.in_edges[alts[0].first.front()],
                      graph.out_edges[alts[0].first.back()]))))
        score = (al.alignments[0].score + al.alignments[1].score) / 2;
    } else {
      // Two alternate alleles
      assert(alts.size() == 2);
      int score1 = 0, score2 = 0;
      bool is_ins = alts[0].second;
      bool is_covered_1_by_1 =
          is_ins ? check_ins(al.alignments[0].path, alts[0].first)
                 : check_del(al.alignments[0].path,
                             graph.in_edges[alts[0].first.front()],
                             graph.out_edges[alts[0].first.back()]);

      bool is_covered_1_by_2 =
          is_ins ? check_ins(al.alignments[1].path, alts[0].first)
                 : check_del(al.alignments[1].path,
                             graph.in_edges[alts[0].first.front()],
                             graph.out_edges[alts[0].first.back()]);
      bool is_covered_2_by_1 =
          is_ins ? check_ins(al.alignments[0].path, alts[1].first)
                 : check_del(al.alignments[0].path,
                             graph.in_edges[alts[1].first.front()],
                             graph.out_edges[alts[1].first.back()]);
      bool is_covered_2_by_2 =
          is_ins ? check_ins(al.alignments[1].path, alts[1].first)
                 : check_del(al.alignments[1].path,
                             graph.in_edges[alts[1].first.front()],
                             graph.out_edges[alts[1].first.back()]);

      // CHECKME: assuming that the same path cannot cover both alleles
      if ((is_covered_1_by_1 || is_covered_1_by_2) &&
          (is_covered_2_by_1 || is_covered_2_by_2)) {
        // arbitrary
        score1 = al.alignments[0].score / 2;
        score2 = al.alignments[1].score / 2;
      } else {
        if (is_covered_1_by_1 || is_covered_1_by_2) {
          if (is_covered_1_by_1)
            score1 = al.alignments[0].score / 2;
          else
            score1 = al.alignments[1].score / 2;
        } else if (is_covered_2_by_1 || is_covered_2_by_2) {
          if (is_covered_2_by_1)
            score2 = al.alignments[0].score / 2;
          else
            score2 = al.alignments[1].score / 2;
        }
      }
      score = score1 + score2;
    }
    cerr << "C " << idx << " " << score << endl;
  }

  bcf_sr_destroy(vcf);
  // bcf_destroy(rec); // CHECKME: why we don't need this?

  free(region);
  free(region_seq);
  free(tvcf_path);
  return 0;
}
