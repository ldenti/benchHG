#include <algorithm>
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
#include <htslib/vcf.h>
#include <region.hpp> // vg
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "aligner.hpp"
#include "argument_parser.hpp"
#include "consenser.hpp"
#include "graph.hpp"
#include "locator.hpp"
#include "scorer.hpp"

#include "ksw2.h"

using namespace std;

string seq_align(const char *tseq, const char *qseq, int sc_mch, int sc_mis,
                 int gapo, int gape) {
  int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a,
                    b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  int tl = strlen(tseq), ql = strlen(qseq);
  uint8_t *ts, *qs, c[256];
  ksw_extz_t ez;

  memset(&ez, 0, sizeof(ksw_extz_t));
  memset(c, 4, 256);
  c['A'] = c['a'] = 0;
  c['C'] = c['c'] = 1;
  c['G'] = c['g'] = 2;
  c['T'] = c['t'] = 3; // build the encoding table
  ts = (uint8_t *)malloc(tl);
  qs = (uint8_t *)malloc(ql);
  for (i = 0; i < tl; ++i)
    ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
  for (i = 0; i < ql; ++i)
    qs[i] = c[(uint8_t)qseq[i]];
  ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
  string cigar = "";
  for (i = 0; i < ez.n_cigar; ++i) {
    cigar += to_string(ez.cigar[i] >> 4) + "MID"[ez.cigar[i] & 0xf];
    // printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
  }
  // putchar('\n');
  free(ez.cigar);
  free(ts);
  free(qs);
  return cigar;
}

int main(int argc, char *argv[]) {
  spdlog::set_default_logger(spdlog::stderr_color_st("stderr"));
  parse_arguments(argc, argv);

  spdlog::info("Initializing..");
  faidx_t *fai =
      fai_load3_format(opt::fa_path.c_str(), NULL, NULL, FAI_CREATE, FAI_FASTA);

  string filtered_tvcf_path;
  string filtered_cvcf_path;
  vector<string> regions;
  ofstream regions_file(opt::out + ".regions.bed");
  if (opt::region.compare("") == 0) {
    spdlog::info("Creating trees..");
    Locator l(opt::k, opt::w, opt::W);
    l.add_conf(opt::conf);
    l.add_trf(opt::trf);
    spdlog::info("Locating regions of interest from truth..");
    filtered_tvcf_path = opt::out + ".truth.vcf.gz";
    l.parse_truth(fai, opt::tvcf_path, filtered_tvcf_path);
    spdlog::info("Locating regions of interest from call..");
    filtered_cvcf_path = opt::out + ".predictions.vcf.gz";
    l.parse_call(fai, opt::cvcf_path, filtered_cvcf_path);
    spdlog::info("Intersecting regions..");
    l.intersect(regions_file);
    regions = l.get_regions(regions_file);
  } else {
    regions.push_back(opt::region);
    filtered_tvcf_path = opt::tvcf_path;
    filtered_cvcf_path = opt::cvcf_path;
  }
  fai_destroy(fai);
  regions_file.close();

  if (opt::regions_only)
    return 0;

  spdlog::info("Preallocating data..");
  map<string, pair<string, float>> truth_results;
  map<string, pair<string, float>> call_results;
  map<string, float> regions_results;
  vector<faidx_t *> FAI(opt::threads);
  vector<char *> REGION(opt::threads);
  // vector<char *> REGION_SEQ(opt::threads);

  for (uint i = 0; i < regions.size(); ++i) {
    regions_results[regions[i]] = -1.0;
  }

  htsFile *vcf = bcf_open(filtered_tvcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  bcf1_t *vcf_record = bcf_init();
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    truth_results[vcf_record->d.id] = make_pair("", -1.0);
  }
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);

  vcf = bcf_open(filtered_cvcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  vcf_record = bcf_init();
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    call_results[vcf_record->d.id] = make_pair("", -1.0);
  }
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);

  for (uint i = 0; i < opt::threads; ++i) {
    FAI[i] = fai_load3_format(opt::fa_path.c_str(), NULL, NULL, FAI_CREATE,
                              FAI_FASTA);
    REGION[i] = (char *)malloc(512);
  }
  char *tvcf_path = (char *)malloc(filtered_tvcf_path.size() + 1);
  strcpy(tvcf_path, filtered_tvcf_path.c_str());
  tvcf_path[filtered_tvcf_path.size()] = '\0';

  spdlog::info("Starting analysis of {} regions..", regions.size());
  omp_set_dynamic(0);
  omp_set_num_threads(opt::threads);

#pragma omp parallel for
  for (uint i = 0; i < regions.size(); ++i) {
    strcpy(REGION[omp_get_thread_num()], regions[i].c_str());
    REGION[omp_get_thread_num()][regions[i].size()] = '\0';

    string seq_name;
    int64_t start_pos = -1, stop_pos = -1;
    vg::parse_region(REGION[omp_get_thread_num()], seq_name, start_pos,
                     stop_pos);

    // if (stop_pos - start_pos + 1 > 75000)
    //   continue;

#ifdef PDEBUG
    cerr << omp_get_thread_num() << " " << regions[i] << " ("
         << stop_pos - start_pos + 1 << ")" << endl;
#endif

    // Extract region from .fai
    hts_pos_t seq_len;
    char *region_seq = fai_fetch64(FAI[omp_get_thread_num()],
                                   REGION[omp_get_thread_num()], &seq_len);
    // spdlog::info("Building consensus..");
    string hap1, hap2;
    Consenser c1(tvcf_path, 1);
    hap1 = c1.build(REGION[omp_get_thread_num()], region_seq);
    c1.destroy_data();
    Consenser c2(tvcf_path, 2);
    hap2 = c2.build(REGION[omp_get_thread_num()], region_seq);
    c2.destroy_data();

    if (hap1 == "" or hap2 == "") {
      spdlog::warn("Skipping {} due to error in consensus", regions[i]);
      continue;
    }
    if (!c1.has_alts && !c2.has_alts) {
      spdlog::warn("Skipping {} since no alt alleles in truth", regions[i]);
      continue;
    }

#ifdef PDEBUG
    cerr << region_seq << endl;
    cerr << hap1 << endl;
    cerr << hap2 << endl;
#endif

    // spdlog::info("Building graph..");
    Graph graph(opt::fa_path, filtered_cvcf_path,
                string(REGION[omp_get_thread_num()]));
    graph.build();
    graph.analyze();
#ifdef PDEBUG
    // graph.to_gfa();
    for (const int &n : graph.ref_path)
      cerr << n << " ";
    cerr << endl;
#endif

    float gscore1 = -1, gscore2 = -1, rscore1 = -1, rscore2 = -1;
    vector<string> haps;
    if (c1.has_alts)
      haps.push_back(hap1);
    if (c2.has_alts)
      haps.push_back(hap2);
    Aligner al(graph.hg, haps, 2);
    al.align();
    vector<int> path1;
    vector<int> path2;
    if (c1.has_alts && c2.has_alts) {
      gscore1 = al.alignments[0].score;
      path1 = al.alignments[0].path;
      gscore2 = al.alignments[1].score;
      path2 = al.alignments[1].path;
    } else if (c1.has_alts) {
      gscore1 = al.alignments[0].score;
      path1 = al.alignments[0].path;
      gscore2 = gscore1;
      path2 = path1;
    } else {
      gscore2 = al.alignments[0].score;
      path2 = al.alignments[0].path;
      gscore1 = gscore2;
      path1 = path2;
    }

    //     Alignment refal1, refal2;
    //     if (c1.has_alts) {
    //       string cigar = seq_align(region_seq, hap1.c_str(), 1, -2, 2, 1);
    //       refal1 = {-1, hap1.size(), 0, cigar, vector<int>(), 0.0};
    //       refal1.set_score();
    //       rscore1 = refal1.score;
    // #ifdef PDEBUG
    //       cerr << hap1.size() << " " << refal1.cigar << " " << rscore1 <<
    //       endl;
    // #endif
    //     } else {
    // #ifdef PDEBUG
    //       cerr << "No ref alignment 1" << endl;
    // #endif
    //     }
    //     if (c2.has_alts) {
    //       string cigar = seq_align(region_seq, hap2.c_str(), 1, -2, 2, 1);
    //       refal2 = {-1, hap2.size(), 0, cigar, vector<int>(), 0.0};
    //       refal2.set_score();
    //       rscore2 = refal2.score;
    // #ifdef PDEBUG
    //       cerr << hap2.size() << " " << refal2.cigar << " " << rscore2 <<
    //       endl;
    // #endif
    //     } else {
    // #ifdef PDEBUG
    //       cerr << "No ref alignment 2" << endl;
    // #endif
    //     }
    //     rscore1 = rscore1 == -1 ? rscore2 : rscore1;
    //     rscore2 = rscore2 == -1 ? rscore1 : rscore2;

#ifdef PDEBUG
    for (const auto &a : al.alignments) {
      cerr << a.id << " " << a.cigar << " ";
      for (const auto &i : a.path)
        cerr << " " << i;
      cerr << " " << a.score << endl;
    }
#endif

    // assert(rscore1 != 1 && rscore2 != 1);
    assert(gscore1 != -1 &&
           gscore2 != -1 /**&& rscore1 != -1 && rscore2 != -1**/);

    float score1 = !equal(path1.begin(), path1.end(), graph.ref_path.begin())
                       ? gscore1
                       : 0;
    float score2 = !equal(path2.begin(), path2.end(), graph.ref_path.begin())
                       ? gscore2
                       : 0;
    float score = (score1 + score2) / 2.0;

#ifdef PDEBUG
    cerr << "Score 1: " << score1 << endl;
    cerr << "Score 2: " << score2 << endl;
    cerr << "Score: " << score << endl;
#endif
    Scorer ts(filtered_tvcf_path, seq_name, start_pos, stop_pos, score);
    ts.compute();

    Scorer cs(filtered_cvcf_path, seq_name, start_pos, stop_pos, score);
    cs.compute();

    for (const pair<string, float> res : ts.results)
      truth_results.at(res.first) = make_pair(regions[i], res.second);
    for (const pair<string, float> res : cs.results)
      call_results.at(res.first) = make_pair(regions[i], res.second);
    regions_results.at(regions[i]) = score;

    free(region_seq);
  }
  for (uint i = 0; i < opt::threads; ++i) {
    fai_destroy(FAI[i]);
    free(REGION[i]);
  }
  free(tvcf_path);

  // OUTPUT
  ofstream ofile(opt::out + ".results.txt");
  spdlog::info("Outputting regions scores..");
  for (const auto &region : regions_results)
    ofile << "R "
          << "."
          << " " << region.first << " " << region.second << endl;

  spdlog::info("Outputting truth scores..");
  vcf = bcf_open(filtered_tvcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  vcf_record = bcf_init();
  string idx;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    idx = vcf_record->d.id;
    ofile << "T " << idx << " " << truth_results[idx].first << " "
          << truth_results[idx].second << endl;
  }
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);

  spdlog::info("Outputting call scores..");
  vcf = bcf_open(filtered_cvcf_path.c_str(), "r");
  vcf_header = bcf_hdr_read(vcf);
  vcf_record = bcf_init();
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_STR);
    idx = vcf_record->d.id;
    ofile << "C " << idx << " " << call_results[idx].first << " "
          << call_results[idx].second << endl;
  }
  ofile.close();
  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);
  bcf_destroy(vcf_record);

  return 0;
}
