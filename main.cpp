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

  fai_destroy(fai);
  regions_file.close();
  if (opt::regions_only)
    return 0;

  spdlog::info("Cleaning regions..");

  map<string, pair<string, float>> truth_results;
  map<string, pair<string, float>> call_results;
  map<string, float> regions_results;
  vector<string> clean_regions;

  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);
  bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
  bcf_sr_add_reader(sr, filtered_tvcf_path.c_str());
  bcf_sr_add_reader(sr, filtered_cvcf_path.c_str());
  bcf_hdr_t *thdr = bcf_sr_get_header(sr, 0);
  bcf_hdr_t *chdr = bcf_sr_get_header(sr, 1);

  filtered_cvcf_path = opt::out + ".predictions.clean.vcf.gz";

  htsFile *cvcf_o = bcf_open(filtered_cvcf_path.c_str(), "wz");
  bcf_hdr_write(cvcf_o, chdr);

  bcf1_t *trec, *crec;
  int ndst = 8;
  char *svtype = (char *)malloc(ndst);

  for (const string &region : regions) {
    regions_results[region] = -1.0;

    string seq_name;
    int64_t start_pos = -1, stop_pos = -1;
    vg::parse_region(region, seq_name, start_pos, stop_pos);

    int true_vtypes = 0;
    bool tend = false, cend = false;
    bool is_clean_region = false;

    bcf_sr_seek(sr, seq_name.c_str(), start_pos);
    while (bcf_sr_next_line(sr)) {
      if (tend && cend)
        break;
      // truth
      if (bcf_sr_has_line(sr, 0)) {
        trec = bcf_sr_get_line(sr, 0);
        if (strcmp(seq_name.c_str(), bcf_hdr_id2name(thdr, trec->rid)) == 0 &&
            trec->pos <= stop_pos) {
          // FIXME: assuming SVTYPE to be in the INFO
          // svtype_info = bcf_get_info(thdr, trec, "SVTYPE");
          // svtype_idx = svtype_info->key;
          // const char *key = thdr->id[BCF_DT_ID][svtype_idx].key;
          bcf_get_info_string(thdr, trec, "SVTYPE", &svtype, &ndst);
          true_vtypes |= strcmp(svtype, "INS") == 0 ? 1 : 2;
          truth_results[trec->d.id] = make_pair(region, -1.0);
        } else
          tend = true;
      }
      // predictions
      if (bcf_sr_has_line(sr, 1)) {
        crec = bcf_sr_get_line(sr, 1);
        if (strcmp(seq_name.c_str(), bcf_hdr_id2name(chdr, crec->rid)) == 0 &&
            crec->pos <= stop_pos) {
          bcf_get_info_string(chdr, crec, "SVTYPE", &svtype, &ndst);
          if (((strcmp(svtype, "INS") == 0 ? 1 : 2) & true_vtypes)) {
            bcf_write1(cvcf_o, chdr, crec);
            is_clean_region = true;
          }
          call_results[crec->d.id] = make_pair(region, -1.0);
        } else
          cend = true;
      }
    }
    if (is_clean_region)
      // push the region if we have a type match
      clean_regions.push_back(region);
  }
  free(svtype);
  bcf_sr_destroy(sr);
  bcf_close(cvcf_o);
  tbx_index_build2(filtered_cvcf_path.c_str(),
                   (filtered_cvcf_path + ".tbi").c_str(), 14, &tbx_conf_vcf);

  spdlog::info("Preallocating data..");
  vector<faidx_t *> FAI(opt::threads);
  vector<char *> REGION(opt::threads);
  // vector<char *> REGION_SEQ(opt::threads);

  for (uint i = 0; i < opt::threads; ++i) {
    FAI[i] = fai_load3_format(opt::fa_path.c_str(), NULL, NULL, FAI_CREATE,
                              FAI_FASTA);
    REGION[i] = (char *)malloc(512);
  }

  char *tvcf_path = (char *)malloc(filtered_tvcf_path.size() + 1);
  strcpy(tvcf_path, filtered_tvcf_path.c_str());
  tvcf_path[filtered_tvcf_path.size()] = '\0';

  spdlog::info("Starting analysis of {} regions..", clean_regions.size());
  omp_set_dynamic(0);
  omp_set_num_threads(opt::threads);

#pragma omp parallel for
  for (uint i = 0; i < clean_regions.size(); ++i) {
    strcpy(REGION[omp_get_thread_num()], clean_regions[i].c_str());
    REGION[omp_get_thread_num()][clean_regions[i].size()] = '\0';

    string seq_name;
    int64_t start_pos = -1, stop_pos = -1;
    vg::parse_region(REGION[omp_get_thread_num()], seq_name, start_pos,
                     stop_pos);

#ifdef PDEBUG
    cerr << omp_get_thread_num() << " " << clean_regions[i] << " ("
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
      spdlog::warn("Skipping {} due to error in consensus", clean_regions[i]);
      continue;
    }
    if (!c1.has_alts && !c2.has_alts) {
      spdlog::warn("Skipping {} since no alt alleles in truth",
                   clean_regions[i]);
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
    graph.to_gfa();
    for (const int &n : graph.ref_path)
      cerr << n << " ";
    cerr << endl;
#endif

    float gscore1 = -1, gscore2 = -1;
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

#ifdef PDEBUG
    for (const auto &a : al.alignments) {
      cerr << a.id << " " << a.cigar << " ";
      for (const auto &i : a.path)
        cerr << " " << i;
      cerr << " " << a.score << endl;
    }
#endif

    assert(gscore1 != -1 && gscore2 != -1);

    float score1 = !is_subpath(graph.ref_path, path1) ? gscore1 : 0;
    float score2 = !is_subpath(graph.ref_path, path2) ? gscore2 : 0;
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
      truth_results.at(res.first) = make_pair(clean_regions[i], res.second);
    for (const pair<string, float> res : cs.results)
      call_results.at(res.first) = make_pair(clean_regions[i], res.second);
    regions_results.at(clean_regions[i]) = score;

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
  for (const auto &r : truth_results)
    ofile << "T " << r.first << " " << r.second.first << " " << r.second.second
          << endl;
  spdlog::info("Outputting call scores..");
  for (const auto &r : call_results)
    ofile << "C " << r.first << " " << r.second.first << " " << r.second.second
          << endl;

  return 0;
}
