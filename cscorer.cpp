#include "cscorer.hpp"

CScorer::CScorer(const string &vcf_path, const string &_seq_name,
                 const int start_pos, const int stop_pos) {
  vcf = bcf_sr_init();
  vcf->require_index = 1;
  if (!bcf_sr_add_reader(vcf, vcf_path.c_str()))
    cerr << "Failed to read from " << vcf_path << endl;
  hdr = vcf->readers[0].header;
  seq_name = _seq_name;
  bcf_sr_seek(vcf, seq_name.c_str(), start_pos);
  rec = bcf_init();
  stop = stop_pos;
}

void CScorer::compute(const vector<Alignment> &alignments, const Graph &graph) {
  float score = 0.0;
  while (bcf_sr_next_line(vcf)) {
    rec = vcf->readers[0].buffer[0];
    if (seq_name.compare(bcf_hdr_id2name(hdr, rec->rid)) != 0 || rec->pos > stop)
      break;
    bcf_unpack(rec, BCF_UN_ALL);
    char *idx(rec->d.id);
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
    score = 0.0;
    if (rec->n_allele == 2) {
      // Just one alternate allele
      assert(alts.size() == 1);
      // we assign the score if the allele is covered by at least one path
      if ((alts[0].second &&
           (check_ins(alignments.at(0).path, alts[0].first) ||
            check_ins(alignments.at(1).path, alts[0].first))) ||
          (!alts[0].second &&
           (check_del(alignments.at(0).path,
                      graph.in_edges.at(alts.at(0).first.front()),
                      graph.out_edges.at(alts.at(0).first.back())) ||
            check_del(alignments.at(1).path,
                      graph.in_edges.at(alts.at(0).first.front()),
                      graph.out_edges.at(alts.at(0).first.back())))))
        score = (alignments.at(0).score + alignments.at(1).score) / 2.0;
    } else {
      // Two alternate alleles
      assert(alts.size() == 2);
      int score1 = 0, score2 = 0;
      bool is_ins_1 = alts.at(0).second;
      bool is_ins_2 = alts.at(1).second;
      bool is_covered_1_by_1 =
          is_ins_1 ? check_ins(alignments.at(0).path, alts.at(0).first)
                   : check_del(alignments.at(0).path,
                               graph.in_edges.at(alts.at(0).first.front()),
                               graph.out_edges.at(alts.at(0).first.back()));
      bool is_covered_1_by_2 =
          is_ins_1 ? check_ins(alignments.at(1).path, alts.at(0).first)
                   : check_del(alignments.at(1).path,
                               graph.in_edges.at(alts.at(0).first.front()),
                               graph.out_edges.at(alts.at(0).first.back()));
      bool is_covered_2_by_1 =
          is_ins_2 ? check_ins(alignments.at(0).path, alts.at(1).first)
                   : check_del(alignments.at(0).path,
                               graph.in_edges.at(alts.at(1).first.front()),
                               graph.out_edges.at(alts.at(1).first.back()));

      bool is_covered_2_by_2 =
          is_ins_2 ? check_ins(alignments.at(1).path, alts.at(1).first)
                   : check_del(alignments.at(1).path,
                               graph.in_edges.at(alts.at(1).first.front()),
                               graph.out_edges.at(alts.at(1).first.back()));
      // CHECKME: assuming that the same path cannot cover both alleles
      if ((is_covered_1_by_1 || is_covered_1_by_2) &&
          (is_covered_2_by_1 || is_covered_2_by_2)) {
        // arbitrary
        score1 = alignments[0].score / 2.0;
        score2 = alignments[1].score / 2.0;
      } else {
        if (is_covered_1_by_1 || is_covered_1_by_2) {
          if (is_covered_1_by_1)
            score1 = alignments[0].score / 2.0;
          else
            score1 = alignments[1].score / 2.0;
        } else if (is_covered_2_by_1 || is_covered_2_by_2) {
          if (is_covered_2_by_1)
            score2 = alignments[0].score / 2.0;
          else
            score2 = alignments[1].score / 2.0;
        }
      }
      score = score1 + score2;
    }
    results[idx] = score;
  }

  // TODO: move this somewhere else
  bcf_sr_destroy(vcf);
  // bcf_destroy(rec); // CHECKME: we HAVE to do this (accordingly to valgrind
  // too). But why do we have a double free or corruption (!prev) with this?
}

bool CScorer::is_subpath(const vector<int> &A, const vector<int> &B) {
  // FIXME: is this the best way to do this?
  string sA = "";
  for (const int &a : A)
    sA += to_string(a);
  string sB = "";
  for (const int &b : B)
    sB += to_string(b);
  return sA.find(sB) != sA.npos;
}

bool CScorer::check_ins(const vector<int> &true_path,
                        const vector<int> &alt_path) {
  return is_subpath(true_path, alt_path);
}

bool CScorer::check_del(const vector<int> &true_path,
                        const vector<int> &pre_nodes,
                        const vector<int> &post_nodes) {
  for (const int &pre_n : pre_nodes)
    for (const int &post_n : post_nodes)
      if (is_subpath(true_path, {pre_n, post_n}))
        return true; // CHECKME: is this correct? One is enough?
  return false;
}
