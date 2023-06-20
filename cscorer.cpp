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
  start = start_pos;
  stop = stop_pos;
}

void CScorer::compute(const vector<Alignment> &alignments, const Graph &graph) {
  float score = 0.0;
  while (bcf_sr_next_line(vcf)) {
    rec = vcf->readers[0].buffer[0];
    if (seq_name.compare(bcf_hdr_id2name(hdr, rec->rid)) != 0 ||
        rec->pos > stop)
      break;
    bcf_unpack(rec, BCF_UN_ALL);
    char *idx(rec->d.id);
    string refall = rec->d.allele[0];
    transform(refall.begin(), refall.end(), refall.begin(), ::toupper);

    vector<tuple<vector<int>, vector<int>, bool>> alts;
    // retrieve path by matching the sequences
    assert(rec->n_allele <= 3); // FIXME: assuming diploid
    for (int i = 1; i < rec->n_allele; ++i) {
      string altall = rec->d.allele[i];
      transform(altall.begin(), altall.end(), altall.begin(), ::toupper);
      bool is_ins = refall.size() < altall.size();
      for (const tuple<string, vector<int>, string> &path : graph.alt_paths) {
        if (is_ins) {
          if (get<0>(path).back() == '0')
            // we look for alternate paths
            continue;
          if (altall.find(get<2>(path)) != string::npos) {
            alts.push_back(make_tuple(get<1>(path), vector<int>(), is_ins));
            break;
          }
        } else {
          if (get<0>(path).back() != '0')
            // we look for reference path
            continue;
          if (refall.find(get<2>(path)) != string::npos) {
            // Check if we have the alternate as path
            vector<int> apath;
            for (const tuple<string, vector<int>, string> &_path :
                 graph.alt_paths) {
              string _path_name = get<0>(_path);
              if (get<0>(path).compare(_path_name) == 0)
                // same path
                continue;
              _path_name.pop_back();
              _path_name += '0';
              if (get<0>(path).compare(_path_name) == 0) {
                apath = get<1>(_path);
                break;
              }
            }
            alts.push_back(make_tuple(get<1>(path), apath, is_ins));
            break;
          }
        }
      }
    }
    // string quality = record->qual == "nan" ? "." : record->qual;

    // We assume that we must have found the paths since graph is built from
    // same VCF. If not, issue could be in graph construction
    score = 0.0;
    if (rec->n_allele == 2) {
      // Just one alternate allele
      if (alts.size() != 1)
        cerr << alts.size() << " " << seq_name << ":" << start << "-" << stop
             << endl;
      ;
      assert(alts.size() == 1);
      // we assign the score if the allele is covered by at least one path
      bool found = false;
      if (get<2>(alts[0])) {
        // INS
        found = check_ins(alignments.at(0).path, get<0>(alts[0])) ||
                check_ins(alignments.at(1).path, get<0>(alts[0]));
      } else {
        // DEL
        if (get<1>(alts[0]).empty())
          // no alternate path in the graph
          found = check_del(alignments.at(0).path,
                            graph.in_edges.at(get<0>(alts.at(0)).front()),
                            graph.out_edges.at(get<0>(alts.at(0)).back())) ||
                  check_del(alignments.at(1).path,
                            graph.in_edges.at(get<0>(alts.at(0)).front()),
                            graph.out_edges.at(get<0>(alts.at(0)).back()));
        else
          found = check_ins(alignments.at(0).path, get<1>(alts[0])) ||
                  check_ins(alignments.at(1).path, get<1>(alts[0]));
      }
      if (found)
        score = (alignments.at(0).score + alignments.at(1).score) / 2.0;
      else
        score = -2;
    } else {
      if (alts.size() != 2)
        cerr << seq_name << ":" << start << "-" << stop << " " << alts.size()
             << endl;
      // Two alternate alleles
      assert(alts.size() == 2);

      bool is_covered_1_by_1 = false, is_covered_1_by_2 = false,
           is_covered_2_by_1 = false, is_covered_2_by_2 = false;
      if (get<2>(alts.at(0))) {
        is_covered_1_by_1 =
            check_ins(alignments.at(0).path, get<0>(alts.at(0)));
        is_covered_1_by_2 =
            check_ins(alignments.at(1).path, get<0>(alts.at(0)));
      } else {
        if (get<1>(alts.at(0)).empty()) {
          is_covered_1_by_1 =
              check_del(alignments.at(0).path,
                        graph.in_edges.at(get<0>(alts.at(0)).front()),
                        graph.out_edges.at(get<0>(alts.at(0)).back()));
          is_covered_1_by_2 =
              check_del(alignments.at(1).path,
                        graph.in_edges.at(get<0>(alts.at(0)).front()),
                        graph.out_edges.at(get<0>(alts.at(0)).back()));
        } else {
          is_covered_1_by_1 =
              check_ins(alignments.at(0).path, get<1>(alts.at(0)));
          is_covered_1_by_2 =
              check_ins(alignments.at(1).path, get<1>(alts.at(0)));
        }
      }
      if (get<2>(alts.at(1))) {
        is_covered_2_by_1 =
            check_ins(alignments.at(0).path, get<0>(alts.at(1)));
        is_covered_2_by_2 =
            check_ins(alignments.at(1).path, get<0>(alts.at(1)));
      } else {
        if (get<1>(alts.at(1)).empty()) {
          is_covered_2_by_1 =
              check_del(alignments.at(0).path,
                        graph.in_edges.at(get<0>(alts.at(1)).front()),
                        graph.out_edges.at(get<0>(alts.at(1)).back()));
          is_covered_2_by_2 =
              check_del(alignments.at(1).path,
                        graph.in_edges.at(get<0>(alts.at(1)).front()),
                        graph.out_edges.at(get<0>(alts.at(1)).back()));
        } else {
          is_covered_2_by_1 =
              check_ins(alignments.at(0).path, get<1>(alts.at(1)));
          is_covered_2_by_2 =
              check_ins(alignments.at(1).path, get<1>(alts.at(1)));
        }
      }
      // CHECKME: assuming that the same path cannot cover both alleles
      float score1 = 0.0, score2 = 0.0;
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
    // if (score != 1)
    //   cerr << seq_name << ":" << start << "-" << stop << " " << score << endl;
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
