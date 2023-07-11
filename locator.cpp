#include "locator.hpp"

Locator::Locator(const int _k, const int _w, const int _W) {
  k = _k;
  w = _w;
  W = _W;
}

map<string, interval_tree_t<int>> Locator::build_tree(const string &bed_path,
                                                      const int w = 0) {
  map<string, interval_tree_t<int>> trees;
  string line;
  ifstream bed_f(bed_path);
  if (bed_f.is_open()) {
    while (getline(bed_f, line)) {
      uint pos = 0;
      string delimiter = "\t";
      pos = line.find(delimiter);
      string seq_name = line.substr(0, pos);
      line.erase(0, pos + delimiter.length());
      pos = line.find(delimiter);
      int start_pos = stoi(line.substr(0, pos));
      line.erase(0, pos + delimiter.length());
      int stop_pos = stoi(line);

      if (trees.find(seq_name) == trees.end())
        trees[seq_name] = interval_tree_t<int>();
      trees[seq_name].insert({start_pos - w - 1, stop_pos + w});
    }
    bed_f.close();
  }
  for (auto it = trees.begin(); it != trees.end(); ++it)
    it->second.deoverlap();
  return trees;
}

void Locator::add_conf(const string &bed_path) {
  conf_trees = build_tree(bed_path);
}

void Locator::add_trf(const string &bed_path) {
  trf_trees = build_tree(bed_path, 0 /*w / 2*/);
}

void Locator::parse_truth(faidx_t *fai, const string &vcf_path,
                          const string &ovcf_path) {
  true_trees = parse(fai, vcf_path, ovcf_path);
}

void Locator::parse_call(faidx_t *fai, const string &vcf_path,
                         const string &ovcf_path) {
  call_trees = parse(fai, vcf_path, ovcf_path);
}

map<string, interval_tree_t<int>>
Locator::parse(faidx_t *fai, const string &vcf_path, const string &ovcf_path) {

  map<string, interval_tree_t<int>> trees;

  htsFile *vcf = bcf_open(vcf_path.c_str(), "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  bcf1_t *vcf_record = bcf_init();
  bcf_info_t *l_info;

  htsFile *ovcf = bcf_open(ovcf_path.c_str(), "wz");
  bcf_hdr_write(ovcf, vcf_header);

  string seq_name;
  int key_idx;
  int pos, stop;
  int istart = -1;
  int istop = -1;
  string last_seq_name;
  int last_trf_s = -1, last_trf_e = -1;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_INFO);
    seq_name = bcf_hdr_id2name(vcf_header, vcf_record->rid);
    pos = vcf_record->pos;
    stop = pos + 1;

    bcf_fmt_t *fmt = bcf_get_fmt(vcf_header, vcf_record, "GT");
    uint8_t *ptr = fmt->p + fmt->size * 0;
    if (bcf_gt_allele(ptr[0]) == 0 && bcf_gt_allele(ptr[0]) == 0)
      continue;

    l_info = bcf_get_info(vcf_header, vcf_record, "SVLEN");
    key_idx = l_info->key;
    const char *key = vcf_header->id[BCF_DT_ID][key_idx].key;
    int *v = NULL;
    int ndst = 0;
    bcf_get_info_int32(vcf_header, vcf_record, key, &v, &ndst);
    int l = *v;
    if (l < 0)
      stop = pos - l;

    // Variant is from pos to stop

    if (conf_trees.find(seq_name) != conf_trees.end() &&
        conf_trees[seq_name].overlap_find({pos, stop + 1}) ==
            end(conf_trees[seq_name]))
      // variant outside confidence regions
      continue;

    bcf_write1(ovcf, vcf_header, vcf_record);

    if (last_seq_name != seq_name || istart == -1) {
      if (istart != -1)
        // Let's add last interval on previous chromosome
        trees[last_seq_name].insert({istart, istop + 1});
      last_seq_name = seq_name;
      istart = pos;
      istop = stop;
      last_trf_s = -1;
      last_trf_e = -1;
    }
    int trf_s = -1, trf_e = -1;
    if (trf_trees.find(seq_name) != trf_trees.end()) {
      auto overlap = trf_trees[seq_name].overlap_find({pos, stop});
      // CHECKME: do we always have one overlap?
      if (overlap != end(trf_trees[seq_name])) {
        trf_s = overlap->low();
        trf_e = overlap->high();
      }
      if (last_trf_s != -1 && last_trf_e != -1 && last_trf_s == trf_s &&
          last_trf_e == trf_e) {
        if (pos - (istop + 1) > W) {
          // we are too far anyway, so let's add the region
          trees[seq_name].insert({istart, istop + 1});
          last_seq_name = seq_name;
          istart = pos;
          istop = stop;
        } else {
          istop = max(istop, stop);
        }
        continue;
      } else {
        // current trf region becomes last analyzed
        last_trf_s = trf_s;
        last_trf_e = trf_e;
      }
    }
    if (pos - (istop + 1) > w) {
      // cerr << "I: " << seq_name << " " << istart << " " << istop + 1 << endl;
      // if (seq_name.compare("chr10") == 0)
      //     cerr << seq_name << ":" << istart << "-" << istop << endl;
      trees[seq_name].insert({istart, istop + 1});
      last_seq_name = seq_name;
      istart = pos;
      istop = stop;
      // cerr << "(new) Setting istart to " << istart << endl;
    }
    istop = max(istop, stop);
    // cerr << "(iter) Setting istop to " << istop << endl;
  }
  if (istart != -1) {
    // cerr << "I: " << seq_name << " " << istart << " " << istop + 1 << endl;
    trees[last_seq_name].insert({istart, istop + 1});
  }

  for (auto it = trees.begin(); it != trees.end(); ++it) {
    it->second.deoverlap();

    interval_tree_t<int> copy = it->second;
    // for (const auto &i : it->second) {
    //   copy.insert({i.low(), i.high()});
    // }

    for (const auto &i : copy) {
      // cerr << "Parsing " << it->first << " " << i.low() << " " << i.high() <<
      // endl; cerr << it->first << " " << copy.size() << " vs " <<
      // it->second.size() << endl;
      auto overlap = it->second.overlap_find(i);
      assert(overlap != end(it->second));
      int s = overlap->low(), e = overlap->high();
      assert(s == i.low() && e == i.high());
      it->second.erase(overlap);
      string region =
          it->first + ":" + to_string(s - w + 1) + "-" + to_string(e + w);
      hts_pos_t seq_len;
      string region_seq = fai_fetch64(fai, region.c_str(), &seq_len);
      transform(region_seq.begin(), region_seq.end(), region_seq.begin(),
                ::toupper);
      assert(region_seq.size() > 0);
      pair<int, int> uniques = get_unique_kmers(region_seq);
      // cerr << "Inserting " << s - w + 1 + uniques.first << " " << s - w + 1 +
      // uniques.second << endl;
      it->second.insert(
          {s - w + 1 + uniques.first, s - w + 1 + uniques.second});
    }
    it->second.deoverlap();
  }

  bcf_hdr_destroy(vcf_header);
  bcf_close(vcf);
  bcf_close(ovcf);
  bcf_destroy(vcf_record);

  tbx_index_build2(ovcf_path.c_str(), (ovcf_path + ".tbi").c_str(), 14,
                   &tbx_conf_vcf);

  return trees;
}

void Locator::intersect(ofstream &of) {
  for (auto it = true_trees.begin(); it != true_trees.end(); ++it)
    for (const auto &i : it->second)
      of << it->first << "\t" << i.low() << "\t" << i.high() << "\t"
         << "T"
         << "\t" << i.high() - i.low() + 1 << endl;
  for (auto it = call_trees.begin(); it != call_trees.end(); ++it)
    for (const auto &i : it->second)
      of << it->first << "\t" << i.low() << "\t" << i.high() << "\t"
         << "C"
         << "\t" << i.high() - i.low() + 1 << endl;

  for (auto it = call_trees.begin(); it != call_trees.end(); ++it) {
    for (const auto &i : it->second) {
      if (true_trees[it->first].overlap_find(i) == end(true_trees[it->first]))
        continue;
      int s = i.low(), e = i.high();
      true_trees[it->first].overlap_find_all(i, [&s, &e](auto iter) {
        s = min(s, iter->low());
        e = max(e, iter->high());
        return true;
      });
      regions[it->first].insert({s, e});
    }
    regions[it->first].deoverlap();
  }
}

vector<string> Locator::get_regions(ofstream &of) const {
  vector<string> Rs;
  for (auto it = regions.begin(); it != regions.end(); ++it)
    for (const auto &i : it->second) {
      of << it->first << "\t" << i.low() << "\t" << i.high() << "\t"
         << "M"
         << "\t" << i.high() - i.low() + 1 << endl;
      Rs.push_back(it->first + ":" + to_string(i.low()) + "-" +
                   to_string(i.high()));
    }
  return Rs;
}

pair<int, int> Locator::get_unique_kmers(const string &seq) {
  map<string, int> KMERS;
  for (int i = 0; i < seq.size() - k + 1; ++i) {
    string kmer(seq, i, k);
    // if (KMERS.find(kmer) == KMERS.end())
    //   KMERS[kmer] = 1;
    // else
    ++KMERS[kmer];
  }
  int first = seq.size(), last = -1;
  for (int i = 0; i < seq.size() - k + 1; ++i) {
    string kmer(seq, i, k);
    if (KMERS[kmer] == 1) {
      if (i < first)
        first = i;
      if (i > last)
        last = i;
    }
  }
  first = first == seq.size() || first > w / 2 ? 0 : first;
  last = last == -1 || last < seq.size() - w / 2 ? seq.size() : last + k;
  return make_pair(first, last);
}
