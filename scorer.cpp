#include "scorer.hpp"

Scorer::Scorer(const string &vcf_path, const string &_seq_name,
               const int _start, const int _stop, const float _score) {
  vcf = bcf_sr_init();
  vcf->require_index = 1;
  if (!bcf_sr_add_reader(vcf, vcf_path.c_str()))
    cerr << "Failed to read from " << vcf_path << endl;
  hdr = vcf->readers[0].header;
  seq_name = _seq_name;
  bcf_sr_seek(vcf, seq_name.c_str(), _start);
  rec = bcf_init();
  start = _start;
  stop = _stop;
  score = _score;
}

void Scorer::compute() {
  while (bcf_sr_next_line(vcf)) {
    rec = vcf->readers[0].buffer[0];
    if (seq_name.compare(bcf_hdr_id2name(hdr, rec->rid)) != 0 ||
        rec->pos > stop)
      break;
    bcf_unpack(rec, BCF_UN_ALL);
    char *idx(rec->d.id);
    results[idx] = score;
  }
  // TODO: move this somewhere else
  bcf_sr_destroy(vcf);
  // bcf_destroy(rec); // CHECKME: we HAVE to do this (accordingly to valgrind
  // too). But why do we have a double free or corruption (!prev) with this?
}
