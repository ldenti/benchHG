#include "tscorer.hpp"

TScorer::TScorer(const string &vcf_path, const string &seq_name,
                 const int start_pos, const int stop_pos) {
  vcf = bcf_sr_init();
  vcf->require_index = 1;
  if (!bcf_sr_add_reader(vcf, vcf_path.c_str()))
    cerr << "Failed to read from " << vcf_path << endl;
  hdr = vcf->readers[0].header;
  bcf_sr_seek(vcf, seq_name.c_str(), start_pos);
  rec = bcf_init();
  stop = stop_pos;
}

void TScorer::compute(const vector<Alignment> &alignments) {
  int score = 0;
  while (bcf_sr_next_line(vcf)) {
    rec = vcf->readers[0].buffer[0];
    if (rec->pos > stop)
      break;
    bcf_unpack(rec, BCF_UN_ALL);
    char *idx(rec->d.id);
    bcf_fmt_t *fmt = bcf_get_fmt(hdr, rec, "GT");
    if (!fmt) {
      cerr << "No GT format in VCF. Halting.." << endl;
      exit(EXIT_FAILURE);
    }
    uint8_t *ptr = fmt->p; // First sample
    int a1 = bcf_gt_allele(ptr[0]);
    int a2 = bcf_gt_allele(ptr[1]);
    // FIXME: assuming we have always two alignments
    score = 0;
    if (a1 == 1 && a2 == 0)
      score = alignments[0].score;
    else if (a1 == 0 && a2 == 1)
      score = alignments[1].score;
    else if (a1 == 0 && a2 == 0)
      score = -1;
    else if (a1 == a2)
      score = max(alignments[0].score, alignments[1].score);
    else
      score = (alignments[0].score + alignments[1].score) / 2;
    cerr << "T " << idx << " " << score << endl;
    // TRUTHS[idx] = score
  }
  // TODO: move this somewhere else
  bcf_sr_destroy(vcf);
  // bcf_destroy(rec); // CHECKME: we HAVE to do this (accordingly to valgrind
  // too). But why do we have a double free or corruption (!prev) with this?
}
