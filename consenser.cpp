#include "consenser.hpp"

Consenser::Consenser(char *vcf_path, int hap) {
  args = (args_t *)calloc(1, sizeof(args_t));
  args->haplotype = hap;
  args->fname = vcf_path;
  init_data();
}

// FIXME: is this better?
// Consenser::~Consenser(void) { destroy_data(); }

string Consenser::build(char *region, char *seq) {
  init_region(args, region);

  args->fa_length = strlen(seq);
  args->fa_src_pos = strlen(seq);

  // determine if uppercase or lowercase is used in this fasta file
  args->fa_case = toupper(seq[0]) == seq[0] ? 1 : 0;

  kputs(seq, &args->fa_buf);

  bcf1_t **rec_ptr = NULL;
  while (args->rid >= 0 && (rec_ptr = next_vcf_line(args))) {
    bcf1_t *rec = *rec_ptr;
    if (strcmp(args->chr, bcf_hdr_id2name(args->hdr, rec->rid)) != 0 || (args->fa_end_pos && rec->pos > args->fa_end_pos))
      break;

    // clang-format off
    /**
       ==57670== Conditional jump or move depends on uninitialised value(s)
       ==57670==    at 0x484ED28: strlen (in /usr/libexec/valgrind/vgpreload_memcheck-amd64-linux.so)
       ==57670==    by 0x4F23DB3: std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >::basic_string(char const*, std::allocator<char> co$
       ==57670==    by 0x159297: main (main.cpp:232)
       ==57670==  Uninitialised value was created by a heap allocation
       ==57670==    at 0x484DCD3: realloc (in /usr/libexec/valgrind/vgpreload_memcheck-amd64-linux.so)
       ==57670==    by 0x154FC2: ks_resize (kstring.h:160)
       ==57670==    by 0x154FC2: ks_resize (kstring.h:155)
       ==57670==    by 0x154FC2: apply_variant (consensus.c:1002)
       ==57670==    by 0x158DA2: get_consensus(char*, char*, args_t*) (main.cpp:164)
       ==57670==    by 0x15920A: main (main.cpp:227)
    **/
    // clang-format on
    apply_variant(args, rec);
  }

  // cerr << "Applied " << args->napplied << " variants" << endl;
  string hap(args->fa_buf.s, args->fa_buf.l);
  transform(hap.begin(), hap.end(), hap.begin(), ::toupper);
  return hap;
}

void Consenser::init_data() {
  args->files = bcf_sr_init();
  args->files->require_index = 1;
  if (!bcf_sr_add_reader(args->files, args->fname))
    exit(1); // error("Failed to read from %s: %s\n", !strcmp("-", args->fname)
             // ? "standard input" : args->fname,
             // bcf_sr_strerror(args->files->errnum));
  args->hdr = args->files->readers[0].header;
  args->smpl = smpl_ilist_init(args->hdr, NULL, 0, SMPL_NONE | SMPL_VERBOSE);
  args->isample = args->smpl->idx[0];
  rbuf_init(&args->vcf_rbuf, 100);
  args->vcf_buf = (bcf1_t **)calloc(args->vcf_rbuf.m, sizeof(bcf1_t *));
  args->rid = -1;
}

void Consenser::destroy_data() {
  if (args->smpl)
    smpl_ilist_destroy(args->smpl);
  bcf_sr_destroy(args->files);
  int i;
  for (i = 0; i < args->vcf_rbuf.m; i++)
    if (args->vcf_buf[i])
      bcf_destroy1(args->vcf_buf[i]);
  free(args->vcf_buf);
  free(args->fa_buf.s);
  free(args->chr);
  free(args);
}
