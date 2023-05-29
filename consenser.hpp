#ifndef CONSENSER_H_
#define CONSENSER_H_

#include <algorithm>
#include <ctype.h>
#include <iostream>
#include <string>

#include <bcftools/filter.h>
#include <bcftools/rbuf.h>
#include <bcftools/regidx.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
extern "C" {
#include <bcftools/smpl_ilist.h>
}

using namespace std;

typedef struct {
  int num;            // number of ungapped blocks in this chain
  int *block_lengths; // length of the ungapped blocks in this chain
  int *ref_gaps; // length of the gaps on the reference sequence between blocks
  int *
      alt_gaps; // length of the gaps on the alternative sequence between blocks
  int ori_pos;
  int ref_last_block_ori; // start position on the reference sequence of the
                          // following ungapped block (0-based)
  int alt_last_block_ori; // start position on the alternative sequence of the
                          // following ungapped block (0-based)
} chain_t;

#define MASK_LC 1
#define MASK_UC 2
#define MASK_SKIP(x) (((x)->with != MASK_LC && (x)->with != MASK_UC) ? 1 : 0)
typedef struct {
  char *fname, with;
  regidx_t *idx;
  regitr_t *itr;
} mask_t;

typedef struct {
  kstring_t fa_buf; // buffered reference sequence
  int fa_ori_pos;   // start position of the fa_buffer (wrt original sequence)
  int fa_frz_pos; // protected position to avoid conflicting variants (last pos
                  // for SNPs/ins)
  int fa_mod_off; // position difference of fa_frz_pos in the ori and modified
                  // sequence (ins positive)
  int fa_frz_mod; // the fa_buf offset of the protected fa_frz_pos position,
                  // includes the modified sequence
  int fa_end_pos; // region's end position in the original sequence
  int fa_length;  // region's length in the original sequence (in case end_pos
                  // not provided in the FASTA header)
  int fa_case;    // output upper case or lower case: TO_UPPER|TO_LOWER
  int fa_src_pos; // last genomic coordinate read from the input fasta (0-based)
  char prev_base; // this is only to validate the REF allele in the VCF - the
                  // modified fa_buf cannot be used for inserts following
                  // deletions, see 600#issuecomment-383186778
  int prev_base_pos; // the position of prev_base
  int prev_is_insert;

  rbuf_t vcf_rbuf;
  bcf1_t **vcf_buf;
  int nvcf_buf, rid;
  char *chr, *chr_prefix;

  mask_t *mask;
  int nmask;

  int chain_id;   // chain_id, to provide a unique ID to each chain in the chain
                  // output
  chain_t *chain; // chain structure to store the sequence of ungapped blocks
                  // between the ref and alt sequences Note that the chain is
                  // re-initialised for each chromosome/seq_region

  filter_t *filter;
  char *filter_str;
  int filter_logic; // include or exclude sites which match the filters? One of
                    // FLT_INCLUDE/FLT_EXCLUDE

  bcf_srs_t *files;
  bcf_hdr_t *hdr;
  FILE *fp_out;
  FILE *fp_chain;
  char **argv;
  int argc, output_iupac, iupac_GTs, haplotype, allele, isample, napplied;
  uint8_t *iupac_bitmask, *iupac_als;
  int miupac_bitmask, miupac_als;
  char *fname, *ref_fname, *sample, *sample_fname, *output_fname, *mask_fname,
      *chain_fname, missing_allele, absent_allele;
  char mark_del, mark_ins, mark_snv;
  smpl_ilist_t *smpl;
} args_t;

extern "C" {
void init_region(args_t *args, char *line);
void apply_variant(args_t *args, bcf1_t *rec);
bcf1_t **next_vcf_line(args_t *args);
}

class Consenser {
private:
  args_t *args;
  void init_data();

public:
  Consenser(char *, int);
  // ~Consenser(void);
  string build(char *, char *);
  void destroy_data();
};
#endif // CONSENSER_H_
