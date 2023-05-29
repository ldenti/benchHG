#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>

#include <bcftools/filter.h>
#include <bcftools/rbuf.h>
#include <bcftools/regidx.h>
extern "C" {
#include <bcftools/smpl_ilist.h>
}

#include "bdsg/hash_graph.hpp"
#include "constructor.hpp"
#include "gfa.hpp"
#include "region.hpp"

#include "align.hpp"
#include "graphLoad.hpp"

using namespace std;

struct Alignment {
  int id;
  int l;
  int il;
  string cigar;
  vector<int> path;
  float score;

  void set_score() {
    map<char, int> OPs;
    OPs['='] = 0;
    OPs['X'] = 0;
    OPs['I'] = 0;
    OPs['D'] = 0;
    regex word_regex("([0-9]+[=XID])");
    auto words_begin = sregex_iterator(cigar.begin(), cigar.end(), word_regex);
    auto words_end = std::sregex_iterator();
    for (sregex_iterator i = words_begin; i != words_end; ++i) {
      smatch match = *i;
      string match_str = match.str();
      char op = match_str.back();
      match_str.pop_back();
      int l = stoi(match_str);
      OPs[op] += l;
    }
    il = OPs['X'] + OPs['='] + OPs['I'];
    score = 1 - (OPs['I'] + OPs['D'] + OPs['X'] + abs(l - il)) / l; // CHECKME
  }
};

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

void get_consensus(char *region, char *seq, args_t *args) {
  init_region(args, region);

  args->fa_length = strlen(seq);
  args->fa_src_pos = strlen(seq);

  // determine if uppercase or lowercase is used in this fasta file
  args->fa_case = toupper(seq[0]) == seq[0] ? 1 : 0;

  kputs(seq, &args->fa_buf);

  bcf1_t **rec_ptr = NULL;
  while (args->rid >= 0 && (rec_ptr = next_vcf_line(args))) {
    bcf1_t *rec = *rec_ptr;
    if (args->fa_end_pos && rec->pos > args->fa_end_pos)
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
}

static void init_data(args_t *args) {
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

static void destroy_data(args_t *args) {
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
}

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
  char *fa_path = argv[1];
  char *tvcf_path = argv[2];
  char *cvcf_path = argv[3];
  char *region = argv[4];

  // Create or load the .fai
  faidx_t *fai = fai_load3_format(fa_path, NULL, NULL, FAI_CREATE, FAI_FASTA);

  // Extract region from .fa
  hts_pos_t seq_len;
  char *region_seq = fai_fetch64(fai, region, &seq_len);
  fai_destroy(fai);

  // Extract subhaplotypes
  args_t *args = (args_t *)calloc(1, sizeof(args_t));
  args->haplotype = 1;
  args->fname = tvcf_path;
  init_data(args);
  get_consensus(region, region_seq, args);
  cerr << "Applied " << args->napplied << " variants" << endl;
  // char *hap1 = (char *)malloc(args->fa_buf.l + 1);
  // strcpy(hap1, args->fa_buf.s);
  // hap1[args->fa_buf.l] = '\0';
  string hap1(args->fa_buf.s);
  destroy_data(args);
  free(args);

  args = (args_t *)calloc(1, sizeof(args_t));
  args->haplotype = 2;
  args->fname = tvcf_path;
  init_data(args);
  get_consensus(region, region_seq, args);
  cerr << "Applied " << args->napplied << " variants" << endl;
  // char *hap2 = (char *)malloc(args->fa_buf.l + 1);
  // strcpy(hap2, args->fa_buf.s);
  // hap2[args->fa_buf.l] = '\0';
  string hap2(args->fa_buf.s);
  destroy_data(args);
  free(args);

  free(region_seq);

  // Build the graph
  // vg construct -N -a -r {input.fa} -v {input.cvcf} -R
  // {wildcards.region} >
  // {output.vg} # -S -i -f ?
  bdsg::HashGraph graph;

  string seq_name;
  int64_t start_pos = -1, stop_pos = -1;
  vg::parse_region(region, seq_name, start_pos, stop_pos);

  vg::Constructor constructor;
  constructor.alt_paths = true;
  constructor.allowed_vcf_names.insert(seq_name);
  constructor.allowed_vcf_regions[seq_name] =
      make_pair(start_pos - 1, stop_pos);

  vector<string> fasta_filenames;
  fasta_filenames.push_back(fa_path);
  vector<string> vcf_filenames;
  vcf_filenames.push_back(cvcf_path);
  vector<string> insertion_filenames;

  constructor.max_node_size = 32;
  // clang-format off
  /**
     ==57670==    at 0x4848899: malloc (in /usr/libexec/valgrind/vgpreload_memcheck-amd64-linux.so)
     ==57670==    by 0x492FF9B: ??? (in /usr/lib/x86_64-linux-gnu/libhts.so.1.13+ds)
     ==57670==    by 0x49302B9: bgzf_hopen (in /usr/lib/x86_64-linux-gnu/libhts.so.1.13+ds)
     ==57670==    by 0x494CA9C: hts_hopen (in /usr/lib/x86_64-linux-gnu/libhts.so.1.13+ds)
     ==57670==    by 0x494CCE7: hts_open_format (in /usr/lib/x86_64-linux-gnu/libhts.so.1.13+ds)
     ==57670==    by 0x286AA7: Tabix::Tabix(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&) (tabix.cpp:46)
     ==57670==    by 0x18D415: openTabix (Variant.h:124)
     ==57670==    by 0x18D415: open (Variant.h:109)
     ==57670==    by 0x18D415: vg::Constructor::construct_graph(std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::$
     ==57670==    by 0x1C27DB: operator() (std_function.h:590)
     ==57670==    by 0x1C27DB: vg::io::load_proto_to_graph(handlegraph::MutablePathMutableHandleGraph*, std::function<void (std::function<void (vg::Graph&)> const$
     ==57670==    by 0x17E072: vg::Constructor::construct_graph(std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::$
     ==57670==    by 0x15954E: main (main.cpp:261)
  **/
  // clang-format on

  constructor.construct_graph(fasta_filenames, vcf_filenames,
                              insertion_filenames, &graph);

  // Dump .gfa to stdout (optional)
  stringstream graph_ss;
  graph.serialize(graph_ss);
  // vg convert --gfa-out {output.vg} > {output.gfa}
  const vg::PathHandleGraph *graph_to_write =
      dynamic_cast<const vg::PathHandleGraph *>(&graph);
  set<string> rgfa_paths;
  bool rgfa_pline = false;
  bool wline = true;
  vg::graph_to_gfa(graph_to_write, std::cout, rgfa_paths, rgfa_pline, wline);
  // ---------------------

  // Align to the graph
  psgl::graphLoader g;
  g.loadFromHG(graph);

  vector<string> haps;
  haps.push_back(hap1);
  haps.push_back(hap2);

  vector<psgl::ContigInfo> qmetadata;
  qmetadata.push_back(psgl::ContigInfo{"hap1", (int)hap1.size()});
  qmetadata.push_back(psgl::ContigInfo{"hap2", (int)hap2.size()});

  vector<psgl::BestScoreInfo> bestScoreVector;
  psgl::Parameters parameters = {"", "", "", "", 1, 1, 1, 1, 1};

  // clang-format off
  /**
     ==57670==    at 0x48487A9: malloc (in /usr/libexec/valgrind/vgpreload_memcheck-amd64-linux.so)
     ==57670==    by 0x50F58AC: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
     ==57670==    by 0x5106ACC: ??? (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
     ==57670==    by 0x50FCA10: GOMP_parallel (in /usr/lib/x86_64-linux-gnu/libgomp.so.1.0.0)
     ==57670==    by 0x15750A: psgl::alignToDAGLocal_Phase1_scalar(std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, st$
     ==57670==    by 0x157DC3: psgl::alignToDAGLocal(std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<s$
     ==57670==    by 0x15984C: main (main.cpp:287)
   **/
  // clang-format on

  alignToDAGLocal(haps, g.diCharGraph, parameters, bestScoreVector);

  vector<Alignment> alignments;
  for (auto &e : bestScoreVector) {
    vector<int> path;
    path.push_back(g.diCharGraph.originalVertexId[e.refColumnStart].first);
    for (const int32_t c : e.refColumns) {
      if (c >= e.refColumnStart && c <= e.refColumnEnd) {
        int32_t n = g.diCharGraph.originalVertexId[c].first;
        if (n != path.back())
          path.push_back(n);
      }
    }
    Alignment a = {e.qryId, (int)haps[e.qryId].size(), 0, e.cigar, path};
    a.set_score();
    alignments.push_back(a);
  }

  // for (const auto &a : alignments) {
  //   cerr << a.id << " " << a.cigar << " " << a.score << " |";
  //   for (const auto &v : a.path)
  //     cerr << " " << v;
  //   cerr << endl;
  // }

  // get info from graph
  vector<pair<vector<int>, string>> alt_paths;
  map<int, vector<int>>
      in_edges; // for each source of alt_paths, store in-nodes
  map<int, vector<int>>
      out_edges; // for each sink of alt_paths, store out-nodes
  vector<string> path_names;
  graph.for_each_path_handle([&](const bdsg::path_handle_t &p) {
    path_names.emplace_back(graph.get_path_name(p));
  });
  for (int i = 0; i < graph.get_path_count(); ++i) {
    if (path_names[i].front() != '_')
      continue;
    bdsg::path_handle_t ph = graph.get_path_handle(path_names[i]);
    vector<int> path;
    string path_seq;
    for (bdsg::handle_t handle : graph.scan_path(ph)) {
      path.push_back(graph.get_id(handle));
      path_seq += graph.get_sequence(handle);
    }
    graph.follow_edges(graph.get_handle(path.front()), true,
                       [&](const bdsg::handle_t &n) {
                         in_edges[path.front()].push_back(graph.get_id(n));
                       });
    graph.follow_edges(graph.get_handle(path.back()), false,
                       [&](const bdsg::handle_t &n) {
                         out_edges[path.back()].push_back(graph.get_id(n));
                       });
    transform(path_seq.begin(), path_seq.end(), path_seq.begin(), ::toupper);
    alt_paths.push_back(make_pair(path, path_seq));
  }

  // clang-format off
  // in_edges = {path[0]: [] for path, _ in paths}
  //   out_edges = {path[-1]: [] for path, _ in paths}
  //   for line in open(gfa_path):
  //       if line.startswith("L"):
  //           _, id1, _, id2, _, _ = line.strip("\n").split("\t")
  //           if id1 in out_edges:
  //               out_edges[id1].append(id2)
  //           if id2 in in_edges:
  //               in_edges[id2].append(id1)
  // clang-format on
  // Iterating over truth
  bcf_srs_t *vcf = bcf_sr_init();
  vcf->require_index = 1;
  bcf_hdr_t *hdr;
  if (!bcf_sr_add_reader(vcf, tvcf_path))
    cerr << "Failed to read from " << tvcf_path << endl;
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
  bcf_sr_destroy(vcf);

  // Assign scores to calls
  vcf = bcf_sr_init();
  vcf->require_index = 1;
  if (!bcf_sr_add_reader(vcf, cvcf_path))
    cerr << "Failed to read from " << cvcf_path << endl;
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
      for (const pair<vector<int>, string> &path : alt_paths) {
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

    if (rec->n_allele == 2) {
      // Just one alternate allele
      assert(alts.size() == 1);
      // we assign the score if the allele is covered by at least one path
      if ((alts[0].second && (check_ins(alignments[0].path, alts[0].first) ||
                              check_ins(alignments[1].path, alts[0].first))) ||
          (!alts[0].second &&
           (check_del(alignments[0].path, in_edges[alts[0].first.front()],
                      out_edges[alts[0].first.back()]) ||
            check_del(alignments[1].path, in_edges[alts[0].first.front()],
                      out_edges[alts[0].first.back()]))))
        score = (alignments[0].score + alignments[1].score) / 2;
    } else {
      // Two alternate alleles
      assert(alts.size() == 2);
      int score1 = 0, score2 = 0;
      bool is_ins = alts[0].second;
      bool is_covered_1_by_1 =
          is_ins
              ? check_ins(alignments[0].path, alts[0].first)
              : check_del(alignments[0].path, in_edges[alts[0].first.front()],
                          out_edges[alts[0].first.back()]);

      bool is_covered_1_by_2 =
          is_ins
              ? check_ins(alignments[1].path, alts[0].first)
              : check_del(alignments[1].path, in_edges[alts[0].first.front()],
                          out_edges[alts[0].first.back()]);
      bool is_covered_2_by_1 =
          is_ins
              ? check_ins(alignments[0].path, alts[1].first)
              : check_del(alignments[0].path, in_edges[alts[1].first.front()],
                          out_edges[alts[1].first.back()]);
      bool is_covered_2_by_2 =
          is_ins
              ? check_ins(alignments[1].path, alts[1].first)
              : check_del(alignments[1].path, in_edges[alts[1].first.front()],
                          out_edges[alts[1].first.back()]);

      // CHECKME: assuming that the same path cannot cover both alleles
      if ((is_covered_1_by_1 || is_covered_1_by_2) &&
          (is_covered_2_by_1 || is_covered_2_by_2)) {
        // arbitrary
        score1 = alignments[0].score / 2;
        score2 = alignments[1].score / 2;
      } else {
        if (is_covered_1_by_1 || is_covered_1_by_2) {
          if (is_covered_1_by_1)
            score1 = alignments[0].score / 2;
          else
            score1 = alignments[1].score / 2;
        } else if (is_covered_2_by_1 || is_covered_2_by_2) {
          if (is_covered_2_by_1)
            score2 = alignments[0].score / 2;
          else
            score2 = alignments[1].score / 2;
        }
      }
      score = score1 + score2;
    }

    cerr << "C " << idx << " " << score << endl;
  }

  bcf_sr_destroy(vcf);
  // bcf_destroy(rec); // CHECKME: why we don't need this?

  cerr << "END" << endl;
  return 0;
}
