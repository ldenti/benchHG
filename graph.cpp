#include "graph.hpp"

Graph::Graph(const string &fa_path, const string &vcf_path,
             const string &region) {
  start_pos = -1;
  stop_pos = -1;
  vg::parse_region(region, seq_name, start_pos, stop_pos);
  fasta_filenames.push_back(fa_path);
  vcf_filenames.push_back(vcf_path);
}

void Graph::build() {
  // clang-format off
  // vg construct -N -a -r {input.fa} -v {input.cvcf} -R {wildcards.region} > {output.vg} # -S -i -f ?
  // clang-format on
  vg::Constructor constructor;
  constructor.alt_paths = true;
  constructor.allowed_vcf_names.insert(seq_name);
  constructor.allowed_vcf_regions[seq_name] =
      make_pair(start_pos - 1, stop_pos);
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

     FIXED BUT NOT MERGED YET
  **/
  // clang-format on
  constructor.construct_graph(fasta_filenames, vcf_filenames,
                              insertion_filenames, &hg);
}

void Graph::analyze() {
  vector<string> path_names;
  hg.for_each_path_handle([&](const bdsg::path_handle_t &p) {
    path_names.emplace_back(hg.get_path_name(p));
  });
  for (int i = 0; i < hg.get_path_count(); ++i) {
    if (path_names[i].front() != '_') {
      bdsg::path_handle_t ph = hg.get_path_handle(path_names[i]);
      vector<int> path;
      for (bdsg::handle_t handle : hg.scan_path(ph))
        ref_path.push_back(hg.get_id(handle));
    } else {
      bdsg::path_handle_t ph = hg.get_path_handle(path_names[i]);
      string path_name = hg.get_path_name(ph);
      vector<int> path;
      string path_seq;
      for (bdsg::handle_t handle : hg.scan_path(ph)) {
        path.push_back(hg.get_id(handle));
        path_seq += hg.get_sequence(handle);
      }
      hg.follow_edges(hg.get_handle(path.front()), true,
                      [&](const bdsg::handle_t &n) {
                        in_edges[path.front()].push_back(hg.get_id(n));
                      });
      hg.follow_edges(hg.get_handle(path.back()), false,
                      [&](const bdsg::handle_t &n) {
                        out_edges[path.back()].push_back(hg.get_id(n));
                      });
      transform(path_seq.begin(), path_seq.end(), path_seq.begin(), ::toupper);
      alt_paths.push_back(make_tuple(path_name, path, path_seq));
    }
  }
}

void Graph::to_gfa() const {
  stringstream graph_ss;
  hg.serialize(graph_ss);
  // vg convert --gfa-out {output.vg} > {output.gfa}
  const vg::PathHandleGraph *graph_to_write =
      dynamic_cast<const vg::PathHandleGraph *>(&hg);
  set<string> rgfa_paths;
  bool rgfa_pline = false;
  bool wline = true;
  vg::graph_to_gfa(graph_to_write, std::cerr, rgfa_paths, rgfa_pline, wline);
}
