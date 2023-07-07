#ifndef GRAPH_H_
#define GRAPH_H_

#include <string>

#include <bdsg/hash_graph.hpp>
#include <constructor.hpp>
#include <gfa.hpp>
#include <region.hpp>

using namespace std;

class Graph {
public:
  bdsg::HashGraph hg;
  string seq_name;
  int64_t start_pos;
  int64_t stop_pos;
  vector<string> fasta_filenames;
  vector<string> vcf_filenames;
  vector<string> insertion_filenames;
  vector<int> ref_path;
  vector<tuple<string, vector<int>, string>> alt_paths;
  map<int, vector<int>>
      in_edges; // for each source of alt_paths, store in-nodes
  map<int, vector<int>>
      out_edges; // for each sink of alt_paths, store out-nodes

public:
  Graph(const string &, const string &, const string &);
  void build();
  void analyze();
  void to_gfa() const;
};
#endif // GRAPH_H_
