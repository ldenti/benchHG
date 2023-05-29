#include "aligner.hpp"

Aligner::Aligner(const bdsg::HashGraph &g, const vector<string> &_queries,
                 int threads = 1) {
  graph.loadFromHG(g);
  queries = _queries;
  parameters = {"", "", "", "", threads, 1, 1, 1, 1};
}

void Aligner::align() {
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
  alignToDAGLocal(queries, graph.diCharGraph, parameters, bestScoreVector);

  for (auto &e : bestScoreVector) {
    vector<int> path;
    path.push_back(graph.diCharGraph.originalVertexId[e.refColumnStart].first);
    for (const int32_t c : e.refColumns) {
      if (c >= e.refColumnStart && c <= e.refColumnEnd) {
        int32_t n = graph.diCharGraph.originalVertexId[c].first;
        if (n != path.back())
          path.push_back(n);
      }
    }
    Alignment a = {e.qryId, (int)queries[e.qryId].size(), 0, e.cigar, path};
    a.set_score();
    alignments.push_back(a);
  }
}
// for (const auto &a : alignments) {
//   cerr << a.id << " " << a.cigar << " " << a.score << " |";
//   for (const auto &v : a.path)
//     cerr << " " << v;
//   cerr << endl;
// }}
