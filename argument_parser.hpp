#ifndef ARG_PARSER_H_
#define ARG_PARSER_H_

#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

static const char *VERSION = "0.0.1";

// clang-format off
static const char *USAGE_MESSAGE =
    "\n"
    "Usage: sveval <reference.fa> <truth.vcf.gz> <call.vcf.gz>\n"
    "\n"
    "      -o, --out <prefix>                output prefix (default: OUT)\n"
    "      -k, --ksize <int>                 kmers size (default: 13)\n"
    "      -w, --window <int>                window size (default: 250)\n"
    "      -W, --window-trf <int>            window size in TRF region (default: 1000)\n"
    "      -t, --trf <trf.bed>               tandem repeats .bed file\n"
    "      -c, --conf <conf.bed>             confindence regions .bed file\n"
    "      -r, --region <chr:start-end>      analyze only this region\n"
    "      -@, --threads <int>               threads (default: 1)\n"
    "      -v, --version                     print version information\n"
    "      -h, --help                        display this help and exit\n"
    "\n";
// clang-format on

namespace opt {
static bool regions_only = false;
static int k = 13;
static int w = 250;
static int W = 1000;
static string out = "OUT";
static string trf = "";
static string conf = "";
static string region = "";
static int threads = 1;
static string fa_path;
static string tvcf_path;
static string cvcf_path;
} // namespace opt

static const char *shortopts = "o:k:w:W:t:c:r:R@:h";

static const struct option longopts[] = {
    {"out", required_argument, NULL, 'o'},
    {"ksize", required_argument, NULL, 'k'},
    {"window", required_argument, NULL, 'w'},
    {"window-trf", required_argument, NULL, 'W'},
    {"trf", required_argument, NULL, 't'},
    {"conf", required_argument, NULL, 'c'},
    {"region", required_argument, NULL, 'r'},
    {"regions", no_argument, NULL, 'R'},
    {"threads", required_argument, NULL, '@'},
    {"version", no_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};

void parse_arguments(int argc, char **argv) {
  bool die = false;
  for (char c;
       (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'o':
      arg >> opt::out;
      break;
    case 'k':
      arg >> opt::k;
      break;
    case 'w':
      arg >> opt::w;
      break;
    case 'W':
      arg >> opt::W;
      break;
    case 't':
      arg >> opt::trf;
      break;
    case 'c':
      arg >> opt::conf;
      break;
    case 'r':
      arg >> opt::region;
      break;
    case 'R':
      opt::regions_only = true;
      break;
    case '@':
      arg >> opt::threads;
      break;
    case 'v':
      cerr << "sveval " << VERSION << "\n";
      exit(EXIT_SUCCESS);
    case '?':
      die = true;
      break;
    case 'h':
      cout << USAGE_MESSAGE;
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind < 3) {
    cerr << "sveval : missing arguments\n";
    die = true;
  } else if (argc - optind > 3) {
    cerr << "sveval : too many arguments\n";
    die = true;
  }
  if (die) {
    cerr << "\n" << USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  opt::fa_path = argv[optind++];
  opt::tvcf_path = argv[optind++];
  opt::cvcf_path = argv[optind++];
}

#endif // ARG_PARSER_H_
