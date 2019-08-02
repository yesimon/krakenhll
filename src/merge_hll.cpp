#include "hyperloglogplus.hpp"
#include "kraken_headers.hpp"
#include "readcounts.hpp"
#include "taxdb.hpp"
#include "krakenutil.hpp"
// #include "krakendb.hpp"
// #include "quickfile.hpp"
// #include "krakendb.hpp"
#include <iostream>
#include <fstream>
#include <random>

using namespace std;
using namespace kraken;

using READCOUNTS = ReadCounts<HyperLogLogPlusMinus<uint64_t>>;
string TaxDB_file;
vector<string> Input_filenames;
string Output_filename;

void parse_command_line(int argc, char **argv);
void usage(int exit_code = EX_USAGE);
unordered_map<uint32_t, READCOUNTS>read_hll(string filename);

unordered_map<uint32_t, READCOUNTS>read_hll(string filename) {
  std::ifstream file;
  file.open(filename.c_str());
  if (file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", filename.c_str());
  }

  unordered_map<uint32_t, READCOUNTS> taxon_counts;
  string line;
  int i = 0;
  while(file.good()) {
    getline(file, line);
    if (line == "") break;
    istringstream iss(line);
    uint32_t taxid;
    iss >> taxid;
    iss >> ws;
    string iss_left;
    getline(iss, iss_left);
    cout << i << "\t" << "taxid: " << taxid << endl;
    ++i;
    READCOUNTS rc(iss_left);
    int j = 0;
    taxon_counts.insert(make_pair(taxid, move(rc)));
  }
  return taxon_counts;
}

auto key_selector = [](auto pair){return pair.first;};

int main(int argc, char **argv) {
  parse_command_line(argc, argv);
  unordered_map<uint32_t, READCOUNTS> taxon_counts;
  TaxonomyDB<uint32_t> taxdb;
  unordered_set<uint32_t> taxids;

  vector<unordered_map<uint32_t, READCOUNTS>> ind_taxon_counts;
  for (auto const& filename: Input_filenames) {
    auto tc = read_hll(filename);
    ind_taxon_counts.push_back(std::move(tc));
    // cout << tc.size() << endl;
    for(auto kv : tc) {
      // cout << kv.first << endl;
      taxids.insert(kv.first);
    }
    // transform(tc.begin(), tc.end(), taxids.begin(), key_selector);
  }

  for (auto const& taxid: taxids) {
    for (auto const& tc: ind_taxon_counts) {
      auto it = tc.find(taxid);
      if (it != tc.end()) {
        // cout << taxid << endl;
        taxon_counts[taxid] += it->second;
      }
    }
  }
  if (!TaxDB_file.empty()) {
    taxdb = TaxonomyDB<uint32_t>(TaxDB_file, false);
    // taxdb.readGenomeSizes(TaxDB_file);
  } else {
    cerr << "TaxDB argument is required!" << endl;
    return 1;
  }
  managed_ostream report_output(Output_filename, true, false);

  TaxReport<uint32_t, READCOUNTS> rep = TaxReport<uint32_t, READCOUNTS>(*report_output, taxdb, taxon_counts, false);
  rep.setReportCols(vector<string>{"%", "reads", "taxReads", "kmers", "dup", "cov", "taxID", "rank", "taxName"});
  rep.printReport("kraken");
  // gettimeofday(&tv2, NULL);
  // fprintf(stderr, "Report finished in %.3f seconds.\n",
  //         get_seconds(tv1, tv2));
  return 0;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "i:d:o:")) !=
         -1) {
    switch (opt) {
    case 'd':
      TaxDB_file = optarg;
      break;
    case 'i':
      Input_filenames.push_back(optarg);
      break;
    case 'o':
      Output_filename = optarg;
      break;
    }
  }
}


void usage(int exit_code) {
  cerr << "Usage: merge_hll [options] [-i <hll file>]*" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "  -D filename      Kraken catted DB stored in a single file" << endl
       << "  -o filename      Output file for Kraken output" << endl
       << "  -r filename      Output file for Kraken report output" << endl
       << "  -a filename      TaxDB" << endl
       << "  -I filename      UID to TaxId map" << endl
       << "  -l               Memory lock DB files" << endl
       << "  -p #             Precision for unique k-mer counting, between 10 "
    "and 18"
       << endl
       << "  -t #             Number of threads" << endl
       << "  -u #             Thread work unit size (in bp)" << endl
       << "  -O               Order output matching input" << endl
       << "  -q               Quick operation" << endl
       << "  -m #             Minimum hit count (ignored w/o -q)" << endl
       << "  -C filename      Print classified sequences" << endl
       << "  -U filename      Print unclassified sequences" << endl
       << "  -H filename      Print HLL counters per-taxon" << endl
       << "  -f               Input is in FASTQ format" << endl
       << "  -c               Only include classified reads in output" << endl
       << "  -M               Preload database files" << endl
       << "  -s               Print read sequence in Kraken output" << endl
       << "  -h               Print this message" << endl
       << endl
       << "At least one FASTA or FASTQ file must be specified." << endl
       << "Kraken output is to standard output by default." << endl;
  exit(exit_code);
}
