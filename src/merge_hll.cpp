#include "hyperloglogplus.hpp"
#include "kraken_headers.hpp"
#include "readcounts.hpp"
#include "taxdb.hpp"
#include "krakenutil.hpp"
#include <iostream>
#include <fstream>

using namespace std;
using namespace kraken;

using READCOUNTS = ReadCounts<HyperLogLogPlusMinus<uint64_t>>;
string TaxDB_filename;
string CountsDB_filename;
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
  while(getline(file, line)) {
    istringstream iss(line);
    uint32_t taxid;
    iss >> taxid;
    iss >> ws;
    string iss_left;
    getline(iss, iss_left);
    READCOUNTS rc(iss_left);
    taxon_counts.insert(make_pair(taxid, move(rc)));
  }
  return taxon_counts;
}

int main(int argc, char **argv) {
  parse_command_line(argc, argv);

  unordered_map<uint32_t, READCOUNTS> taxon_counts;
  TaxonomyDB<uint32_t> taxdb;
  unordered_set<uint32_t> taxids;

  unordered_map<uint32_t, vector<READCOUNTS>> ind_taxon_counts;
  for (auto const& filename: Input_filenames) {
    auto tc = read_hll(filename);
    for(auto& kv : tc) {
      auto& rc = ind_taxon_counts[kv.first];
      rc.push_back(move(kv.second));
    }
  }

  for (auto const& vtc: ind_taxon_counts) {
      for (auto const& rc: vtc.second) {
        taxon_counts[vtc.first] += rc;
      }
  }
  if (!TaxDB_filename.empty()) {
    taxdb = TaxonomyDB<uint32_t>(TaxDB_filename, false);
  } else {
    cerr << "TaxDB argument is required!" << endl;
    return 1;
  }
  if (!CountsDB_filename.empty()) {
    taxdb.readGenomeSizes(CountsDB_filename);
  }

  managed_ostream report_output(Output_filename, true, false);

  TaxReport<uint32_t, READCOUNTS> rep = TaxReport<uint32_t, READCOUNTS>(*report_output, taxdb, taxon_counts, false);
  rep.setReportCols(vector<string>{"%", "reads", "taxReads", "kmers", "dup", "cov", "taxID", "rank", "taxName"});
  rep.printReport("kraken");
  return 0;
}

void parse_command_line(int argc, char **argv) {
  int opt;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "i:d:o:c:")) != -1) {
    switch (opt) {
    case 'd':
      TaxDB_filename = optarg;
      break;
    case 'c':
      CountsDB_filename = optarg;
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
       << "* -d filename      Kraken taxDB filename" << endl
       << "* -i filename      Input hll filename. May be repeated" << endl
       << "* -o filename      Output file for Kraken merged report" << endl
       << "  -c filename      Kraken *.kdb.counts filename" << endl;
  exit(exit_code);
}
