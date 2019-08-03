#include "hyperloglogplus.hpp"
#include "kraken_headers.hpp"
#include "readcounts.hpp"
// #include "krakendb.hpp"
// #include "quickfile.hpp"
// #include "krakendb.hpp"
#include <iostream>
#include <fstream>
#include <random>

using namespace std;
using namespace kraken;

int main(int argc, char **argv) {

  HyperLogLogPlusMinus<uint64_t> hll; // unique k-mer count per taxon
  hll.use_n_observed = false;

  for (uint64_t i = 0; i < 1000; ++i) {
    hll.insert(i);
  }

  string ser = hll.serialize();

  HyperLogLogPlusMinus<uint64_t> hll2(ser);

  ReadCounts<HyperLogLogPlusMinus<uint64_t>> rc;
  for (uint64_t i = 0; i < 1000; ++i) {
    rc.add_kmer(i);
  }

  ostringstream os;
  rc.serialize(os);
  cout << os.str() << endl;
}
