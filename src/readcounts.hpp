/*
 * Copyright 2017-2018, Florian Breitwieser
 *
 * This file is part of the KrakenUniq taxonomic sequence classification system.
 *
 * KrakenUniq is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KrakenUniq is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef READCOUNTS_HPP
#define READCOUNTS_HPP

#include "kraken_headers.hpp"
#include "hyperloglogplus.hpp"
#include "khset.h"
#include <unordered_set>

namespace kraken {
  static size_t HLL_PRECISION = 14;

  template <typename CONTAINER>
  class ReadCounts {

  public:
    uint64_t readCount() const { return n_reads; }
    void incrementReadCount() { ++n_reads; }
    uint64_t kmerCount() const { return n_kmers; }
    uint64_t uniqueKmerCount() const; // to be implemented for each CONTAINER
    void serialize(ostringstream& os) const;
    string serialize() const;
    void deserialize(string& serialized);

    ReadCounts() : n_reads(0), n_kmers(0) {
    }


    ReadCounts(uint64_t _n_reads, uint64_t _n_kmers, const CONTAINER& _kmers) :
            n_reads(_n_reads), n_kmers(_n_kmers), kmers(_kmers) {
      kmers.set_nObserved(n_kmers);
    }

    //ReadCounts(const ReadCounts& other) = delete;
    
    ReadCounts(const ReadCounts& other) : n_reads(other.n_reads), n_kmers(other.n_kmers), kmers(other.kmers) {
    }

    ReadCounts(ReadCounts&& other) : n_reads(other.n_reads), n_kmers(other.n_kmers), kmers(std::move(other.kmers)) {
    }

    ReadCounts(string& serialized) {
      deserialize(serialized);
    }

    ReadCounts& operator=(const ReadCounts& other) {
      n_reads = other.n_reads;
      n_kmers =other.n_kmers;
      kmers = other.kmers;
      return *this;
    }


    ReadCounts& operator=(ReadCounts&& other) {
      n_reads = other.n_reads;
      n_kmers =other.n_kmers;
      kmers = std::move(other.kmers);
      return *this;
    }
  
    void add_kmer(uint64_t kmer) {
      ++n_kmers;
      kmers.insert(kmer);
    }

    ReadCounts& operator+=(const ReadCounts& other) {
      n_reads += other.n_reads;
      n_kmers += other.n_kmers;
      kmers += other.kmers;
      return *this;
    }

    ReadCounts& operator+=(ReadCounts&& other) {
      n_reads += other.n_reads;
      n_kmers += other.n_kmers;
      kmers += std::move(other.kmers);
      return *this;
    }

    bool operator<(const ReadCounts& other) {
      if (n_reads < other.n_reads) {
        return true;
      }
      if (n_reads == other.n_reads && n_kmers < other.n_kmers) {
        return true;
      }
      return false;
    }


  private:
    uint64_t n_reads;
    uint64_t n_kmers; 
    CONTAINER kmers; // unique k-mer count per taxon

  };

  // Overload operator += for set, so that it can be used for merging
  template <typename T>
  unordered_set<T>& operator+=(unordered_set<T>& left, const unordered_set<T>& right) {
    left.insert(right.begin(), right.end());
    return left;
  }

  template <typename T>
  set<T>& operator+=(set<T>& left, const set<T>& right) {
    left.insert(right.begin(), right.end());
    return left;
  }
  
  template<>
  uint64_t ReadCounts< HyperLogLogPlusMinus<uint64_t> >::uniqueKmerCount() const {
    return(kmers.cardinality());
  }

  template<typename T>
  uint64_t ReadCounts< T >::uniqueKmerCount() const {
    return(kmers.size());
  }

  template<typename T>
  void ReadCounts< T >::serialize(ostringstream& os) const {
    os << n_reads << "\t";
    os << n_kmers << "\t";
    os << kmers.serialize();
  }

  template<typename T>
  string ReadCounts< T >::serialize() const {
    string out;
    out += to_string(n_reads);
    out += "\t";
    out += to_string(n_kmers);
    out += "\t";
    out += kmers.serialize();
    return out;
  }

  template<typename T>
  void ReadCounts< T >::deserialize(string& serialized) {
    stringstream iss(serialized);
    iss >> n_reads;
    iss >> n_kmers;
    iss >> ws;
    string iss_left;
    getline(iss, iss_left);
    kmers.deserialize(iss_left);
    kmers.set_nObserved(n_kmers);
  }

}
#endif
