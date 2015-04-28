/*
  Copyright (c) 2015 Genome Research Ltd.

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <sstream>
#include <string>

#include "gcsa.h"

using namespace gcsa;

//------------------------------------------------------------------------------

struct KMer;

size_type readKMers(const std::string& base_name, std::vector<size_type>& keys, std::vector<KMer>& kmers);
void filterKMers(std::vector<size_type>& kmers);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc != 2)
  {
    std::cerr << "Usage: build_gcsa base_name" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "GCSA builder" << std::endl;
  std::cout << std::endl;
  std::string base_name = argv[1];
  std::cout << "Input: " << base_name << std::endl;
  std::cout << std::endl;

  std::vector<size_type> keys;
  std::vector<KMer> kmers;
  readKMers(base_name, keys, kmers);
  filterKMers(keys);
  GCSA index(keys);
  store_to_file(index, base_name + GCSA::EXTENSION);

  /*
    FIXME Sort the kmers and replace the keys with ranks based on the kmer part of the old key.
    Then verify that everything went fine:
      - unpack the kmer from 'keys'
      - search for it in the GCSA
      - the returned range should contain only the index of the kmer in 'keys'
  */

  std::cout << "Nodes: " << index.size() << ", edges: " << index.edge_count() << std::endl;
  std::cout << "GCSA size: " << size_in_bytes(index) << " bytes" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

struct KMer
{
  size_type key;
  size_type from, to;

  const static size_type OFFSET_BITS = 6;
  const static size_type OFFSET_MASK = 63;

  KMer() {}

  KMer(const std::vector<std::string>& tokens, const Alphabet& alpha, size_type successor)
  {
    uint8_t predecessors = chars(tokens[2], alpha);
    uint8_t successors = chars(tokens[3], alpha);
    this->key = GCSA::encode(alpha, tokens[0], predecessors, successors);

    this->from = encodePosition(tokens[1]);
    this->to = encodePosition(tokens[successor]);
  }

  inline bool
  operator< (const KMer& another) const
  {
    return (GCSA::kmer(this->key) < GCSA::kmer(another.key));
  }

  inline static size_type node(size_type value) { return value >> OFFSET_BITS; }
  inline static size_type offset(size_type value) { return value & OFFSET_MASK; }

  static bool
  tokenize(const std::string& line, std::vector<std::string>& tokens)
  {
    {
      std::string token;
      std::istringstream ss(line);
      while(std::getline(ss, token, '\t'))
      {
        tokens.push_back(token);
      }
      if(tokens.size() < 4 || tokens.size() > 5)
      {
        std::cerr << "KMer::tokenize(): The kmer line must contain 4 or 5 tokens." << std::endl;
        std::cerr << "KMer::tokenize(): The line was: " << line << std::endl;
        return false;
      }
    }

    // Split the list of successor positions into separate tokens.
    if(tokens.size() == 5)
    {
      std::string destinations = tokens[4], token;
      std::istringstream ss(destinations);
      tokens.resize(4);
      while(std::getline(ss, token, ','))
      {
        tokens.push_back(token);
      }
    }
    else  // Use the source node as the destination.
    {
      tokens.push_back(tokens[1]);
    }

    return true;
  }

  static size_type
  encodePosition(const std::string& token)
  {
    size_t separator = 0;
    size_type _node = std::stoul(token, &separator);
    if(separator >= token.length())
    {
      std::cerr << "KMer::encodePosition(): Invalid position token " << token << std::endl;
      return 0;
    }

    std::string temp = token.substr(separator + 1);
    size_type _offset = std::stoul(temp);
    if(_offset > OFFSET_MASK)
    {
      std::cerr << "KMer::encodePosition(): Offset " << _offset << " too large!" << std::endl;
      return 0;
    }

    return (_node << 6) | _offset;
  }

  static uint8_t
  chars(const std::string& token, const Alphabet& alpha)
  {
    uint8_t val = 0;
    for(size_type i = 0; i < token.length(); i += 2) { val |= 1 << alpha.char2comp[token[i]]; }
    return val;
  }
};

std::ostream&
operator<< (std::ostream& out, const KMer& kmer)
{
  out << "(key " << GCSA::kmer(kmer.key)
      << ", in " << (size_type)(GCSA::predecessors(kmer.key))
      << ", out " << (size_type)(GCSA::successors(kmer.key))
      << ", from " << range_type(KMer::node(kmer.from), KMer::offset(kmer.from))
      << ", to " << range_type(KMer::node(kmer.to), KMer::offset(kmer.to)) << ")";
  return out;
}

//------------------------------------------------------------------------------

size_type
readKMers(const std::string& base_name, std::vector<size_type>& keys, std::vector<KMer>& kmers)
{
  std::string filename = base_name + ".gcsa2";
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "build_gcsa: readKMers(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  Alphabet alpha;
  size_type kmer_length = 0, sink_node = 0;
  while(input)
  {
    std::string line;
    std::getline(input, line);
    if(line.length() == 0) { continue; }

    std::vector<std::string> tokens;
    if(!(KMer::tokenize(line, tokens))) { continue; }
    if(kmer_length > 0 && tokens[0].length() != kmer_length)
    {
      std::cerr << "build_gcsa: readKMers(): kmer length changed from " << kmer_length
                << " to " << tokens[0].length() << std::endl;
    }
    kmer_length = tokens[0].length();

    for(size_type successor = 4; successor < tokens.size(); successor++)
    {
      KMer kmer(tokens, alpha, successor); kmers.push_back(kmer);
      if(successor == 4)
      {
        if(kmer.from == kmer.to) { sink_node = KMer::node(kmer.from); }
        keys.push_back(kmer.key);
      }
    }
  }
  input.close();

  // If the kmer includes one or more endmarkers, the successor position is past
  // the GCSA sink node. We set to = from for such nodes to denote that the kmer
  // does not have successors.
  for(size_type i = 0; i < kmers.size(); i++)
  {
    if(KMer::node(kmers[i].to) == sink_node && KMer::offset(kmers[i].to) > 0)
    {
      kmers[i].to = kmers[i].from;
    }
  }

  std::cout << "Read " << kmers.size() << " kmers of length " << kmer_length << std::endl;
  std::cout << std::endl;
  return kmer_length;
}

//------------------------------------------------------------------------------

void
filterKMers(std::vector<size_type>& kmers)
{
  if(kmers.empty()) { return; }
  parallelQuickSort(kmers.begin(), kmers.end());

  size_type size = 0;
  for(size_type i = 1; i < kmers.size(); i++)
  {
    if(GCSA::kmer(kmers[size]) == GCSA::kmer(kmers[i]))
    {
      kmers[size] = GCSA::merge(kmers[size], kmers[i]);
    }
    else { size++; }
  }
  size++;
  kmers.resize(size);

  std::cout << "Unique kmers: " << kmers.size() << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
