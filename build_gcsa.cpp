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

#include <map>
#include <sstream>
#include <string>

#include "gcsa.h"

using namespace gcsa;

//------------------------------------------------------------------------------

#define VERIFY_GRAPH
#define VERIFY_INDEX

struct KMer;

size_type readKMers(const std::string& base_name, std::vector<KMer>& kmers, std::vector<size_type>& keys);
void filterKeys(std::vector<size_type>& keys);
void sortKMers(std::vector<KMer>& kmers, const std::vector<size_type>& keys);

bool verifyGraph(const std::string& base_name);
bool verifyIndex(const GCSA& index, const std::vector<size_type>& keys, size_type kmer_length);

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

#ifdef VERIFY_GRAPH
  if(!(verifyGraph(base_name))) { return 2; }
#endif

  std::vector<size_type> keys;
  std::vector<KMer> kmers;
  size_type kmer_length = readKMers(base_name, kmers, keys);
  filterKeys(keys);
  sortKMers(kmers, keys);
  GCSA index(keys, kmer_length);
  sdsl::store_to_file(index, base_name + GCSA::EXTENSION);
  std::cout << "Nodes: " << index.size() << ", edges: " << index.edge_count() << std::endl;
  std::cout << "GCSA size: " << sdsl::size_in_bytes(index) << " bytes" << std::endl;
  std::cout << std::endl;

#ifdef VERIFY_INDEX
  verifyIndex(index, keys, kmer_length);
#endif

  std::cout << "Memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
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
    byte_type predecessors = chars(tokens[2], alpha);
    byte_type successors = chars(tokens[3], alpha);
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

  static byte_type
  chars(const std::string& token, const Alphabet& alpha)
  {
    byte_type val = 0;
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

inline bool
operator< (size_type key, const KMer& kmer)
{
  return (GCSA::kmer(key) < GCSA::kmer(kmer.key));
}

//------------------------------------------------------------------------------

size_type
readKMers(const std::string& base_name, std::vector<KMer>& kmers, std::vector<size_type>& keys)
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
  return kmer_length;
}

//------------------------------------------------------------------------------

void
filterKeys(std::vector<size_type>& keys)
{
  if(keys.empty()) { return; }
  parallelQuickSort(keys.begin(), keys.end());

  size_type tail = 0;
  for(size_type i = 1; i < keys.size(); i++)
  {
    if(GCSA::kmer(keys[i]) == GCSA::kmer(keys[tail]))
    {
      keys[tail] = GCSA::merge(keys[tail], keys[i]);
    }
    else
    {
      tail++; keys[tail] = keys[i];
    }
  }
  keys.resize(tail + 1);

  std::cout << "Unique keys: " << keys.size() << std::endl;
}

//------------------------------------------------------------------------------

void
sortKMers(std::vector<KMer>& kmers, const std::vector<size_type>& keys)
{
  if(kmers.empty()) { return; }
  parallelQuickSort(kmers.begin(), kmers.end());

  // Replace the keys with ranks in the sorted key array.
  for(size_type kmer = 0, key = 0; kmer < kmers.size(); kmer++)
  {
    while(keys[key] < kmers[kmer]) { key++; }
    kmers[kmer].key = key;
  }
}

//------------------------------------------------------------------------------

bool
verifyGraph(const std::string& base_name)
{
  std::string filename = base_name + ".gcsa2";
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "build_gcsa: verifyGraph(): Cannot open input file " << filename << std::endl;
    return false;
  }

  Alphabet alpha;
  std::map<std::string, std::pair<byte_type, byte_type>> graph;
  size_type kmer_length = 0;
  bool ok = true;
  while(input)
  {
    std::string line;
    std::getline(input, line);
    if(line.length() == 0) { continue; }

    std::vector<std::string> tokens;
    if(!(KMer::tokenize(line, tokens))) { continue; }
    if(kmer_length > 0 && tokens[0].length() != kmer_length)
    {
      std::cerr << "build_gcsa: verifyGraph(): kmer length changed from " << kmer_length
                << " to " << tokens[0].length() << std::endl;
    }
    kmer_length = tokens[0].length();

    // We don't verify the edge from the sink to the source.
    byte_type pred = (tokens[0][kmer_length - 1] == '#' ? 0 : KMer::chars(tokens[2], alpha));
    byte_type succ = (tokens[0][0] == '$' ? 0 : KMer::chars(tokens[3], alpha));
    if(graph.find(tokens[0]) == graph.end())
    {
      graph[tokens[0]] = std::make_pair(pred, succ);
    }
    else
    {
      graph[tokens[0]].first |= pred;
      graph[tokens[0]].second |= succ;
    }
  }
  input.close();

  for(auto iter = graph.begin(); iter != graph.end(); ++iter)
  {
    for(size_type i = 1; i < alpha.sigma; i++)
    {
      if(iter->second.first & (1 << i))
      {
        std::string backward_pattern = std::string(1, (char)(alpha.comp2char[i])) + iter->first.substr(0, kmer_length - 1);
        if(graph.find(backward_pattern) == graph.end())
        {
          std::cerr << "Node " << iter->first << " is missing predecessor("
                    << (char)(alpha.comp2char[i]) << "): "
                    << backward_pattern << std::endl;
          ok = false;
        }
        else
        {
          char_type last = alpha.char2comp[iter->first[kmer_length - 1]];
          if((graph[backward_pattern].second & (1 << last)) == 0)
          {
            std::cerr << "Reverse: Node " << backward_pattern << " is missing successor("
                      << (char)(alpha.comp2char[last]) << "): "
                      << iter->first << std::endl;
            ok = false;
          }
        }
      }
      if(iter->second.second & (1 << i))
      {
        std::string forward_pattern = iter->first.substr(1) + (char)(alpha.comp2char[i]);
        if(graph.find(forward_pattern) == graph.end())
        {
          std::cerr << "Node " << iter->first << " is missing successor("
                    << (char)(alpha.comp2char[i]) << "): "
                    << forward_pattern << std::endl;
          ok = false;
        }
        else
        {
          char_type first = alpha.char2comp[iter->first[0]];
          if((graph[forward_pattern].first & (1 << first)) == 0)
          {
            std::cerr << "Reverse: Node " << forward_pattern << " is missing predecessor("
                      << (char)(alpha.comp2char[first]) << "): "
                      << iter->first << std::endl;
            ok = false;
          }
        }
      }
    }
  }

  std::cout << "Verification " << (ok ? "completed." : "failed.") << std::endl;
  std::cout << std::endl;
  return ok;
}

//------------------------------------------------------------------------------

bool
verifyIndex(const GCSA& index, const std::vector<size_type>& keys, size_type kmer_length)
{
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    std::string kmer = GCSA::decode(index.alpha, keys[i], kmer_length);
    range_type range = index.find(kmer);
    if(range != range_type(i, i))
    {
      std::cerr << "build_gcsa: find(" << kmer << ") = " << range
                << ", expected " << range_type(i, i) << std::endl;
      ok = false;
    }
    std::cout << kmer << std::endl;
  }
  std::cout << "Verification " << (ok ? "completed." : "failed.") << std::endl;
  std::cout << std::endl;
  return ok;
}

//------------------------------------------------------------------------------
