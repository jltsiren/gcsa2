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
#include <string>

#include "gcsa.h"

using namespace gcsa;

//------------------------------------------------------------------------------

#define VERIFY_GRAPH
#define VERIFY_MAPPER
#define VERIFY_INDEX

size_type readKMers(const std::string& base_name, std::vector<KMer>& kmers, bool print = false);

bool verifyGraph(const std::vector<KMer>& kmers, size_type kmer_length);
bool verifyMapper(const GCSA& mapper, const std::vector<key_type>& keys, size_type kmer_length);

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
  {
    std::vector<KMer> kmers;
    size_type kmer_length = readKMers(base_name, kmers, true);
    if(!(verifyGraph(kmers, kmer_length))) { return 2; }
  }
#endif


#ifdef VERIFY_MAPPER
  {
    std::vector<KMer> kmers;
    std::vector<key_type> keys;
    sdsl::int_vector<0> last_chars;
    size_type kmer_length = readKMers(base_name, kmers, true);
    uniqueKeys(kmers, keys, last_chars, true);
    GCSA mapper(keys, kmer_length);
    std::cout << "Nodes: " << mapper.size() << ", edges: " << mapper.edge_count() << std::endl;
    std::cout << "GCSA size: " << sdsl::size_in_bytes(mapper) << " bytes" << std::endl;
    verifyMapper(mapper, keys, kmer_length);
  }
#endif

  std::vector<KMer> kmers;
  size_type kmer_length = readKMers(base_name, kmers);
  GCSA index(kmers, kmer_length);
  sdsl::store_to_file(index, base_name + GCSA::EXTENSION);

#ifdef VERIFY_INDEX
#endif

  std::cout << "Memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

size_type
readKMers(const std::string& base_name, std::vector<KMer>& kmers, bool print)
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
  sdsl::util::clear(kmers);
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
        if(kmer.sorted()) { sink_node = Node::id(kmer.from); }
      }
    }
  }
  input.close();

  // If the kmer includes one or more endmarkers, the successor position is past
  // the GCSA sink node. Those kmers are marked as sorted, as they cannot be
  // extended.
  for(size_type i = 0; i < kmers.size(); i++)
  {
    if(Node::id(kmers[i].to) == sink_node && Node::offset(kmers[i].to) > 0)
    {
      kmers[i].makeSorted();
    }
  }

  if(print)
  {
    std::cout << "Read " << kmers.size() << " kmers of length " << kmer_length << std::endl;
  }
  return kmer_length;
}

//------------------------------------------------------------------------------

bool
verifyGraph(const std::vector<KMer>& kmers, size_type kmer_length)
{
  Alphabet alpha;
  std::map<std::string, std::pair<byte_type, byte_type>> graph;
  bool ok = true;
  for(size_type i = 0; i < kmers.size(); i++)
  {
    // We don't verify the edge from the sink to the source.
    std::string kmer = Key::decode(alpha, kmers[i].key, kmer_length);
    byte_type pred = (kmer[kmer_length - 1] == '#' ? 0 : Key::predecessors(kmers[i].key));
    byte_type succ = (kmer[0] == '$' ? 0 : Key::successors(kmers[i].key));
    if(graph.find(kmer) == graph.end())
    {
      graph[kmer] = std::make_pair(pred, succ);
    }
    else
    {
      graph[kmer].first |= pred;
      graph[kmer].second |= succ;
    }
  }

  for(auto iter = graph.begin(); iter != graph.end(); ++iter)
  {
    for(size_type i = 1; i < alpha.sigma; i++)
    {
      if(iter->second.first & (1 << i))
      {
        std::string backward_pattern =
          std::string(1, (char)(alpha.comp2char[i])) + iter->first.substr(0, kmer_length - 1);
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

  std::cout << "Graph verification " << (ok ? "complete." : "failed.") << std::endl;
  std::cout << std::endl;
  return ok;
}

//------------------------------------------------------------------------------

bool
verifyMapper(const GCSA& mapper, const std::vector<key_type>& keys, size_type kmer_length)
{
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    std::string kmer = Key::decode(mapper.alpha, keys[i], kmer_length);
    size_type endmarker_pos = kmer.find('$'); // The actual kmer ends at the first endmarker.
    if(endmarker_pos != std::string::npos) { kmer = kmer.substr(0, endmarker_pos + 1); }
    range_type range = mapper.find(kmer);
    if(range != range_type(i, i))
    {
      std::cerr << "build_gcsa: find(" << kmer << ") = " << range
                << ", expected " << range_type(i, i) << std::endl;
      ok = false;
    }
  }
  std::cout << "Mapper verification " << (ok ? "complete." : "failed.") << std::endl;
  std::cout << std::endl;
  return ok;
}

//------------------------------------------------------------------------------
