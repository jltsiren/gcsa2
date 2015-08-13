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
#include <unistd.h>

#include "files.h"
#include "gcsa.h"

using namespace gcsa;

//------------------------------------------------------------------------------

/*
  Various verification options for debugging purposes.
*/

#ifdef VERIFY_CONSTRUCTION
//#define VERIFY_GRAPH
//#define VERIFY_MAPPER
//#define LOAD_INDEX
#define VERIFY_INDEX
#endif

//------------------------------------------------------------------------------

bool verifyGraph(const std::vector<KMer>& kmers, size_type kmer_length);
bool verifyMapper(const DeBruijnGraph& mapper, const std::vector<key_type>& keys, size_type kmer_length);
void verifyIndex(const GCSA& index, std::vector<KMer>& kmers, size_type kmer_length);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: build_gcsa [options] base_name [base_name2 ..]" << std::endl;
    std::cerr << "  -b    Read the input in binary format" << std::endl;
    std::cerr << "  -d N  Doubling steps (default and max " << GCSA::DOUBLING_STEPS << ")" << std::endl;
    std::cerr << "  -l N  Limit the size of the graph to N gigabytes (default 200)" << std::endl;
    std::cerr << "  -o X  Use X as the base name for output (default: the first input)" << std::endl;
    std::cerr << "  -t    Read the input in text format (default)" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  size_type doubling_steps = GCSA::DOUBLING_STEPS, size_limit = GCSA::SIZE_LIMIT;
  int c = 0; bool binary = false;
  std::string output_file;
  while((c = getopt(argc, argv, "bd:l:o:t")) != -1)
  {
    switch(c)
    {
    case 'b':
      binary = true; break;
    case 'd':
      doubling_steps =  std::stoul(optarg);
      if(doubling_steps > GCSA::DOUBLING_STEPS)
      {
        std::cerr << "build_gcsa: Number of doubling steps is too high: " << doubling_steps << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'l':
      size_limit = std::stoul(optarg);
      break;
    case 'o':
      output_file = std::string(optarg) + GCSA::EXTENSION;
      break;
    case 't':
      binary = false; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(optind >= argc)
  {
    std::cerr << "build_gcsa: No input files specified" << std::endl;
    std::exit(EXIT_SUCCESS);
  }
  if(output_file.length() == 0)
  {
    output_file = std::string(argv[optind]) + GCSA::EXTENSION;
  }

  std::cout << "GCSA builder" << std::endl;
  std::cout << std::endl;
  for(int i = optind; i < argc; i++)
  {
    std::cout << "Input: " << argv[i];
    if(binary) { std::cout << BINARY_EXTENSION << " (binary format)" << std::endl; }
    else { std::cout << TEXT_EXTENSION << " (text format)" << std::endl; }
  }
  std::cout << "Output: " << output_file << std::endl;
  std::cout << "Doubling steps: " << doubling_steps << std::endl;
  std::cout << "Size limit: " << size_limit << " GB" << std::endl;
  std::cout << std::endl;

#ifdef VERIFY_GRAPH
  {
    std::vector<KMer> kmers;
    size_type kmer_length = readKMers(argc - optind, argv + optind, kmers, binary, true);
    if(!(verifyGraph(kmers, kmer_length))) { std::exit(EXIT_FAILURE); }
  }
#endif


#ifdef VERIFY_MAPPER
  {
    std::vector<KMer> kmers;
    std::vector<key_type> keys;
    sdsl::int_vector<0> last_chars;
    size_type kmer_length = readKMers(argc - optind, argv + optind, kmers, binary, true);
    uniqueKeys(kmers, keys, last_chars, true);
    DeBruijnGraph mapper(keys, kmer_length);
    std::cout << "Nodes: " << mapper.size() << ", edges: " << mapper.edge_count() << std::endl;
    std::cout << "Mapper size: " << sdsl::size_in_bytes(mapper) << " bytes" << std::endl;
    verifyMapper(mapper, keys, kmer_length);
  }
#endif

  GCSA index;
#ifdef LOAD_INDEX
  sdsl::load_from_file(index, output_file);
#else
  {
    std::vector<KMer> kmers;
    size_type kmer_length = readKMers(argc - optind, argv + optind, kmers, binary, true);
    double start = readTimer();
    GCSA temp(kmers, kmer_length, doubling_steps, size_limit); index.swap(temp);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << std::endl;
    sdsl::store_to_file(index, output_file);
  }
#endif

  printHeader("Paths"); std::cout << index.size() << std::endl;
  printHeader("Edges"); std::cout << index.edge_count() << std::endl;
  printHeader("Samples");
  std::cout << index.sample_count() << " (at " << index.sampled_positions() << " positions, "
            << index.sample_bits() << " bits each)" << std::endl;
  printHeader("Max query"); std::cout << index.order() << std::endl;

  size_type index_bytes = sdsl::size_in_bytes(index);
  size_type sample_bytes = sdsl::size_in_bytes(index.stored_samples);
  printHeader("Index size"); std::cout << inMegabytes(index_bytes) << " MB" << std::endl;
  printHeader("Without samples"); std::cout << inMegabytes(index_bytes - sample_bytes) << " MB" << std::endl;
  std::cout << std::endl;

#ifdef VERIFY_INDEX
  {
    std::vector<KMer> kmers;
    size_type kmer_length = readKMers(argc - optind, argv + optind, kmers, binary, true);
    verifyIndex(index, kmers, kmer_length);
  }
#endif

  std::cout << "Memory usage: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
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
    std::string kmer = Key::decode(kmers[i].key, kmer_length, alpha);
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
          std::cerr << "build_gcsa: verifyGraph(): Node " << iter->first << " is missing predecessor("
                    << (char)(alpha.comp2char[i]) << "): "
                    << backward_pattern << std::endl;
          ok = false;
        }
        else
        {
          comp_type last = alpha.char2comp[iter->first[kmer_length - 1]];
          if((graph[backward_pattern].second & (1 << last)) == 0)
          {
            std::cerr << "build_gcsa: verifyGraph(): Reverse: Node " << backward_pattern << " is missing successor("
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
          std::cerr << "build_gcsa: verifyGraph(): Node " << iter->first << " is missing successor("
                    << (char)(alpha.comp2char[i]) << "): "
                    << forward_pattern << std::endl;
          ok = false;
        }
        else
        {
          comp_type first = alpha.char2comp[iter->first[0]];
          if((graph[forward_pattern].first & (1 << first)) == 0)
          {
            std::cerr << "build_gcsa: verifyGraph(): Reverse: Node " << forward_pattern << " is missing predecessor("
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
verifyMapper(const DeBruijnGraph& mapper, const std::vector<key_type>& keys, size_type kmer_length)
{
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    std::string kmer = Key::decode(keys[i], kmer_length, mapper.alpha);
    size_type endmarker_pos = kmer.find('$'); // The actual kmer ends at the first endmarker.
    if(endmarker_pos != std::string::npos) { kmer = kmer.substr(0, endmarker_pos + 1); }
    range_type range = mapper.find(kmer);
    if(range != range_type(i, i))
    {
      std::cerr << "build_gcsa: verifyMapper(): find(" << kmer << ") = " << range
                << ", expected " << range_type(i, i) << std::endl;
      ok = false;
    }
  }

  std::cout << "Mapper verification " << (ok ? "complete." : "failed.") << std::endl;
  std::cout << std::endl;
  return ok;
}

//------------------------------------------------------------------------------

std::ostream&
printOccs(const std::vector<node_type>& occs, std::ostream& out)
{
  out << "{";
  for(size_type i = 0; i < occs.size(); i++)
  {
    out << (i == 0 ? " " : ", ") << Node::decode(occs[i]);
  }
  out << " }";
  return out;
}

void
verifyIndex(const GCSA& index, std::vector<KMer>& kmers, size_type kmer_length)
{
  parallelQuickSort(kmers.begin(), kmers.end());

  size_type i = 0, fails = 0;
  while(i < kmers.size())
  {
    size_type next = i + 1;
    while(next < kmers.size() && Key::label(kmers[next].key) == Key::label(kmers[i].key)) { next++; }

    std::string kmer = Key::decode(kmers[i].key, kmer_length, index.alpha);
    size_type endmarker_pos = kmer.find('$'); // The actual kmer ends at the first endmarker.
    if(endmarker_pos != std::string::npos) { kmer = kmer.substr(0, endmarker_pos + 1); }

    range_type range = index.find(kmer);
    if(Range::empty(range))
    {
      std::cerr << "build_gcsa: verifyIndex(): find(" << kmer << ") returned empty range" << std::endl;
      i = next; fails++; continue;
    }

    std::vector<node_type> expected;
    for(size_type j = i; j < next; j++) { expected.push_back(kmers[j].from); }
    removeDuplicates(expected, false);
    std::vector<node_type> occs;
    index.locate(range, occs);

    bool failed = false;
    if(occs.size() != expected.size())
    {
      std::cerr << "build_gcsa: verifyIndex(): Expected " << expected.size()
                << " occurrences, got " << occs.size() << std::endl;
      failed = true;
    }
    else
    {
      for(size_type j = 0; j < occs.size(); j++)
      {
        if(occs[j] != expected[j])
        {
          std::cerr << "build_gcsa: verifyIndex(): Failure at " << j << ": "
                    << "expected " << Node::decode(expected[j])
                    << ", got " << Node::decode(occs[j]) << std::endl;
          failed = true; break;
        }
      }
    }
    if(failed)
    {
      std::cerr << "build_gcsa: verifyIndex(): locate(" << kmer << ") failed" << std::endl;
      std::cerr << "build_gcsa: verifyIndex(): Expected ";
      printOccs(expected, std::cerr) << std::endl;
      std::cerr << "build_gcsa: verifyIndex(): Got ";
      printOccs(occs, std::cerr) << std::endl;
      fails++;
    }

    i = next;
  }

  if(fails == 0)
  {
    std::cout << "Index verification complete." << std::endl;
  }
  else
  {
    std::cout << "Index verification failed for " << fails << " patterns." << std::endl;
  }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
