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

#ifdef VERIFY_MAPPER
#include "dbg.h"
#include "path_graph.h"
#endif

bool verifyGraph(const InputGraph& input_graph);
void verifyMapper(const InputGraph& graph);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: build_gcsa [options] base_name [base_name2 ..]" << std::endl;
    std::cerr << "  -b    Read the input in binary format (default)" << std::endl;
    std::cerr << "  -d N  Doubling steps (default and max " << GCSA::DOUBLING_STEPS << ")" << std::endl;
    std::cerr << "  -l N  Limit the size of the graph to N gigabytes (default " << GCSA::SIZE_LIMIT << ")" << std::endl;
    std::cerr << "  -o X  Use X as the base name for output (default: the first input)" << std::endl;
    std::cerr << "  -t    Read the input in text format" << std::endl;
    std::cerr << "  -T X  Use X as the directory for temporary files (default: " << DiskIO::DEFAULT_TEMP_DIR << ")" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  size_type doubling_steps = GCSA::DOUBLING_STEPS, size_limit = GCSA::SIZE_LIMIT;
  int c = 0;
  bool binary = true;
  std::string output_file;
  while((c = getopt(argc, argv, "bd:l:o:tT:")) != -1)
  {
    switch(c)
    {
    case 'b':
      binary = true; break;
    case 'd':
      doubling_steps =  std::stoul(optarg);
      if(doubling_steps > GCSA::DOUBLING_STEPS)
      {
        std::cerr << "build_gcsa: The number of doubling steps is too high: " << doubling_steps << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    case 'l':
      size_limit = std::stoul(optarg); break;
    case 'o':
      output_file = std::string(optarg) + GCSA::EXTENSION; break;
    case 't':
      binary = false; break;
    case 'T':
      DiskIO::setTemp(optarg); break;
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
    std::cout << "Input:           " << argv[i];
    if(binary) { std::cout << InputGraph::BINARY_EXTENSION << " (binary format)" << std::endl; }
    else { std::cout << InputGraph::TEXT_EXTENSION << " (text format)" << std::endl; }
  }
  std::cout << "Output:          " << output_file << std::endl;
  std::cout << "Doubling steps:  " << doubling_steps << std::endl;
  std::cout << "Size limit:      " << size_limit << " GB" << std::endl;
  std::cout << "Temp directory:  " << DiskIO::temp_dir << std::endl;
  std::cout << std::endl;

  InputGraph graph(argc - optind, argv + optind, binary);

#ifdef VERIFY_GRAPH
  if(!(verifyGraph(graph))) { std::exit(EXIT_FAILURE); }
#endif

#ifdef VERIFY_MAPPER
  verifyMapper(graph);
#endif

  GCSA index;
#ifdef LOAD_INDEX
  sdsl::load_from_file(index, output_file);
#else
  {
    double start = readTimer();
    GCSA temp(graph, doubling_steps, size_limit); index.swap(temp);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << "Memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
    std::cout << "I/O volume: " << inGigabytes(DiskIO::read_volume) << " GB read, "
              << inGigabytes(DiskIO::write_volume) << " GB write" << std::endl;
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
  index.verifyIndex(graph);
#endif

  std::cout << "Final memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

bool
verifyGraph(const InputGraph& input_graph)
{
  std::vector<KMer> kmers;
  input_graph.read(kmers);

  Alphabet alpha;
  std::map<std::string, std::pair<byte_type, byte_type>> graph;
  bool ok = true;
  for(size_type i = 0; i < kmers.size(); i++)
  {
    // We don't verify the edge from the sink to the source.
    std::string kmer = Key::decode(kmers[i].key, input_graph.k(), alpha);
    byte_type pred = (kmer[input_graph.k() - 1] == '#' ? 0 : Key::predecessors(kmers[i].key));
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
          std::string(1, (char)(alpha.comp2char[i])) + iter->first.substr(0, input_graph.k() - 1);
        if(graph.find(backward_pattern) == graph.end())
        {
          std::cerr << "build_gcsa: verifyGraph(): Node " << iter->first << " is missing predecessor("
                    << (char)(alpha.comp2char[i]) << "): "
                    << backward_pattern << std::endl;
          ok = false;
        }
        else
        {
          comp_type last = alpha.char2comp[iter->first[input_graph.k() - 1]];
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

#ifdef VERIFY_MAPPER
void
verifyMapper(const InputGraph& graph)
{
    std::vector<key_type> keys; graph.read(keys);
    DeBruijnGraph mapper(keys, graph.k());
#ifdef VERBOSE_STATUS_INFO
    std::cerr << "build_gcsa: verifyMapper(): Mapper has " << mapper.size() << " nodes, "
              << mapper.edge_count() << " edges" << std::endl;
    std::cerr << "build_gcsa: verifyMapper(): Mapper size: " << sdsl::size_in_bytes(mapper) << " bytes" << std::endl;
#endif

  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    std::string kmer = Key::decode(keys[i], graph.k(), mapper.alpha);
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


  ok = true;
  Alphabet alpha;
  LCP lcp(keys, graph.k());
  for(size_type i = 1; i < keys.size(); i++)
  {
    std::string left = Key::decode(keys[i - 1], graph.k(), alpha);
    std::string right = Key::decode(keys[i], graph.k(), alpha);
    size_type real_lcp = 0;
    while(real_lcp < graph.k() && left[real_lcp] == right[real_lcp]) { real_lcp++; }
    if(lcp.kmer_lcp[i] != real_lcp)
    {
      std::cerr << "build_gcsa: verifyMapper(): lcp(" << left << ", " << right << ") = " << real_lcp
                << ", got " << lcp.kmer_lcp[i] << std::endl;
      ok = false;
    }
  }
  std::cout << "LCP verification " << (ok ? "complete." : "failed.") << std::endl;

  std::cout << std::endl;
}
#endif

//------------------------------------------------------------------------------
