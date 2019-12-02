/*
  Copyright (c) 2018, 2019 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#include <string>
#include <unistd.h>

#include <gcsa/algorithms.h>

using namespace gcsa;

//------------------------------------------------------------------------------

const size_type INDENT = 20;

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    Version::print(std::cerr, "GCSA2 builder");
    std::cerr << "Usage: build_gcsa [options] base_name [base_name2 ..]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Input/output options:" << std::endl;
    std::cerr << "  -b    Read the input in binary format (default)" << std::endl;
    std::cerr << "  -t    Read the input in text format" << std::endl;
    std::cerr << "  -o X  Use X as the base name for output (default: the first input)" << std::endl;
    std::cerr << "Index construction options:" << std::endl;
    std::cerr << "  -d N  Doubling steps (default " << ConstructionParameters::DOUBLING_STEPS << ", max " << ConstructionParameters::MAX_STEPS << ")" << std::endl;
    std::cerr << "  -m X  Use node mapping from file X" << std::endl;
    std::cerr << "  -s N  Use sample period N (default " << ConstructionParameters::SAMPLE_PERIOD << ")" << std::endl;
    std::cerr << "  -B N  Set LCP branching factor to N (default " << ConstructionParameters::LCP_BRANCHING << ")" << std::endl;
    std::cerr << "  -L    Load the index instead of building it" << std::endl;
    std::cerr << "  -v    Verify the index by querying it with the kmers" << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "  -D X  Use X as the directory for temporary files (default: " << TempFile::DEFAULT_TEMP_DIR << ")" << std::endl;
    std::cerr << "  -l N  Limit disk space usage to N gigabytes (default " << ConstructionParameters::SIZE_LIMIT << ")" << std::endl;
    std::cerr << "  -T N  Set the number of threads to N (default and max " << omp_get_max_threads() << " on this system)" << std::endl;
    std::cerr << "  -V N  Set verbosity level to N (default " << Verbosity::DEFAULT << ")" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  int c = 0;
  bool binary = true, load_index = false, verify = false;
  std::string index_file, lcp_file, mapping_file;
  ConstructionParameters parameters;
  while((c = getopt(argc, argv, "bto:d:m:s:B:LvD:l:T:V:")) != -1)
  {
    switch(c)
    {
    case 'b':
      binary = true; break;
    case 't':
      binary = false; break;
    case 'o':
      index_file = std::string(optarg) + GCSA::EXTENSION;
      lcp_file = std::string(optarg) + LCPArray::EXTENSION;
      break;
    case 'd':
      parameters.setSteps(std::stoul(optarg)); break;
    case 'm':
      mapping_file = optarg; break;
    case 's':
      parameters.setSamplePeriod(std::stoul(optarg)); break;
    case 'B':
      parameters.setLCPBranching(std::stoul(optarg)); break;
    case 'L':
      load_index = true; break;
    case 'v':
      verify = true; break;
    case 'D':
      TempFile::setDirectory(optarg); break;
    case 'l':
      parameters.setLimit(std::stoul(optarg)); break;
    case 'T':
      omp_set_num_threads(Range::bound(std::stoul(optarg), 1, omp_get_max_threads())); break;
    case 'V':
      Verbosity::set(std::stoul(optarg)); break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(optind >= argc)
  {
    std::cerr << "build_gcsa: No input files specified" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(index_file.empty())
  {
    index_file = std::string(argv[optind]) + GCSA::EXTENSION;
    lcp_file = std::string(argv[optind]) + LCPArray::EXTENSION;
  }

  Version::print(std::cout, "GCSA2 builder");
  for(int i = optind; i < argc; i++)
  {
    printHeader("Input", INDENT);
    std::cout << argv[i];
    if(binary) { std::cout << InputGraph::BINARY_EXTENSION << " (binary format)" << std::endl; }
    else { std::cout << InputGraph::TEXT_EXTENSION << " (text format)" << std::endl; }
  }
  if(!(mapping_file.empty()))
  {
    printHeader("Node mapping", INDENT); std::cout << mapping_file << std::endl;
  }
  printHeader("Output", INDENT); std::cout << index_file << ", " << lcp_file << std::endl;
  if(!load_index)
  {
    printHeader("Doubling steps", INDENT); std::cout << parameters.doubling_steps << std::endl;
    printHeader("Sample period", INDENT); std::cout << parameters.sample_period << std::endl;
    printHeader("Branching factor", INDENT); std::cout << parameters.lcp_branching << std::endl;
    printHeader("Temp directory", INDENT); std::cout << TempFile::temp_dir << std::endl;
    printHeader("Size limit", INDENT); std::cout << inGigabytes(parameters.size_limit) << " GB" << std::endl;
    printHeader("Threads", INDENT); std::cout << omp_get_max_threads() << std::endl;
    printHeader("Verbosity", INDENT); std::cout << Verbosity::levelName() << std::endl;
  }
  std::cout << std::endl;

  InputGraph graph(argc - optind, argv + optind, binary, Alphabet(), mapping_file);

  GCSA index;
  LCPArray lcp;
  if(load_index)
  {
    if(!sdsl::load_from_file(index, index_file))
    {
      std::cerr << "build_gcsa: Cannot load the index from " << index_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(!sdsl::load_from_file(lcp, lcp_file))
    {
      std::cerr << "build_gcsa: Cannot load the LCP array from " << lcp_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else
  {
    double start = readTimer();
    index = GCSA(graph, parameters);
    lcp = LCPArray(graph, parameters);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << "Memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
    std::cout << "I/O volume: " << inGigabytes(readVolume()) << " GB read, "
              << inGigabytes(writeVolume()) << " GB write" << std::endl;
    std::cout << std::endl;
    if(!sdsl::store_to_file(index, index_file))
    {
      std::cerr << "build_gcsa: Cannot write the index to " << index_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(!sdsl::store_to_file(lcp, lcp_file))
    {
      std::cerr << "build_gcsa: Cannot write the LCP array to " << lcp_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  printStatistics(index, lcp);

  if(verify) { verifyIndex(index, &lcp, graph); }

  std::cout << "Final memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
