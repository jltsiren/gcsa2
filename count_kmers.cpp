/*
  Copyright (c) 2016 Genome Research Ltd.

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

#include "algorithms.h"

using namespace gcsa;

//------------------------------------------------------------------------------

const size_type DEFAULT_K = 32;

void countKmers(const std::string& base_name, size_type k, bool force, bool include_Ns);
void compareKmers(const std::string& left_name, const std::string& right_name, size_type k, bool force, bool include_Ns);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "usage: count_kmers [options] base_name [base_name2]" << std::endl;
    std::cerr << "  -f    Force counting kmers longer than the order of the index" << std::endl;
    std::cerr << "  -k N  Set the length of the kmers to N (default " << DEFAULT_K << ")" << std::endl;
    std::cerr << "  -N    Include kmers containing Ns" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  int c = 0;
  size_type k = DEFAULT_K;
  bool force = false, include_Ns = false;
  while((c = getopt(argc, argv, "fk:N")) != -1)
  {
    switch(c)
    {
    case 'f':
      force = true; break;
    case 'k':
      k = std::stoul(optarg); break;
    case 'N':
      include_Ns = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  if(optind >= argc)
  {
    std::cerr << "count_kmers: Base name not specified" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  bool compare = (optind + 1 < argc);
  std::string left_name = argv[optind];
  std::string right_name = (compare ? argv[optind + 1] : "");

  std::cout << "Kmer counter" << std::endl;
  std::cout << std::endl;
  if(compare)
  {
    std::cout << "Left name:   " << left_name << std::endl;
    std::cout << "Right name:  " << right_name << std::endl;
  }
  else
  {
    std::cout << "Base name:   " << left_name << std::endl;
  }
  std::cout << "K:           " << k << std::endl;
  std::cout << "Options:    ";
  if(force) { std::cout << " force"; }
  if(include_Ns) { std::cout << " include_Ns"; }
  std::cout << std::endl << std::endl;

  if(compare) { compareKmers(left_name, right_name, k, force, include_Ns); }
  else { countKmers(left_name, k, force, include_Ns); }

  return 0;
}

//------------------------------------------------------------------------------

void
countKmers(const std::string& base_name, size_type k, bool force, bool include_Ns)
{
  GCSA index;
  sdsl::load_from_file(index, base_name + GCSA::EXTENSION);
  std::cout << "GCSA:        " << index.size() << " paths, order " << index.order() << std::endl;

  double start = readTimer();
  size_type kmer_count = countKMers(index, k, include_Ns, force);
  double seconds = readTimer() - start;
  std::cout << "Kmers:       " << kmer_count << std::endl;
  std::cout << std::endl;

  std::cout << "Kmers counted in " << seconds << " seconds (" << (kmer_count / seconds) << " / s)" << std::endl;
  std::cout << std::endl;
}

void
compareKmers(const std::string& left_name, const std::string& right_name, size_type k, bool force, bool include_Ns)
{
  GCSA left, right;
  sdsl::load_from_file(left, left_name + GCSA::EXTENSION);
  std::cout << "Left:        " << left.size() << " paths, order " << left.order() << std::endl;
  sdsl::load_from_file(right, right_name + GCSA::EXTENSION);
  std::cout << "Right:       " << right.size() << " paths, order " << right.order() << std::endl;

  double start = readTimer();
  std::array<size_type, 3> results = countKMers(left, right, k, include_Ns, force);
  double seconds = readTimer() - start;
  std::cout << "Shared:      " << results[0] << " kmers" << std::endl;
  std::cout << "Left:        " << results[1] << " unique kmers" << std::endl;
  std::cout << "Right:       " << results[2] << " unique kmers" << std::endl;
  std::cout << std::endl;

  std::cout << "Kmers counted in " << seconds << " seconds ("
            << ((results[0] + results[1] + results[2]) / seconds) << " / s)" << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
