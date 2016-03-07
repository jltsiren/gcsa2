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

#include "algorithms.h"

using namespace gcsa;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "usage: count_kmers base_name k" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  size_type k = std::stoul(argv[2]);
  std::cout << "Kmer counter" << std::endl;
  std::cout << std::endl;
  std::cout << "Base name:  " << base_name << std::endl;
  std::cout << "K:          " << k << std::endl;
  std::cout << std::endl;

  GCSA index;
  std::string gcsa_name = base_name + GCSA::EXTENSION;
  sdsl::load_from_file(index, gcsa_name);
  std::cout << "GCSA:       " << index.size() << " paths, order " << index.order() << std::endl;

  double start = readTimer();
  size_type kmer_count = countKMers(index, k);
  double seconds = readTimer() - start;
  std::cout << "Kmers:      " << kmer_count << std::endl;
  std::cout << std::endl;

  std::cout << "Kmers counted in " << seconds << " seconds (" << (kmer_count / seconds) << " / s)" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
