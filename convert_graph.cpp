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

#include "files.h"
#include "gcsa.h"

using namespace gcsa;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)  // From stdin to stdout.
  {
    std::vector<KMer> kmers;
    size_type kmer_length = readText(std::cin, kmers);
    writeBinary(std::cout, kmers, kmer_length);
  }
  else
  {
    std::cout << "GCSA input converter" << std::endl;
    std::cout << std::endl;
    for(int i = 1; i < argc; i++)
    {
      std::string base_name = argv[i];
      std::cout << base_name << TEXT_EXTENSION << " -> " << base_name << BINARY_EXTENSION << std::endl;
      std::vector<KMer> kmers;
      size_type kmer_length = readKMers(1, argv + i, kmers, false);
      writeKMers(base_name, kmers, kmer_length);
      std::cout << std::endl;
    }
  }

  return 0;
}

//------------------------------------------------------------------------------
