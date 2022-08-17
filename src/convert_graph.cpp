/*
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

#include <gcsa/files.h>
#include <gcsa/gcsa.h>

using namespace gcsa;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)  // From stdin to stdout.
  {
    std::vector<KMer> kmers;
    size_type kmer_length = readText(std::cin, kmers, Alphabet());
    writeBinary(std::cout, kmers, kmer_length);
  }
  else
  {
    Version::print(std::cout, "GCSA2 input graph converter");
    InputGraph graph(argc - 1, argv + 1, false);
    for(size_type i = 0; i < graph.files(); i++)
    {
      std::string base_name = argv[i + 1];
      std::cout << base_name << InputGraph::TEXT_EXTENSION << " -> "
                << base_name << InputGraph::BINARY_EXTENSION << std::endl;

      std::vector<KMer> kmers;
      graph.read(kmers, i, false);
      writeKMers(base_name, kmers, graph.k());
      std::cout << std::endl;
    }
  }

  return 0;
}

//------------------------------------------------------------------------------
