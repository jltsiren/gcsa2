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

#include <gcsa/utils.h>
#include <sdsl/suffix_arrays.hpp>

using namespace gcsa;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: csa_builder base_name" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  std::string index_name = base_name + ".csa";

  double start = readTimer();
  sdsl::csa_wt<sdsl::wt_huff<>, 17> index;
  sdsl::construct(index, base_name, 1);
  sdsl::store_to_file(index, index_name);
  double seconds = readTimer() - start;

  double megabytes = inMegabytes(sdsl::size_in_bytes(index));
  double bpc = sdsl::size_in_bytes(index) * 8.0 / index.size();
  std::cout << "CSA built in " << seconds << " seconds (" << (inMegabytes(index.size()) / seconds) << " MB/s)" << std::endl;
  std::cout << "Index size: " << megabytes << " MB (" << bpc << " bpc)" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
