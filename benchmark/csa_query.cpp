/*
  Copyright (c) 2019 Jouni Siren
  Copyright (c) 2016, 2017 Genome Research Ltd.

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

typedef sdsl::csa_wt<sdsl::wt_huff<>, 17>                      StandardCSA;
typedef sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<>>, 17> FastCSA;

//------------------------------------------------------------------------------

size_type filter(std::vector<std::string>& patterns);

template<class IndexType>
void benchmark(const std::string& file_name, const std::string& header,
  const std::vector<std::string>& patterns, size_type pattern_total);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    Version::print(std::cerr, "CSA query benchmark");
    std::cerr << "usage: csa_query base_name patterns" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  std::string pattern_name = argv[2];
  Version::print(std::cout, "CSA query benchmark");
  printHeader("Base name"); std::cout << base_name << std::endl;
  printHeader("Pattern file"); std::cout << pattern_name << std::endl;
  std::cout << std::endl;

  std::vector<std::string> patterns;
  size_type pattern_total = readRows(pattern_name, patterns, true);
  pattern_total -= filter(patterns);
  printHeader("Patterns");
  std::cout << patterns.size() << " (total " << inMegabytes(pattern_total) << " MB)" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  std::string csa_name = base_name + ".csa", csa_header = "SSA";
  benchmark<StandardCSA>(csa_name, csa_header, patterns, pattern_total);

  std::string fast_name = base_name + ".fcsa", fast_header = "SSA-fast";
  benchmark<FastCSA>(fast_name, fast_header, patterns, pattern_total);

  return 0;
}

//------------------------------------------------------------------------------

size_type
filter(std::vector<std::string>& patterns)
{
  size_type tail = 0, filtered = 0;
  for(size_type i = 0; i < patterns.size(); i++)
  {
    const std::string& curr = patterns[i];
    bool ok = false;
    for(size_type j = 0; j < curr.length(); j++)
    {
      if(curr[j] != 'N') { ok = true; break; }
    }
    if(ok) { patterns[tail] = curr; tail++; }
    else { filtered += curr.length(); }
  }

  patterns.resize(tail);
  return filtered;
}

//------------------------------------------------------------------------------

template<class IndexType>
void
benchmark(const std::string& file_name, const std::string& header,
  const std::vector<std::string>& patterns, size_type pattern_total)
{
  IndexType index;
  if(!sdsl::load_from_file(index, file_name))
  {
    std::cerr << "benchmark(): Cannot load the index from " << file_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  printHeader(header);
  std::cout << inMegabytes(sdsl::size_in_bytes(index)) << " MB ("
            << inBPC(sdsl::size_in_bytes(index), index.size()) << " bpc)" << std::endl;
  std::cout << std::endl;

  std::vector<range_type> ranges; ranges.reserve(patterns.size());
  {
    double start = readTimer();
    size_type total = 0;
    for(size_type i = 0; i < patterns.size(); i++)
    {
      range_type temp;
      sdsl::backward_search(index, 0, index.size() - 1, patterns[i].begin(), patterns[i].end(), temp.first, temp.second);
      if(!Range::empty(temp)) { ranges.push_back(temp); }
      total += Range::length(temp);
    }
    double seconds = readTimer() - start;
    printTime("find()", patterns.size(), seconds);
    printHeader("find()");
    std::cout << "Found " << ranges.size() << " patterns matching " << total << " positions ("
              << (inMegabytes(pattern_total) / seconds) << " MB/s)" << std::endl;
    std::cout << std::endl;
  }

  {
    double start = readTimer();
    std::vector<size_type> results;
    size_type total = 0;
    for(size_type i = 0; i < ranges.size(); i++)
    {
      for(size_type j = ranges[i].first; j <= ranges[i].second; j++) { results.push_back(index[j]); }
      total += results.size();
      results.clear();
    }
    double seconds = readTimer() - start;
    printTime("locate()", ranges.size(), seconds);
    printHeader("locate()");
    std::cout << total << " occurrences (" <<
              inMicroseconds(seconds / total) << " µs/occurrence)" << std::endl;
    std::cout << std::endl;
  }

  std::cout << std::endl;
}

//------------------------------------------------------------------------------
