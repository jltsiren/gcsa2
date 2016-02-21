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

#include "gcsa.h"

using namespace gcsa;

//------------------------------------------------------------------------------

void filter(std::vector<std::string>& patterns);

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "usage: query_gcsa base_name patterns" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  std::string pattern_name = argv[2];
  std::cout << "GCSA2 query benchmark" << std::endl;
  std::cout << std::endl;
  printHeader("Base name"); std::cout << base_name << std::endl;
  printHeader("Pattern file"); std::cout << pattern_name << std::endl;
  std::cout << std::endl;

  GCSA index;
  std::string gcsa_name = base_name + GCSA::EXTENSION;
  sdsl::load_from_file(index, gcsa_name);
  printHeader("GCSA"); std::cout << inMegabytes(sdsl::size_in_bytes(index)) << " MB" << std::endl;

  std::vector<std::string> patterns;
  size_type pattern_total = readRows(pattern_name, patterns, true);
  filter(patterns);
  printHeader("Patterns");
  std::cout << patterns.size() << " (total " << inMegabytes(pattern_total) << " MB)" << std::endl;
  std::cout << std::endl;

  std::vector<range_type> ranges; ranges.reserve(patterns.size());
  {
    double start = readTimer();
    for(size_type i = 0; i < patterns.size(); i++)
    {
      range_type temp = index.find(patterns[i]);
      if(!Range::empty(temp)) { ranges.push_back(temp); }
    }
    double seconds = readTimer() - start;
    printTime("find()", patterns.size(), seconds);
    printHeader("find()");
    std::cout << "Found " << ranges.size() << " patterns ("
              << (inMegabytes(pattern_total) / seconds) << " MB/s)" << std::endl;
    std::cout << std::endl;
  }

  std::vector<size_type> counts(ranges.size());
  {
    double start = readTimer();
    for(size_type i = 0; i < ranges.size(); i++)
    {
      counts[i] = index.count(ranges[i]);
    }
    double seconds = readTimer() - start;
    printTime("count()", ranges.size(), seconds);
    std::cout << std::endl;
  }

  {
    double start = readTimer();
    std::vector<node_type> results;
    size_type total = 0;
    for(size_type i = 0; i < ranges.size(); i++)
    {
      index.locate(ranges[i], results);
      counts[i] -= results.size();
      total += results.size();
    }
    double seconds = readTimer() - start;
    printTime("locate()", ranges.size(), seconds);
    printHeader("locate()");
    std::cout << "Found " << total << " occurrences (" <<
              (total / seconds) << " / second)" << std::endl;
    std::cout << std::endl;
  }

  bool ok = true;
  for(size_type i = 0; i < counts.size(); i++)
  {
    if(counts[i] != 0) { ok = false; }
  }
  if(!ok)
  {
    std::cout << "Warning: count() and locate() returned inconsistent results" << std::endl;
    std::cout << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

void
filter(std::vector<std::string>& patterns)
{
  size_type tail = 0;
  for(size_type i = 0; i < patterns.size(); i++)
  {
    const std::string& curr = patterns[i];
    for(size_type j = 0; j < curr.length(); j++)
    {
      if(curr[j] != 'N') { patterns[tail] = curr; tail++; break; }
    }
  }
  patterns.resize(tail);
}

//------------------------------------------------------------------------------
