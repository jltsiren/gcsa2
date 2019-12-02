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

#include <gcsa/algorithms.h>

using namespace gcsa;

//------------------------------------------------------------------------------

size_type filter(std::vector<std::string>& patterns);

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    Version::print(std::cerr, "GCSA2 query benchmark");
    std::cerr << "usage: query_gcsa base_name [patterns]" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  std::string pattern_name;
  if(argc > 2) { pattern_name = argv[2]; }
  Version::print(std::cout, "GCSA2 query benchmark");
  printHeader("Base name"); std::cout << base_name << std::endl;
  if(!(pattern_name.empty())) { printHeader("Pattern file"); std::cout << pattern_name << std::endl; }
  std::cout << std::endl;

  GCSA index;
  std::string gcsa_name = base_name + GCSA::EXTENSION;
  if(!sdsl::load_from_file(index, gcsa_name))
  {
    std::cerr << "query_gcsa: Cannot load the index from " << gcsa_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  LCPArray lcp;
  std::string lcp_name = base_name + LCPArray::EXTENSION;
  if(!sdsl::load_from_file(lcp, lcp_name))
  {
    std::cerr << "query_gcsa: Cannot load the LCP array from " << lcp_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(pattern_name.empty())
  {
    printStatistics(index, lcp);
    std::exit(EXIT_SUCCESS);
  }
  else
  {
    printHeader("GCSA"); std::cout << inMegabytes(sdsl::size_in_bytes(index)) << " MB" << std::endl;
    printHeader("LCP"); std::cout << inMegabytes(sdsl::size_in_bytes(lcp)) << " MB" << std::endl;
  }

  std::vector<std::string> patterns;
  size_type pattern_total = readRows(pattern_name, patterns, true);
  pattern_total -= filter(patterns);
  printHeader("Patterns");
  std::cout << patterns.size() << " (total " << inMegabytes(pattern_total) << " MB)" << std::endl;
  std::cout << std::endl;

  std::vector<range_type> ranges; ranges.reserve(patterns.size());
  std::vector<size_type> lengths; lengths.reserve(patterns.size());
  {
    double start = readTimer();
    size_type total = 0;
    for(size_type i = 0; i < patterns.size(); i++)
    {
      range_type temp = index.find(patterns[i]);
      if(!Range::empty(temp)) { ranges.push_back(temp); lengths.push_back(patterns[i].length()); }
      total += Range::length(temp);
    }
    double seconds = readTimer() - start;
    printTime("find()", patterns.size(), seconds);
    printHeader("find()");
    std::cout << "Found " << ranges.size() << " patterns matching " << total << " paths ("
              << (inMegabytes(pattern_total) / seconds) << " MB/s)" << std::endl;
    std::cout << std::endl;
  }

  std::vector<range_type> parents(ranges.size());
  {
    double start = readTimer();
    size_type total = 0;
    for(size_type i = 0; i < ranges.size(); i++)
    {
      STNode temp = lcp.parent(ranges[i]);
      total += lengths[i] - temp.lcp();
      parents[i] = temp.range();
    }
    double seconds = readTimer() - start;
    printTime("parent()", ranges.size(), seconds);
    printHeader("parent()");
    std::cout << "Average distance " << (total / (double)(ranges.size())) << " characters" << std::endl;
    std::cout << std::endl;
  }

  {
    double start = readTimer();
    size_type total = 0;
    for(size_type i = 0; i < parents.size(); i++)
    {
      total += lengths[i] - lcp.depth(parents[i]);
    }
    double seconds = readTimer() - start;
    printTime("depth()", ranges.size(), seconds);
    printHeader("depth()");
    std::cout << "Average distance " << (total / (double)(parents.size())) << " characters" << std::endl;
    std::cout << std::endl;
  }

  std::vector<size_type> counts(ranges.size());
  {
    double start = readTimer();
    size_type total = 0;
    for(size_type i = 0; i < ranges.size(); i++)
    {
      counts[i] = index.count(ranges[i]);
      total += counts[i];
    }
    double seconds = readTimer() - start;
    printTime("count()", ranges.size(), seconds);
    printHeader("count()");
    std::cout << total << " occurrences" << std::endl;
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
    std::cout << total << " occurrences (" <<
              inMicroseconds(seconds / total) << " µs/occurrence)" << std::endl;
    std::cout << std::endl;
  }

  for(size_type i = 0; i < counts.size(); i++)
  {
    if(counts[i] != 0)
    {
      std::cout << "Warning: count() and locate() returned inconsistent results" << std::endl;
      std::cout << std::endl;
      break;
    }
  }

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
