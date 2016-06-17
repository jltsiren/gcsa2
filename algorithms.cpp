/*
  Copyright (c) 2015, 2016 Genome Research Ltd.

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

#include <stack>

#include "algorithms.h"

namespace gcsa
{

//------------------------------------------------------------------------------

std::ostream&
printOccs(const std::vector<node_type>& occs, std::ostream& out)
{
  out << "{";
  for(size_type i = 0; i < occs.size(); i++)
  {
    out << (i == 0 ? " " : ", ") << Node::decode(occs[i]);
  }
  out << " }";
  return out;
}

void
locateFailure(const std::vector<node_type>& expected, const std::vector<node_type>& occs)
{
  std::cerr << "verifyIndex(): Expected ";
  printOccs(expected, std::cerr) << std::endl;
  std::cerr << "verifyIndex(): Got ";
  printOccs(occs, std::cerr) << std::endl;
}

bool
printFailure(size_type& failure_count)
{
  if(failure_count == GCSA::MAX_ERRORS)
  {
    std::cerr << "verifyIndex(): There were further errors" << std::endl;
  }
  failure_count++;
  return (failure_count <= GCSA::MAX_ERRORS);
}

struct KMerSplitComparator
{
  inline bool operator() (const KMer& left, const KMer& right) const
    {
      return (Key::label(left.key) != Key::label(right.key));
    }
};

bool
verifyIndex(const GCSA& index, const LCPArray* lcp, const InputGraph& graph)
{
  std::vector<KMer> kmers; graph.read(kmers);
  return verifyIndex(index, lcp, kmers, graph.k());
}

bool
verifyIndex(const GCSA& index, const LCPArray* lcp, std::vector<KMer>& kmers, size_type kmer_length)
{
  double start = readTimer();

  size_type threads = omp_get_max_threads();
  parallelQuickSort(kmers.begin(), kmers.end());
  KMerSplitComparator k_comp;
  std::vector<range_type> bounds = getBounds(kmers, threads, k_comp);
  assert(bounds.size() == threads);

  size_type fails = 0, unique = 0;
  #pragma omp parallel for schedule(static)
  for(size_type thread = 0; thread < threads; thread++)
  {
    size_type i = bounds[thread].first;
    while(i <= bounds[thread].second)
    {
      size_type next = i + 1;
      while(next <= bounds[thread].second && Key::label(kmers[next].key) == Key::label(kmers[i].key)) {
          next++;
      }
      #pragma omp atomic
      unique++;

      std::string kmer = Key::decode(kmers[i].key, kmer_length, index.alpha);
      size_type endmarker_pos = kmer.find('$'); // The actual kmer ends at the first endmarker.
      if(endmarker_pos != std::string::npos) { kmer = kmer.substr(0, endmarker_pos + 1); }

      // find()
      range_type range = index.find(kmer);
      if(Range::empty(range))
      {
        #pragma omp critical
        {
          if(printFailure(fails))
          {
            std::cerr << "verifyIndex(): find(" << kmer << ") returned empty range" << std::endl;
          }
        }
        i = next; continue;
      }

      // parent() and depth()
      if(lcp != 0)
      {
        STNode parent_res = lcp->parent(range);
        range_type query_res = range;
        std::string::iterator begin = kmer.begin(), end = kmer.end();
        while(query_res == range)
        {
          --end;
          query_res = index.find(begin, end);
        }
        if(parent_res != query_res || parent_res.lcp() != (size_type)(end - begin))
        {
          #pragma omp critical
          {
            if(printFailure(fails))
            {
              std::cerr << "verifyIndex(): parent" << range << " returned " << parent_res
                        << ", expected " << query_res << " at depth " << (end - begin) << std::endl;
            }
          }
          i = next; continue;
        }
        size_type string_depth = lcp->depth(parent_res.range());
        if(parent_res.lcp() != string_depth)
        {
          #pragma omp critical
          {
            if(printFailure(fails))
            {
              std::cerr << "verifyIndex(): depth" << parent_res.range() << " returned " << string_depth
                        << ", expected " << parent_res.lcp() << std::endl;
            }
          }
          i = next; continue;
        }
      }

      // count()
      std::vector<node_type> expected;
      for(size_type j = i; j < next; j++) { expected.push_back(kmers[j].from); }
      removeDuplicates(expected, false);
      size_type unique_count = index.count(range);
      if(unique_count != expected.size())
      {
        #pragma omp critical
        {
          if(printFailure(fails))
          {
            std::cerr << "verifyIndex(): count" << range << " failed: Expected "
                      << expected.size() << " occurrences, got " << unique_count << std::endl;
          }
        }
        i = next; continue;
      }

      // locate()
      std::vector<node_type> occs;
      index.locate(range, occs);
      if(occs.size() != expected.size())
      {
        #pragma omp critical
        {
          if(printFailure(fails))
          {
            std::cerr << "verifyIndex(): locate(" << kmer << ") failed: Expected "
                      << expected.size() << " occurrences, got " << occs.size() << std::endl;
            locateFailure(expected, occs);
          }
        }
      }
      else
      {
        for(size_type j = 0; j < occs.size(); j++)
        {
          if(occs[j] != expected[j])
          {
            #pragma omp critical
            {
              if(printFailure(fails))
              {
                std::cerr << "verifyIndex(): locate(" << kmer << ") failed: Expected "
                          << Node::decode(expected[j]) << ", got " << Node::decode(occs[j]) << std::endl;
                locateFailure(expected, occs);
              }
            }
            break;
          }
        }
      }

      i = next;
    }
  }

  double seconds = readTimer() - start;
  std::cout << "Queried the index with " << unique << " patterns in " << seconds << " seconds ("
            << (unique / seconds) << " patterns / second)" << std::endl;
  if(fails == 0)
  {
    std::cout << "Index verification complete" << std::endl;
  }
  else
  {
    std::cout << "Index verification failed for " << fails << " patterns" << std::endl;
  }
  std::cout << std::endl;

  return fails == 0;
}

//------------------------------------------------------------------------------

struct ReverseTrieNode
{
  range_type range;
  size_type  depth;
  bool       ends_at_sink;

  const static size_type SEED_LENGTH = 5; // Create all patterns of that length before parallelizing.

  ReverseTrieNode() : range(0, 0), depth(0), ends_at_sink(false) {}
  ReverseTrieNode(range_type rng, size_type d, bool sink) : range(rng), depth(d), ends_at_sink(sink) {}
};

struct KMerCounter
{
  size_type count, limit;

  KMerCounter(size_type _limit) : count(0), limit(_limit) {}

  inline void add(const ReverseTrieNode& node)
  {
    if(node.ends_at_sink || node.depth >= this->limit) { this->count++; }
  }
};

struct SeedCollector
{
  size_type count, limit;
  std::vector<ReverseTrieNode> seeds;

  SeedCollector(size_type _limit) : count(0), limit(_limit), seeds() {}

  inline void add(const ReverseTrieNode& node)
  {
    if(node.depth >= this->limit) { this->seeds.push_back(node); }
    else if(node.ends_at_sink) { this->count++; } // Seeds ending at the sink are counted later.
  }
};

template<class Handler>
void
processSubtree(const GCSA& index, std::stack<ReverseTrieNode>& node_stack, Handler& handler)
{
  std::vector<range_type> predecessors(index.alpha.sigma);
  while(!(node_stack.empty()))
  {
    ReverseTrieNode curr = node_stack.top(); node_stack.pop();
    if(Range::empty(curr.range)) { continue; }
    handler.add(curr);
    if(curr.depth < handler.limit)
    {
      index.LF(curr.range, predecessors);
      for(size_type comp = 1; comp + 1 < index.alpha.sigma; comp++)
      {
        node_stack.push(ReverseTrieNode(predecessors[comp], curr.depth + 1, curr.ends_at_sink));
      }
    }
  }
}

size_type
countKMers(const GCSA& index, size_type k, bool force)
{
  if(k == 0) { return 1; }
  if(k > index.order())
  {
    if(force)
    {
      std::cerr << "countKMers(): Warning: The value of k (" << k
                << ") is greater than the order of the graph (" << index.order() << ")" << std::endl;
    }
    else
    {
      std::cerr << "countKMers(): The value of k (" << k
                << ") is greater than the order of the graph (" << index.order() << ")" << std::endl;
      return 0;
    }
  }

  // Create an array of seed kmers of length ReverseTrieNode::SEED_LENGTH.
  size_type result = 0;
  std::vector<ReverseTrieNode> seeds;
  {
    std::stack<ReverseTrieNode> node_stack;
    node_stack.push(ReverseTrieNode(index.charRange(0), 1, true));
    for(size_type comp = 1; comp + 1 < index.alpha.sigma; comp++)
    {
      node_stack.push(ReverseTrieNode(index.charRange(comp), 1, false));
    }
    size_type seed_length = ReverseTrieNode::SEED_LENGTH;  // avoid direct use of static const
    SeedCollector collector(std::min(k, seed_length));
    processSubtree(index, node_stack, collector);
    result = collector.count;
    seeds = collector.seeds;
  }

  // Extend the seeds in parallel.
  #pragma omp parallel for schedule (dynamic, 1)
  for(size_type i = 0; i < seeds.size(); i++)
  {
    std::stack<ReverseTrieNode> node_stack;
    node_stack.push(seeds[i]);
    KMerCounter counter(k);
    processSubtree(index, node_stack, counter);
    #pragma omp atomic
    result += counter.count;
  }

  return result;
}

//------------------------------------------------------------------------------

} // namespace gcsa
