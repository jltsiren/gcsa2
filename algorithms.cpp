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

const size_type MAX_ERRORS = 100; // Suppress further error messages.

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
  if(failure_count == MAX_ERRORS)
  {
    std::cerr << "verifyIndex(): There were further errors" << std::endl;
  }
  failure_count++;
  return (failure_count <= MAX_ERRORS);
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
    if(Range::empty(bounds[thread])) { continue; }
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

const size_type KMER_SEED_LENGTH = 5; // Create all patterns of that length before parallelizing.

struct ReverseTrieNode
{
  range_type range;
  size_type  depth;

  ReverseTrieNode() : range(0, 0), depth(0) {}
  ReverseTrieNode(range_type rng, size_type d) : range(rng), depth(d) {}
};

struct KMerCounter
{
  size_type count, depth_limit;
  bool      all_chars;

  KMerCounter(size_type limit, bool include_Ns) : count(0), depth_limit(limit), all_chars(include_Ns) {}

  inline void add(const ReverseTrieNode& node)
  {
    if(node.depth == this->depth_limit) { this->count++; }
  }
};

struct SeedCollector
{
  size_type depth_limit;
  std::vector<ReverseTrieNode> seeds;
  bool all_chars;

  SeedCollector(size_type limit, bool include_Ns) : depth_limit(limit), seeds(), all_chars(include_Ns) {}

  inline void add(const ReverseTrieNode& node)
  {
    if(node.depth >= this->depth_limit) { this->seeds.push_back(node); }
  }
};

template<class Handler>
void
processSubtree(const GCSA& index, std::stack<ReverseTrieNode>& node_stack, Handler& handler)
{
  std::vector<range_type> predecessors(index.alpha.sigma);
  size_type limit = (handler.all_chars ? index.alpha.sigma : index.alpha.fast_chars + 2);
  while(!(node_stack.empty()))
  {
    ReverseTrieNode curr = node_stack.top(); node_stack.pop();
    if(Range::empty(curr.range)) { continue; }
    handler.add(curr);
    if(curr.depth < handler.depth_limit)
    {
      if(handler.all_chars) { index.LF_all(curr.range, predecessors); }
      else { index.LF_fast(curr.range, predecessors); }
      for(size_type comp = 1; comp + 1 < limit; comp++)
      {
        node_stack.push(ReverseTrieNode(predecessors[comp], curr.depth + 1));
      }
    }
  }
}

size_type
countKMers(const GCSA& index, size_type k, bool include_Ns, bool force)
{
  if(k == 0) { return 1; }
  if(k > index.order() && !force)
  {
    std::cerr << "countKMers(): The value of k is greater than the order of the index" << std::endl;
    return 0;
  }

  // Create an array of seed kmers of length KMER_SEED_LENGTH.
  std::vector<ReverseTrieNode> seeds;
  {
    std::stack<ReverseTrieNode> node_stack;
    size_type limit = (include_Ns ? index.alpha.sigma : index.alpha.fast_chars + 2);
    for(size_type comp = 1; comp + 1 < limit; comp++)
    {
      node_stack.push(ReverseTrieNode(index.charRange(comp), 1));
    }
    SeedCollector collector(std::min(k, KMER_SEED_LENGTH), include_Ns);
    processSubtree(index, node_stack, collector);
    seeds = collector.seeds;
  }

  // Extend the seeds in parallel.
  size_type result = 0;
  #pragma omp parallel for schedule (dynamic, 1)
  for(size_type i = 0; i < seeds.size(); i++)
  {
    std::stack<ReverseTrieNode> node_stack;
    node_stack.push(seeds[i]);
    KMerCounter counter(k, include_Ns);
    processSubtree(index, node_stack, counter);
    #pragma omp atomic
    result += counter.count;
  }

  return result;
}

//------------------------------------------------------------------------------

struct KmerSearchState
{
  range_type  left, right;
  size_type   k;

  KmerSearchState() : left(0, 0), right(0, 0), k(0) {}

  KmerSearchState(range_type left_range, range_type right_range) :
    left(left_range), right(right_range), k(0)
  {
  }

  KmerSearchState(range_type left_range, range_type right_range, const KmerSearchState& successor) :
    left(left_range), right(right_range), k(successor.k + 1)
  {
  }
};

struct KmerSymmetricDifference
{
  size_type k;
  bool      include_Ns;
  size_type shared, left, right;

  KmerSymmetricDifference(size_type depth, bool Ns) :
    k(depth), include_Ns(Ns),
    shared(0), left(0), right(0)
  {
  }

  inline bool allChars() const { return this->include_Ns; }

  inline bool reportCondition(const KmerSearchState& state) const
  {
    return (state.k == this->k);
  }

  inline bool expandCondition(const KmerSearchState& state) const
  {
    return (state.k < this->k);
  }

  inline void report(const KmerSearchState& state)
  {
    if(Range::length(state.left) > 0 && Range::length(state.right) > 0) { this->shared++; }
    else if(Range::length(state.left) > 0 && Range::length(state.right) == 0) { this->left++; }
    else if(Range::length(state.left) == 0 && Range::length(state.right) > 0) { this->right++; }
  }
};

struct KmerSeedCollector
{
  size_type k;
  bool include_Ns;
  std::vector<KmerSearchState> seeds;

  KmerSeedCollector(size_type depth, bool Ns) : k(depth), include_Ns(Ns), seeds() {}

  inline bool allChars() const { return this->include_Ns; }

  inline bool reportCondition(const KmerSearchState& state) const
  {
    return (state.k == this->k);
  }

  inline bool expandCondition(const KmerSearchState& state) const
  {
    return (state.k < this->k);
  }

  inline void report(const KmerSearchState& state)
  {
    seeds.push_back(state);
  }
};

template<class Handler>
void
processSubtrees(const GCSA& left, const GCSA& right, std::stack<KmerSearchState>& state_stack, Handler& handler)
{
  std::vector<range_type> left_pred(left.alpha.sigma), right_pred(right.alpha.sigma);
  size_type limit = (handler.allChars() ? left.alpha.sigma : left.alpha.fast_chars + 2);
  while(!(state_stack.empty()))
  {
    KmerSearchState curr = state_stack.top(); state_stack.pop();
    if(Range::empty(curr.left) && Range::empty(curr.right)) { continue; }
    if(handler.reportCondition(curr)) { handler.report(curr); }
    if(handler.expandCondition(curr))
    {
      if(handler.allChars())
      {
        left.LF_all(curr.left, left_pred); right.LF_all(curr.right, right_pred);
      }
      else
      {
        left.LF_fast(curr.left, left_pred); right.LF_fast(curr.right, right_pred);
      }
      for(size_type comp = 1; comp + 1 < limit; comp++)
      {
        state_stack.push(KmerSearchState(left_pred[comp], right_pred[comp], curr));
      }
    }
  }
}

std::array<size_type, 3>
countKMers(const GCSA& left, const GCSA& right, size_type k, bool include_Ns, bool force)
{
  std::array<size_type, 3> result = { 0, 0, 0 };

  if(k == 0) { result[0] = 1; return result; }
  if(k > left.order() && !force)
  {
    std::cerr << "countKMers(): The value of k is greater than the order of the left index" << std::endl;
    return result;
  }
  if(k > right.order() && !force)
  {
    std::cerr << "countKMers(): The value of k is greater than the order of the right index" << std::endl;
    return result;
  }

  if(left.alpha.sigma != right.alpha.sigma || left.alpha.fast_chars != right.alpha.fast_chars)
  {
    std::cerr << "countKMers(): The indexes use incompatible alphabets" << std::endl;
    return result;
  }

  // Create an array of seed kmers of length KMER_SEED_LENGTH.
  std::vector<KmerSearchState> seeds;
  {
    std::stack<KmerSearchState> state_stack;
    state_stack.push(KmerSearchState(range_type(0, left.size() - 1), range_type(0, right.size() - 1)));
    KmerSeedCollector collector(std::min(k, KMER_SEED_LENGTH), include_Ns);
    processSubtrees(left, right, state_stack, collector);
    seeds = collector.seeds;
  }

  // Extend the seeds in parallel.
  #pragma omp parallel for schedule (dynamic, 1)
  for(size_type i = 0; i < seeds.size(); i++)
  {
    std::stack<KmerSearchState> state_stack;
    state_stack.push(seeds[i]);
    KmerSymmetricDifference counter(k, include_Ns);
    processSubtrees(left, right, state_stack, counter);
    #pragma omp critical
    {
      result[0] += counter.shared;
      result[1] += counter.left;
      result[2] += counter.right;
    }
  }

  return result;
}

//------------------------------------------------------------------------------

} // namespace gcsa
