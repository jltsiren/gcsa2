/*
  Copyright (c) 2018, 2019 Jouni Siren
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

#include <gcsa/algorithms.h>

#include <stack>

#include <gcsa/internal.h>

namespace gcsa
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type KMerSearchParameters::SEED_LENGTH;

//------------------------------------------------------------------------------

// Other class variables.

const std::string KMerSearchParameters::LEFT_EXTENSION = ".left";
const std::string KMerSearchParameters::RIGHT_EXTENSION = ".right";

//------------------------------------------------------------------------------

constexpr size_type MAX_ERRORS = 100; // Suppress further error messages.

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

const size_type RANDOM_LOCATE_SIZE = 10;  // Number of occurrences to locate.

bool
verifyIndex(const GCSA& index, const LCPArray* lcp, const InputGraph& graph)
{
  std::vector<KMer> kmers; graph.read(kmers);
  return verifyIndex(index, lcp, kmers, graph.k(), graph.mapping);
}

bool
verifyIndex(const GCSA& index, const LCPArray* lcp, std::vector<KMer>& kmers, size_type kmer_length, const NodeMapping& mapping)
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
      Node::map(expected, mapping);
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

        // Random occurrences.
        std::vector<node_type> random_occs;
        index.locate(range, RANDOM_LOCATE_SIZE, random_occs);
        size_type expected_occs = std::min(RANDOM_LOCATE_SIZE, static_cast<size_type>(occs.size()));
        if(random_occs.size() != expected_occs)
        {
          #pragma omp critical
          {
            if(printFailure(fails))
            {
              std::cerr << "verifyIndex(): locate(" << kmer << ") failed: Expected "
                        << expected_occs << " random occurrences, got " << random_occs.size() << std::endl;
            }
          }
        }
        else
        {
          bool is_subset = true;
          for(size_type i = 0, j = 0; i < random_occs.size(); i++)
          {
            while(occs[j] < random_occs[i] && j + 1 < occs.size()) { j++; }
            if(random_occs[i] != occs[j]) { is_subset = false; break; }
            j++;
          }
          if(!is_subset)
          {
            #pragma omp critical
            {
              if(printFailure(fails))
              {
                std::cerr << "verifyIndex(): locate(" << kmer << ") failed: ";
                printOccs(random_occs, std::cerr);
                std::cerr << " is not a subset of ";
                printOccs(occs, std::cerr);
                std::cerr << std::endl;
              }
            }
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

template<class State>
struct KMerSeedCollector
{
  size_type          k;
  bool               include_Ns;
  std::vector<State> seeds;

  KMerSeedCollector(size_type depth, bool Ns) : k(depth), include_Ns(Ns), seeds() {}

  inline bool allChars() const { return this->include_Ns; }

  inline bool reportCondition(const State& state) const
  {
    return (state.k == this->k);
  }

  inline bool expandCondition(const State& state) const
  {
    return (state.k < this->k);
  }

  inline void report(const State& state)
  {
    seeds.push_back(state);
  }
};

//------------------------------------------------------------------------------

struct KMerSearchState
{
  range_type range;
  size_type  k;

  KMerSearchState() : range(0, 0), k(0) {}
  explicit KMerSearchState(range_type rng) : range(rng), k(0) {}
  KMerSearchState(range_type rng, const KMerSearchState& successor) : range(rng), k(successor.k + 1) {}
};

struct KMerCounter
{
  size_type k;
  bool      include_Ns;
  size_type count;

  KMerCounter(size_type depth, bool Ns) : k(depth), include_Ns(Ns), count(0) {}

  inline bool allChars() const { return this->include_Ns; }

  inline bool reportCondition(const KMerSearchState& state) const
  {
    return (state.k == this->k);
  }

  inline bool expandCondition(const KMerSearchState& state) const
  {
    return (state.k < this->k);
  }

  inline void report(const KMerSearchState&)
  {
    this->count++;
  }
};

template<class Handler>
void
processSubtree(const GCSA& index, std::stack<KMerSearchState>& state_stack, Handler& handler)
{
  std::vector<range_type> pred(index.alpha.sigma);
  size_type limit = (handler.allChars() ? index.alpha.sigma : index.alpha.fast_chars + 2);
  while(!(state_stack.empty()))
  {
    KMerSearchState curr = state_stack.top(); state_stack.pop();
    if(Range::empty(curr.range)) { continue; }
    if(handler.reportCondition(curr)) { handler.report(curr); }
    if(handler.expandCondition(curr))
    {
      if(handler.allChars()) { index.LF_all(curr.range, pred); }
      else { index.LF_fast(curr.range, pred); }
      for(size_type comp = 1; comp + 1 < limit; comp++)
      {
        state_stack.push(KMerSearchState(pred[comp], curr));
      }
    }
  }
}

size_type
countKMers(const GCSA& index, size_type k, const KMerSearchParameters& parameters)
{
  if(k == 0) { return 1; }
  if(k > index.order() && !parameters.force)
  {
    std::cerr << "countKMers(): The value of k is greater than the order of the index" << std::endl;
    return 0;
  }

  // Create an array of seed kmers.
  std::vector<KMerSearchState> seeds;
  {
    std::stack<KMerSearchState> state_stack;
    state_stack.push(KMerSearchState(range_type(0, index.size() - 1)));
    KMerSeedCollector<KMerSearchState> collector(std::min(k, parameters.seed_length), parameters.include_Ns);
    processSubtree(index, state_stack, collector);
    seeds = collector.seeds;
  }

  // Extend the seeds in parallel.
  size_type result = 0;
  #pragma omp parallel for schedule (dynamic, 1)
  for(size_type i = 0; i < seeds.size(); i++)
  {
    std::stack<KMerSearchState> state_stack;
    state_stack.push(seeds[i]);
    KMerCounter counter(k, parameters.include_Ns);
    processSubtree(index, state_stack, counter);
    #pragma omp atomic
    result += counter.count;
  }

  return result;
}

//------------------------------------------------------------------------------

struct KMerComparisonState
{
  range_type               left, right;
  size_type                k;
  std::array<size_type, 3> kmer;

  KMerComparisonState() : left(0, 0), right(0, 0), k(0), kmer({ 0, 0, 0 }) {}

  KMerComparisonState(range_type left_range, range_type right_range) :
    left(left_range), right(right_range), k(0), kmer({ 0, 0, 0 })
  {
  }

  KMerComparisonState(range_type left_range, range_type right_range,
    const KMerComparisonState& successor, size_type comp) :
    left(left_range), right(right_range), k(successor.k + 1), kmer(successor.kmer)
  {
    this->set(successor.k, comp);
  }

  constexpr static size_type MAX_K     = 64;
  constexpr static size_type COMP_BITS = 3;
  constexpr static size_type OFFSET    = 6;
  constexpr static size_type MASK      = 0x3F;

  inline void set(size_type i, size_type comp)
  {
    size_type offset = (i * COMP_BITS) >> OFFSET;
    size_type bit = (i * COMP_BITS) & MASK;
    this->kmer[offset] |= comp << bit;
    this->kmer[offset + 1] |= comp >> (WORD_BITS - bit);  // This is |= 0 if the comp does not overflow.
  }
};

constexpr size_type KMerComparisonState::MAX_K;
constexpr size_type KMerComparisonState::COMP_BITS;
constexpr size_type KMerComparisonState::OFFSET;
constexpr size_type KMerComparisonState::MASK;

struct KMerSymmetricDifference
{
  size_type k;
  bool      include_Ns;
  size_type shared, left, right;

  std::vector<KMerComparisonState> left_kmers, right_kmers;

  KMerSymmetricDifference(size_type depth, bool Ns) :
    k(depth), include_Ns(Ns),
    shared(0), left(0), right(0),
    left_kmers(), right_kmers()
  {
  }

  inline bool allChars() const { return this->include_Ns; }

  inline bool reportCondition(const KMerComparisonState& state) const
  {
    return (state.k == this->k);
  }

  inline bool expandCondition(const KMerComparisonState& state) const
  {
    return (state.k < this->k);
  }

  inline void report(const KMerComparisonState& state)
  {
    if(Range::length(state.left) > 0 && Range::length(state.right) > 0) { this->shared++; }
    else if(Range::length(state.left) > 0 && Range::length(state.right) == 0)
    {
      this->left_kmers.push_back(state); this->left++;
    }
    else if(Range::length(state.left) == 0 && Range::length(state.right) > 0)
    {
      this->right_kmers.push_back(state); this->right++;
    }
  }
};

template<class Handler>
void
processSubtrees(const GCSA& left, const GCSA& right, std::stack<KMerComparisonState>& state_stack, Handler& handler)
{
  std::vector<range_type> left_pred(left.alpha.sigma), right_pred(right.alpha.sigma);
  size_type limit = (handler.allChars() ? left.alpha.sigma : left.alpha.fast_chars + 2);
  while(!(state_stack.empty()))
  {
    KMerComparisonState curr = state_stack.top(); state_stack.pop();
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
        state_stack.push(KMerComparisonState(left_pred[comp], right_pred[comp], curr, comp));
      }
    }
  }
}

std::array<size_type, 3>
compareKMers(const GCSA& left, const GCSA& right, size_type k, const KMerSearchParameters& parameters)
{
  std::array<size_type, 3> result = { 0, 0, 0 };

  if(k == 0) { result[0] = 1; return result; }
  if(k > left.order() && !parameters.force)
  {
    std::cerr << "compareKMers(): The value of k is greater than the order of the left index" << std::endl;
    return result;
  }
  if(k > right.order() && !parameters.force)
  {
    std::cerr << "compareKMers(): The value of k is greater than the order of the right index" << std::endl;
    return result;
  }
  if(k > KMerComparisonState::MAX_K)
  {
    std::cerr << "compareKMers(): Comparison is only supported for k <= " << KMerComparisonState::MAX_K << std::endl;
    return result;
  }

  if(left.alpha.sigma != right.alpha.sigma || left.alpha.fast_chars != right.alpha.fast_chars)
  {
    std::cerr << "compareKMers(): The indexes use incompatible alphabets" << std::endl;
    return result;
  }

  std::ofstream left_output, right_output;
  bool write = !(parameters.output.empty());
  if(write)
  {
    std::string left_name = parameters.output + KMerSearchParameters::LEFT_EXTENSION;
    left_output.open(left_name.c_str(), std::ios_base::binary);
    if(!left_output)
    {
      std::cerr << "compareKMers(): Cannot open output file " << left_name << std::endl;
      return result;
    }
    std::string right_name = parameters.output + KMerSearchParameters::RIGHT_EXTENSION;
    right_output.open(right_name.c_str(), std::ios_base::binary);
    if(!right_output)
    {
      std::cerr << "compareKMers(): Cannot open output file " << right_name << std::endl;
      left_output.close();
      return result;
    }
  }

  // Create an array of seed kmers.
  std::vector<KMerComparisonState> seeds;
  {
    std::stack<KMerComparisonState> state_stack;
    state_stack.push(KMerComparisonState(range_type(0, left.size() - 1), range_type(0, right.size() - 1)));
    KMerSeedCollector<KMerComparisonState> collector(std::min(k, parameters.seed_length), parameters.include_Ns);
    processSubtrees(left, right, state_stack, collector);
    seeds = collector.seeds;
  }

  // Extend the seeds in parallel.
  #pragma omp parallel for schedule (dynamic, 1)
  for(size_type i = 0; i < seeds.size(); i++)
  {
    std::stack<KMerComparisonState> state_stack;
    state_stack.push(seeds[i]);
    KMerSymmetricDifference counter(k, parameters.include_Ns);
    processSubtrees(left, right, state_stack, counter);
    #pragma omp critical
    {
      result[0] += counter.shared;
      result[1] += counter.left;
      result[2] += counter.right;
      if(write)
      {
        DiskIO::write(left_output, counter.left_kmers.data(), counter.left);
        DiskIO::write(right_output, counter.right_kmers.data(), counter.right);
      }
    }
  }

  if(write) { left_output.close(); right_output.close(); }
  return result;
}

//------------------------------------------------------------------------------

void
printStatistics(const GCSA& gcsa, const LCPArray& lcp_array)
{
  printHeader("Paths"); std::cout << gcsa.size() << std::endl;
  printHeader("Edges"); std::cout << gcsa.edgeCount() << std::endl;
  printHeader("Samples");
  std::cout << gcsa.sampleCount() << " (at " << gcsa.sampledPositions() << " positions, "
            << gcsa.sampleBits() << " bits each)" << std::endl;
  printHeader("Max query"); std::cout << gcsa.order() << std::endl;
  std::cout << std::endl;

  size_type index_bytes = sdsl::size_in_bytes(gcsa);
  size_type sample_bytes =
    sdsl::size_in_bytes(gcsa.sampled_paths) + sdsl::size_in_bytes(gcsa.sampled_path_rank) +
    sdsl::size_in_bytes(gcsa.stored_samples) +
    sdsl::size_in_bytes(gcsa.samples) + sdsl::size_in_bytes(gcsa.sample_select);
  size_type counter_bytes =
    sdsl::size_in_bytes(gcsa.extra_pointers) + sdsl::size_in_bytes(gcsa.redundant_pointers);
  size_type lcp_bytes = sdsl::size_in_bytes(lcp_array);
  printHeader("Core index"); std::cout << inMegabytes(index_bytes - sample_bytes - counter_bytes) << " MB" << std::endl;
  printHeader("Samples"); std::cout << inMegabytes(sample_bytes) << " MB" << std::endl;
  printHeader("Counter"); std::cout << inMegabytes(counter_bytes) << " MB" << std::endl;
  printHeader("LCP array"); std::cout << inMegabytes(lcp_bytes) << " MB" << std::endl;
  printHeader("Total size"); std::cout << inMegabytes(index_bytes + lcp_bytes) << " MB" << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

} // namespace gcsa
