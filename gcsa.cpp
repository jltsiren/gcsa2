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

#include <cstdio>
#include <cstdlib>

#include "gcsa.h"

namespace gcsa
{

//------------------------------------------------------------------------------

const std::string GCSA::EXTENSION = ".gcsa";

GCSA::GCSA()
{
  this->path_node_count = 0;
  this->max_query_length = 0;
}

GCSA::GCSA(const GCSA& g)
{
  this->copy(g);
}

GCSA::GCSA(GCSA&& g)
{
  *this = std::move(g);
}

GCSA::~GCSA()
{
}

void
GCSA::copy(const GCSA& g)
{
  this->path_node_count = g.path_node_count;
  this->max_query_length = g.max_query_length;

  this->bwt = g.bwt;
  this->alpha = g.alpha;

  this->path_nodes = g.path_nodes;
  this->path_rank = g.path_rank;
  this->path_select = g.path_select;

  this->edges = g.edges;
  this->edge_rank = g.edge_rank;
  this->edge_select = g.edge_select;

  this->sampled_paths = g.sampled_paths;
  this->sampled_path_rank = g.sampled_path_rank;

  this->stored_samples = g.stored_samples;
  this->samples = g.samples;
  this->sample_select = g.sample_select;

  this->setVectors();
}

void
GCSA::swap(GCSA& g)
{
  if(this != &g)
  {
    std::swap(this->path_node_count, g.path_node_count);
    std::swap(this->max_query_length, g.max_query_length);

    this->bwt.swap(g.bwt);
    this->alpha.swap(g.alpha);

    this->path_nodes.swap(g.path_nodes);
    sdsl::util::swap_support(this->path_rank, g.path_rank, &(this->path_nodes), &(g.path_nodes));
    sdsl::util::swap_support(this->path_select, g.path_select, &(this->path_nodes), &(g.path_nodes));

    this->edges.swap(g.edges);
    sdsl::util::swap_support(this->edge_rank, g.edge_rank, &(this->edges), &(g.edges));
    sdsl::util::swap_support(this->edge_select, g.edge_select, &(this->edges), &(g.edges));

    this->sampled_paths.swap(g.sampled_paths);
    sdsl::util::swap_support(this->sampled_path_rank, g.sampled_path_rank, &(this->sampled_paths), &(g.sampled_paths));

    this->stored_samples.swap(g.stored_samples);
    this->samples.swap(g.samples);
    sdsl::util::swap_support(this->sample_select, g.sample_select, &(this->samples), &(g.samples));
  }
}

GCSA&
GCSA::operator=(const GCSA& g)
{
  if(this != &g) { this->copy(g); }
  return *this;
}

GCSA&
GCSA::operator=(GCSA&& g)
{
  if(this != &g)
  {
    this->path_node_count = std::move(g.path_node_count);
    this->max_query_length = std::move(g.max_query_length);

    this->bwt = std::move(g.bwt);
    this->alpha = std::move(g.alpha);

    this->path_nodes = std::move(g.path_nodes);
    this->path_rank = std::move(g.path_rank);
    this->path_select = std::move(g.path_select);

    this->edges = std::move(g.edges);
    this->edge_rank = std::move(g.edge_rank);
    this->edge_select = std::move(g.edge_select);

    this->sampled_paths = std::move(g.sampled_paths);
    this->sampled_path_rank = std::move(g.sampled_path_rank);

    this->stored_samples = std::move(g.stored_samples);
    this->samples = std::move(g.samples);
    this->sample_select = std::move(g.sample_select);

    this->setVectors();
  }
  return *this;
}

void
GCSA::setVectors()
{
  this->path_rank.set_vector(&(this->path_nodes));
  this->path_select.set_vector(&(this->path_nodes));

  this->edge_rank.set_vector(&(this->edges));
  this->edge_select.set_vector(&(this->edges));

  this->sampled_path_rank.set_vector(&(this->sampled_paths));

  this->sample_select.set_vector(&(this->samples));
}

GCSA::size_type
GCSA::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->path_node_count, out, child, "path_node_count");
  written_bytes += sdsl::write_member(this->max_query_length, out, child, "max_query_length");

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->alpha.serialize(out, child, "alpha");

  written_bytes += this->path_nodes.serialize(out, child, "path_nodes");
  written_bytes += this->path_rank.serialize(out, child, "path_rank");
  written_bytes += this->path_select.serialize(out, child, "path_select");

  written_bytes += this->edges.serialize(out, child, "edges");
  written_bytes += this->edge_rank.serialize(out, child, "edge_rank");
  written_bytes += this->edge_select.serialize(out, child, "edge_select");

  written_bytes += this->sampled_paths.serialize(out, child, "sampled_paths");
  written_bytes += this->sampled_path_rank.serialize(out, child, "sampled_path_rank");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSA::load(std::istream& in)
{
  sdsl::read_member(this->path_node_count, in);
  sdsl::read_member(this->max_query_length, in);

  this->bwt.load(in);
  this->alpha.load(in);

  this->path_nodes.load(in);
  this->path_rank.load(in, &(this->path_nodes));
  this->path_select.load(in, &(this->path_nodes));

  this->edges.load(in);
  this->edge_rank.load(in, &(this->edges));
  this->edge_select.load(in, &(this->edges));

  this->sampled_paths.load(in);
  this->sampled_path_rank.load(in, &(this->sampled_paths));

  this->stored_samples.load(in);
  this->samples.load(in);
  this->sample_select.load(in, &(this->samples));
}

//------------------------------------------------------------------------------
/*
void
readPathNodes(const std::string& filename,
  stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels)
{
  sdsl::util::clear(paths); sdsl::util::clear(labels);

  std::ifstream in(filename.c_str(), std::ios_base::binary);
  if(!in)
  {
    std::cerr << "readPathNodes(): Cannot open temporary file " << filename << std::endl;
    std::cerr << "readPathNodes(): Construction aborted" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  size_type path_count = 0, rank_count = 0;
  sdsl::read_member(path_count, in); sdsl::read_member(rank_count, in);
  paths.reserve(path_count + 4); labels.reserve(rank_count + 4);
  for(size_type i = 0; i < path_count; i++) { paths.push_back(PathNode(in, labels)); }
  in.close(); remove(filename.c_str());
}
*/

GCSA::GCSA(std::vector<KMer>& kmers, size_type kmer_length,
  size_type doubling_steps, size_type size_limit, const Alphabet& _alpha)
{
  if(kmers.size() == 0) { return; }
  if(doubling_steps > DOUBLING_STEPS)
  {
    std::cerr << "GCSA::GCSA(): The number of doubling steps is too high: " << doubling_steps << std::endl;
    std::cerr << "GCSA::GCSA(): Reverting the number of doubling steps to " << DOUBLING_STEPS << std::endl;
    doubling_steps = DOUBLING_STEPS;
  }
  if(size_limit > ABSOLUTE_SIZE_LIMIT)
  {
    std::cerr << "GCSA::GCSA(): The size limit is suspiciously high: " << size_limit << " GB" << std::endl;
    std::cerr << "GCSA::GCSA(): Reverting the size limit to " << ABSOLUTE_SIZE_LIMIT << " GB" << std::endl;
    size_limit = ABSOLUTE_SIZE_LIMIT;
  }
  size_type bytes_required =
    2 * sizeof(size_type) + kmers.size() * (sizeof(PathNode) + 2 * sizeof(PathNode::rank_type));
  if(bytes_required > size_limit * GIGABYTE)
  {
    std::cerr << "GCSA::GCSA(): The input is too large: " << (bytes_required / GIGABYTE_DOUBLE) << " GB" << std::endl;
    std::cerr << "GCSA::GCSA(): Construction aborted" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Sort the kmers, build the mapper GCSA for generating the edges.
  stxxl::vector<key_type> keys;
  sdsl::int_vector<0> last_char;
  uniqueKeys(kmers, keys, last_char);
  DeBruijnGraph mapper(keys, kmer_length, _alpha);
  LCP lcp(keys, kmer_length);
  sdsl::util::clear(keys);

  // Transform the kmers into PathNodes.
  stxxl::vector<PathNode> paths;
  stxxl::vector<PathNode::rank_type> labels;
  {
      //std::string temp_file = tempFile(EXTENSION);
      //std::ofstream out(temp_file.c_str(), std::ios_base::binary);
      /*
    if(!out)
    {
      std::cerr << "GCSA::GCSA(): Cannot open temporary file " << temp_file << std::endl;
      std::cerr << "GCSA::GCSA(): Construction aborted" << std::endl;
      std::exit(EXIT_FAILURE);
    }
      */

    stxxl::vector<PathNode::rank_type> temp_labels = PathNode::dummyRankVector();
    //size_type kmer_count = kmers.size();
    //size_type rank_count = 2 * kmer_count;
    for(size_type i = 0; i < kmers.size(); i++)
    {
      PathNode temp(kmers[i], temp_labels);
      temp.setPointer(labels.size());
      paths.push_back(temp);
      for (auto& range : temp_labels) {
          labels.push_back(range);
      }
      //temp.serialize(out, temp_labels);
      temp_labels.resize(0);
    }
    //out.close();

    sdsl::util::clear(kmers);
    //readPathNodes(temp_file, paths, labels);
    //labels.reserve(rank_count + 4);
  }
  //std::cerr << "before prefix doublin with " << paths.size() << " " << labels.size() << std::endl;
  //for (auto& path : paths) { path.print(std::cerr, labels); std::cerr << std::endl; }

  // Build the GCSA in PathNodes.
  this->prefixDoubling(paths, labels, kmer_length, doubling_steps, size_limit, lcp);
  sdsl::util::clear(lcp);
  stxxl::vector<range_type> from_nodes;
  this->mergeByLabel(paths, labels, from_nodes);

  this->build(paths, labels, mapper, last_char);
  this->sample(paths, from_nodes);
}

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
printFailure(const std::string& kmer,
             const std::vector<node_type>& expected, const std::vector<node_type>& occs)
{
  std::cerr << "GCSA::verifyIndex(): locate(" << kmer << ") failed" << std::endl;
  std::cerr << "GCSA::verifyIndex(): Expected ";
  printOccs(expected, std::cerr) << std::endl;
  std::cerr << "GCSA::verifyIndex(): Got ";
  printOccs(occs, std::cerr) << std::endl;
}

struct KMerSplitComparator
{
  inline bool operator() (const KMer& left, const KMer& right) const
    {
      return (Key::label(left.key) != Key::label(right.key));
    }
};

bool
GCSA::verifyIndex(std::vector<KMer>& kmers, size_type kmer_length) const
{
  size_type threads = omp_get_max_threads();
  parallelQuickSort(kmers.begin(), kmers.end());
  KMerSplitComparator k_comp;
  std::vector<range_type> bounds = getBounds(kmers, threads, k_comp);
  assert(bounds.size() == threads);

  size_type fails = 0;
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

      std::string kmer = Key::decode(kmers[i].key, kmer_length, alpha);
      size_type endmarker_pos = kmer.find('$'); // The actual kmer ends at the first endmarker.
      if(endmarker_pos != std::string::npos) { kmer = kmer.substr(0, endmarker_pos + 1); }

      range_type range = this->find(kmer);
      if(Range::empty(range))
      {
        #pragma omp critical
        {
          std::cerr << "GCSA::verifyIndex(): find(" << kmer << ") returned empty range" << std::endl;
          fails++;
        }
        i = next; continue;
      }

      std::vector<node_type> expected;
      for(size_type j = i; j < next; j++) { expected.push_back(kmers[j].from); }
      removeDuplicates(expected, false);
      std::vector<node_type> occs;
      this->locate(range, occs);

      if(occs.size() != expected.size())
      {
        #pragma omp critical
        {
          std::cerr << "GCSA::verifyIndex(): Expected " << expected.size()
                    << " occurrences, got " << occs.size() << std::endl;
          printFailure(kmer, expected, occs); fails++;
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
              std::cerr << "GCSA::verifyIndex(): Failure at " << j << ": "
                        << "expected " << Node::decode(expected[j])
                        << ", got " << Node::decode(occs[j]) << std::endl;
              printFailure(kmer, expected, occs); fails++;
            }
            break;
          }
        }
      }

      i = next;
    }
  }

  if(fails == 0)
  {
    std::cout << "Index verification complete." << std::endl;
  }
  else
  {
    std::cout << "Index verification failed for " << fails << " patterns." << std::endl;
  }
  std::cout << std::endl;
  
  return fails == 0;
}

//------------------------------------------------------------------------------

struct FromGetter
{
  inline static size_type get(const PathNode& path)
  {
    return path.from;
  }
};

struct FirstGetter
{
  inline static size_type get(range_type range)
  {
    return range.first;
  }
};

template<class ValueType, class Getter>
struct ValueIndex
{
  sdsl::sd_vector<>               values;     // Marks the values that are present.
  sdsl::sd_vector<>::rank_1_type  value_rank;

  sdsl::bit_vector                first_occ;  // Marks the first occurrence of each rank.
  sdsl::bit_vector::select_1_type first_select;

  ValueIndex(const stxxl::vector<ValueType>& input)
  {
    stxxl::vector<size_type> buffer;
    this->first_occ = sdsl::bit_vector(input.size(), 0);

    size_type prev = ~(size_type)0;
    for(size_type i = 0; i < input.size(); i++)
    {
      size_type curr = Getter::get(input[i]);
      if(curr != prev)
      {
        buffer.push_back(curr);
        this->first_occ[i] = 1;
        prev = curr;
      }
    }

    // Fills in values, but only works if there are any values to fill
    if(buffer.size() > 0)
    {
      sdsl::sd_vector<> temp(buffer.begin(), buffer.end());
      this->values.swap(temp);
      sdsl::util::clear(buffer);
    }

    sdsl::util::init_support(this->value_rank, &(this->values));
    sdsl::util::init_support(this->first_select, &(this->first_occ));
  }

  // Finds the first occurrence of the value.
  size_type find(size_type value) const
  {
    if(value >= this->values.size() || this->values[value] == 0) { return this->first_occ.size(); }
    return this->first_select(this->value_rank(value) + 1);
  }
};

//------------------------------------------------------------------------------

/*
  Helper functions for splitting the work between threads.
*/

std::vector<size_type>
getTails(const std::vector<range_type>& bounds)
{
  std::vector<size_type> tail(bounds.size());
  for(size_type i = 0; i < bounds.size(); i++) { tail[i] = bounds[i].first; }
  return tail;
}

void
removeGaps(stxxl::vector<PathNode>& paths, const std::vector<range_type>& bounds, std::vector<size_type>& tail)
{
  for(size_type thread = 1; thread < bounds.size(); thread++)
  {
    for(size_type i = bounds[thread].first; i < tail[thread]; i++)
    {
      paths[tail[0]] = paths[i]; tail[0]++;
    }
  }
  paths.resize(tail[0]);
}

size_type
sumOf(const std::vector<size_type>& statistics)
{
  size_type res = 0;
  for(size_type i = 0; i < statistics.size(); i++) { res += statistics[i]; }
  return res;
}

//------------------------------------------------------------------------------
/*
void
writePaths(stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels,
  std::ostream& out, size_type size_limit,
  size_type& new_path_count, size_type& new_rank_count)
{
  size_type bytes_required = paths.size() * sizeof(PathNode) + labels.size() * sizeof(PathNode::rank_type);
  #pragma omp critical
  {
    bytes_required += (size_type)(out.tellp());
    if(bytes_required > size_limit * GIGABYTE)
    {
      std::cerr << "joinPaths(): Size limit exceeded, construction aborted" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for(auto& path : paths) { path.serialize(out, labels); }
    new_path_count += paths.size(); new_rank_count += labels.size();
  }
  paths.clear(); labels.clear();
}
*/

/*
  Join paths by left.to == right.from. Pre-/postcondition: paths are sorted by labels.
*/
void
joinPaths(stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels, size_type size_limit)
{
  // Initialization.
  PathFromComparator from_c; // Sort the paths by from.
  parallelQuickSort(paths.begin(), paths.end(), from_c);
  ValueIndex<PathNode, FromGetter> from_index(paths);
  size_type threads = omp_get_max_threads();
#ifdef VERBOSE_STATUS_INFO
  size_type old_path_count = paths.size();
#endif

  // Create a temporary file.
  /*
  std::string temp_file = tempFile(GCSA::EXTENSION);
  std::ofstream out(temp_file.c_str(), std::ios_base::binary);
  if(!out)
  {
    std::cerr << "joinPaths(): Cannot open temporary file " << temp_file
              << ", construction aborted" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  */
  //sdsl::write_member(new_path_count, out); sdsl::write_member(new_rank_count, out);

  // Create the next generation.
  stxxl::vector<PathNode> temp_nodes[threads];
  stxxl::vector<PathNode::rank_type> temp_labels[threads];
  #pragma omp parallel for schedule(static)
  for(size_type i = 0; i < paths.size(); i++)
  {
    size_type thread = omp_get_thread_num();
    if(paths[i].sorted())
    {
        // broken for multithreaded case
      temp_nodes[thread].push_back(PathNode(paths[i], labels, temp_labels[thread]));
    }
    else
    {
      size_type first = from_index.find(paths[i].to);
      for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++)
      {
        temp_nodes[thread].push_back(PathNode(paths[i], paths[j], labels, temp_labels[thread]));
      }
    }
    /*
    if(temp_nodes[thread].size() >= GCSA::WRITE_BUFFER_SIZE)
    {
      writePaths(temp_nodes[thread], temp_labels[thread], out, size_limit, new_path_count, new_rank_count);
    }
    */
  }
  
  size_type new_path_count = 0, new_rank_count = 0;
  for(size_type thread = 0; thread < threads; thread++) {
    new_path_count += temp_nodes[thread].size();
    new_rank_count += temp_labels[thread].size();
  }

  // now drop the results into our next iteration
  paths.clear(); paths.reserve(new_path_count+4);
  labels.clear(); labels.reserve(new_rank_count+4);
  
  for(size_type thread = 0; thread < threads; thread++)
  {
    for (auto& path_node : temp_nodes[thread]) {
      paths.push_back(PathNode(path_node, temp_labels[thread], labels));
    }
    /*
    for (auto& label : temp_labels[thread]) {
      labels.push_back(label);
    }
    */
    //writePaths(temp_nodes[thread], temp_labels[thread], out, size_limit, new_path_count, new_rank_count);
  }

  // Write the path/rank counts.
  /*
  out.seekp(0);
  sdsl::write_member(new_path_count, out); sdsl::write_member(new_rank_count, out);
  out.close();
  */

  // Replace the current generation with the next generation.
  //readPathNodes(temp_file, paths, labels);
  //labels.clear();
  
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "  joinPaths(): " << old_path_count << " -> " << paths.size() << " paths ("
            << labels.size() << " ranks)" << std::endl;
#endif

  // Restore the sorted order.
  PathFirstComparator first_c(labels);
  parallelQuickSort(paths.begin(), paths.end(), first_c);
}

//------------------------------------------------------------------------------

/*
  Iterates through ranges of paths with the same label. Empty range means start from
  the beginning, while range.first > bounds.second means the end.
*/
range_type
nextRange(range_type range, const stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels,
  range_type bounds)
{
  if(Range::empty(range)) {
      range.first = bounds.first; range.second = bounds.first; }
  else { range.first = range.second = range.second + 1; }
  PathFirstComparator pfc(labels);
  while(range.second + 1 <= bounds.second && !pfc(paths[range.first], paths[range.second + 1]))
  {
    range.second++;
  }
  return range;
}

range_type
firstRange(const stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels, range_type bounds)
{
  return nextRange(range_type(1, 0), paths, labels, bounds);
}

inline bool
sameFrom(range_type range, const stxxl::vector<PathNode>& paths)
{
  for(size_type i = range.first + 1; i <= range.second; i++)
  {
    if(paths[i].from != paths[range.first].from) { return false; }
  }
  return true;
}

//------------------------------------------------------------------------------

/*
  Tells whether it is safe to split the path array (sorted by labels) between
  'left' and 'right' (which should be adjacent). "Safe" means that the labels
  are different (so they will not be included in the same range) and the 'from'
  nodes are different (so a range cannot be extended to include both of them).
*/
struct SafeSplitComparator
{
  const stxxl::vector<PathNode::rank_type>& labels;

  explicit SafeSplitComparator(const stxxl::vector<PathNode::rank_type>& _labels) : labels(_labels) { }

  inline bool operator() (const PathNode& left, const PathNode& right) const
  {
    return ((left.firstLabel(0, labels) != right.firstLabel(0, labels))
      && (left.from != right.from));
  }
};

/*
  Extends the range forward into a maximal range of paths starting from the same node
  and sharing a common prefix that no other path has. Assumes that the input range only
  contains paths starting from the same node. Returns the lcp of the range as a pair
  (a,b), where a is the lcp of the labels and b is the lcp of the first diverging kmers.

  Looking at the previous PathNode may seem hairy, because it may belong to another thread.
  However, the last PathNode in a range never gets overwritten, so this should be safe.
*/
range_type
extendRange(range_type& range, const stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels,
  const LCP& lcp, range_type bounds)
{
  range_type min_lcp(0, 0);
  if(range.first > 0)
  {
    min_lcp = lcp.increment(lcp.max_lcp(paths[range.first - 1], paths[range.first], labels));
  }
  range_type range_lcp(paths[range.first].order(), 0);

  /*
    Iterate over one range at a time, doing the following tests:
    1. Is the from node still the same? Stop if not.
    2. Is the LCP still high enough? Stop if not.
    3. Is the LCP between the end of the range and the next path lower than the range lcp? Extend if true.
  */
  for(range_type next_range = nextRange(range, paths, labels, bounds); next_range.first <= bounds.second;
    next_range = nextRange(next_range, paths, labels, bounds))
  {
    if(paths[next_range.first].from != paths[range.first].from || !sameFrom(next_range, paths)) { break; }
    range_type next_lcp = lcp.min_lcp(paths[range.first], paths[next_range.second], labels);
    if(next_lcp < min_lcp) { break; }
    if(range.second + 1 <= bounds.second)
    {
      range_type border_lcp = lcp.max_lcp(paths[range.second], paths[range.second + 1], labels);
      if(border_lcp >= next_lcp) { continue; }
    }
    range.second = next_range.second; range_lcp = next_lcp;
  }

  return range_lcp;
}

/*
  Merges the path nodes into paths[range.first].
*/
void
mergePathNodes(stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels,
  range_type range, range_type range_lcp, const LCP& lcp)
{
  if(Range::length(range) == 1)
  {
    paths[range.first].makeSorted();
    return;
  }

  size_type order = range_lcp.first;
  if(range_lcp.second > 0)
  {
    LCP::rank_range
      ranks(paths[range.first].firstLabel(order, labels), paths[range.second].lastLabel(order, labels));
    ranks = lcp.extendRange(ranks, range_lcp.second);
    size_type ptr = paths[range.first].pointer();
    labels[ptr + order] = ranks.first;
    labels[ptr + order + 1] = ranks.second;
    order++;
  }

  paths[range.first].makeSorted();
  paths[range.first].setOrder(order); paths[range.first].setLCP(range_lcp.first);
  for(size_type i = range.first + 1; i <= range.second; i++)
  {
    paths[range.first].addPredecessors(paths[i]);
  }
}

/*
  Merges paths with adjacent labels and the same from node. Marks paths with unique
  labels sorted. Returns the number of unsorted path nodes.

  Note: The graph is not reverse deterministic. If a suffix of the path is unique, the
  entire path can be marked sorted, even though its label is non-unique.
*/
size_type
mergePaths(stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels, const LCP& lcp)
{
#ifdef VERBOSE_STATUS_INFO
  size_type old_path_count = paths.size();
#endif

  // Split the work between the threads so that the first rank is different at thread boundaries.
  size_type threads = omp_get_max_threads();
  SafeSplitComparator split_c(labels);
  std::vector<range_type> bounds = getBounds(paths, threads, split_c);
  std::vector<size_type> tail = getTails(bounds);
  std::vector<size_type> unique(threads, 0), unsorted(threads, 0), nondeterministic(threads, 0);

  #pragma omp parallel for schedule(static)
  for(size_type thread = 0; thread < threads; thread++)
  {
    for(range_type range = firstRange(paths, labels, bounds[thread]); range.first <= bounds[thread].second;
      range = nextRange(range, paths, labels, bounds[thread]))
    {
      if(sameFrom(range, paths))
      {
        range_type range_lcp = extendRange(range, paths, labels, lcp, bounds[thread]);
        mergePathNodes(paths, labels, range, range_lcp, lcp);
        paths[tail[thread]] = paths[range.first];
        tail[thread]++; unique[thread]++;
      }
      else
      {
        for(size_type i = range.first; i <= range.second; i++)
        {
          if(paths[i].sorted()) { nondeterministic[thread]++; }
          else { unsorted[thread]++; }
          paths[tail[thread]] = paths[i];
          tail[thread]++;
        }
      }
    }
  }

  // Remove the gaps and combine the statistics.
  removeGaps(paths, bounds, tail);
  unique[0] = sumOf(unique);
  unsorted[0] = sumOf(unsorted);
  nondeterministic[0] = sumOf(nondeterministic);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "  mergePaths(): " << old_path_count << " -> " << paths.size() << " paths" << std::endl;
  std::cerr << "  mergePaths(): " << unique[0] << " unique, " << unsorted[0] << " unsorted, "
            << nondeterministic[0] << " nondeterministic paths" << std::endl;
#endif
  return unsorted[0];
}

//------------------------------------------------------------------------------

void
GCSA::prefixDoubling(stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels,
  size_type kmer_length, size_type doubling_steps, size_type size_limit,
  const LCP& lcp)
{
  size_type path_order = 1;
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::prefixDoubling(): Initial path length " << kmer_length << std::endl;
#endif
  size_type unsorted = mergePaths(paths, labels, lcp);

  for(size_type step = 1; step <= doubling_steps && unsorted > 0; step++)
  {
#ifdef VERBOSE_STATUS_INFO
    std::cerr << "GCSA::prefixDoubling(): Step " << step << " (path length " << (path_order * kmer_length) << " -> "
              << (2 * path_order * kmer_length) << ")" << std::endl;
#endif
    joinPaths(paths, labels, size_limit); path_order *= 2;
    unsorted = mergePaths(paths, labels, lcp);
  }
  this->max_query_length = (unsorted == 0 ? ~(size_type)0 : kmer_length << doubling_steps);
}

//------------------------------------------------------------------------------

/*
  This version assumes that the paths have identical labels.
*/
void
mergePathNodes(stxxl::vector<PathNode>& paths, range_type range)
{
  paths[range.first].makeSorted();
  for(size_type i = range.first + 1; i <= range.second; i++)
  {
    paths[range.first].addPredecessors(paths[i]);
  }
}

void
GCSA::mergeByLabel(stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels,
  stxxl::vector<range_type>& from_nodes)
{
  sdsl::util::clear(from_nodes);
#ifdef VERBOSE_STATUS_INFO
  size_type old_path_count = paths.size();
#endif

  size_type threads = omp_get_max_threads();
  PathFirstComparator first_c(labels);
  std::vector<range_type> bounds = getBounds(paths, threads, first_c);
  std::vector<size_type> tail = getTails(bounds);
  std::vector<range_type> from_buffer[threads];

  #pragma omp parallel for schedule(static)
  for(size_type thread = 0; thread < threads; thread++)
  {
    for(range_type range = firstRange(paths, labels, bounds[thread]); range.first <= bounds[thread].second;
      range = nextRange(range, paths, labels, bounds[thread]))
    {
      mergePathNodes(paths, range);
      paths[tail[thread]] = paths[range.first];
      for(size_type i = range.first + 1; i <= range.second; i++)
      {
        from_buffer[thread].push_back(range_type(tail[thread], paths[i].from));
      }
      tail[thread]++;
    }
  }

  // Combine and correct the additional samples.
  size_type total_from_nodes = 0;
  for(size_type thread = 0; thread < threads; thread++) { total_from_nodes += from_buffer[thread].size(); }
  from_nodes.reserve(total_from_nodes);
  for(size_type thread = 0, gap = 0; thread < threads; thread++)
  {
    for(size_type i = 0; i < from_buffer[thread].size(); i++)
    {
      range_type temp = from_buffer[thread][i];
      temp.first -= gap;
      from_nodes.push_back(temp);
    }
    sdsl::util::clear(from_buffer[thread]);
    if(thread + 1 < threads) { gap += bounds[thread + 1].first - tail[thread]; }
  }

  removeGaps(paths, bounds, tail);
  this->path_node_count = paths.size();

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::mergeByLabel(): " << old_path_count << " -> " << paths.size() << " paths" << std::endl;
#endif
}

//------------------------------------------------------------------------------

std::pair<PathLabel, PathLabel>
predecessor(const PathNode& curr, stxxl::vector<PathNode::rank_type>& labels,
  comp_type comp, const DeBruijnGraph& mapper, const sdsl::int_vector<0>& last_char)
{
  size_type i = 0, j = curr.pointer();
  PathLabel first, last;
  first.is_first = true; last.is_first = false;

  // Handle the common prefix of the labels.
  while(i < curr.lcp())
  {
    first.label[i] = last.label[i] = mapper.LF(labels[j], comp);
    comp = last_char[labels[j]];
    i++; j++;
  }

  // Handle the diverging suffixes of the labels.
  comp_type first_comp = comp, last_comp = comp;
  if(i < curr.order())
  {
    first.label[i] = mapper.LF(labels[j], first_comp);
    first_comp = last_char[labels[j]];
    last.label[i] = mapper.LF(labels[j + 1], last_comp);
    last_comp = last_char[labels[j + 1]];
    i++;
  }
  if(i < PathLabel::LABEL_LENGTH)
  {
    first.label[i] = mapper.node_rank(mapper.alpha.C[first_comp]);
    last.label[i] = mapper.node_rank(mapper.alpha.C[last_comp + 1]) - 1;
    i++;
  }

  first.length = last.length = i;
  return std::make_pair(first, last);
}

/*
  Path nodes correspond to lexicographic ranges of path labels. The predecessor found using
  mapper.LF() matches all path nodes with ranges intersecting the predecessor range.

  FIXME Later: parallelize
*/
void
GCSA::build(stxxl::vector<PathNode>& paths, stxxl::vector<PathNode::rank_type>& labels,
  DeBruijnGraph& mapper, sdsl::int_vector<0>& last_char)
{
  for(size_type i = 0; i < paths.size(); i++) { paths[i].initDegree(); }

  // Pointers to the next path nodes with labels starting with the given comp value.
  size_type sigma = mapper.alpha.sigma;
  size_type next[sigma + 1];
  for(size_type comp = 0; comp < sigma; comp++) { next[comp] = mapper.charRange(comp).first; }
  next[sigma] = ~(size_type)0;
  for(size_type i = 0, comp = 0; i < paths.size(); i++)
  {
    while(paths[i].firstLabel(0, labels) >= next[comp]) { next[comp] = i; comp++; }
  }

  size_type total_edges = 0;
  sdsl::int_vector<64> counts(sigma, 0);
  sdsl::int_vector<8> bwt_buffer(paths.size() + paths.size() / 2, 0);
  this->path_nodes = bit_vector(bwt_buffer.size(), 0);
  for(size_type i = 0; i < paths.size(); i++)
  {
    for(size_type comp = 0; comp < sigma; comp++)
    {
      if(!(paths[i].hasPredecessor(comp))) { continue; }

      // Find the predecessor of paths[i] with comp and the first path intersecting it.
      std::pair<PathLabel, PathLabel> pred = predecessor(paths[i], labels, comp, mapper, last_char);
      while(!(paths[next[comp]].intersect(pred.first, pred.second, labels))) {
          next[comp]++;
      }
      paths[i].incrementIndegree(); paths[next[comp]].incrementOutdegree(); total_edges++;
      if(total_edges > bwt_buffer.size()) { bwt_buffer.resize(bwt_buffer.size() + paths.size() / 2); }
      bwt_buffer[total_edges - 1] = comp; counts[comp]++;

      // Create additional edges if the next path nodes also intersect with the predecessor.
      while(next[comp] + 1 < paths.size() && paths[next[comp] + 1].intersect(pred.first, pred.second, labels))
      {
        next[comp]++;
        paths[i].incrementIndegree(); paths[next[comp]].incrementOutdegree(); total_edges++;
        if(total_edges > bwt_buffer.size()) { bwt_buffer.resize(bwt_buffer.size() + paths.size() / 2); }
        bwt_buffer[total_edges - 1] = comp; counts[comp]++;
      }
    }

    if(total_edges > this->path_nodes.size()) { this->path_nodes.resize(bwt_buffer.size()); }
    this->path_nodes[total_edges - 1] = 1;
  }

  // Init alpha and bwt; clear unnecessary structures.
  this->alpha = Alphabet(counts, mapper.alpha.char2comp, mapper.alpha.comp2char);
  sdsl::util::clear(mapper); sdsl::util::clear(last_char);
  bwt_buffer.resize(total_edges); this->path_nodes.resize(total_edges);
  directConstruct(this->bwt, bwt_buffer); sdsl::util::clear(bwt_buffer);

  // Init edges and rank/select support.
  this->edges = bit_vector(total_edges, 0); total_edges = 0;
  for(size_type i = 0; i < paths.size(); i++)
  {
    total_edges += paths[i].outdegree();
    this->edges[total_edges - 1] = 1;
  }
  this->initSupport();

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::build(): " << total_edges << " edges" << std::endl;
#endif

  // This is a useful test when something goes wrong.
/*  for(size_type i = 0; i < paths.size(); i++)
  {
    if(paths[i].outdegree() == 0)
    {
      std::cout << "Path " << i << ": ";
      paths[i].print(std::cout, labels);
      std::cout << " has outdegree 0" << std::endl;
    }
  }*/
}

//------------------------------------------------------------------------------

void
GCSA::initSupport()
{
  sdsl::util::init_support(this->path_rank, &(this->path_nodes));
  sdsl::util::init_support(this->path_select, &(this->path_nodes));
  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->edge_select, &(this->edges));
}

//------------------------------------------------------------------------------

std::vector<node_type>
fromNodes(size_type path, const stxxl::vector<PathNode>& paths,
  size_type& additional, const stxxl::vector<range_type>& from_nodes)
{
  std::vector<node_type> res;
  res.push_back(paths[path].from);

  while(additional < from_nodes.size() && from_nodes[additional].first < path) { additional++; }
  while(additional < from_nodes.size() && from_nodes[additional].first == path)
  {
    res.push_back(from_nodes[additional].second); additional++;
  }

  removeDuplicates(res, false);
  return res;
}

void
GCSA::sample(stxxl::vector<PathNode>& paths, stxxl::vector<range_type>& from_nodes)
{
  this->sampled_paths = bit_vector(paths.size(), 0);
  this->samples = bit_vector(paths.size() + from_nodes.size(), 0);

  size_type sample_bits = 0;
  std::vector<node_type> sample_buffer;
  ValueIndex<range_type, FirstGetter> from_index(from_nodes);
  for(size_type i = 0, j = 0; i < paths.size(); i++)
  {
    bool sample_this = false;
    std::vector<node_type> curr = fromNodes(i, paths, j, from_nodes);
    if(paths[i].indegree() > 1) { sample_this = true; }
    if(paths[i].hasPredecessor(ENDMARKER_COMP)) { sample_this = true; }
    for(size_type k = 0; k < curr.size(); k++)
    {
      if(Node::offset(curr[k]) == 0) { sample_this = true; break; }
    }

    if(!sample_this)  // Compare to the from nodes at the predecessor.
    {
      size_type pred = this->LF(i);
      size_type temp = from_index.find(pred);
      std::vector<node_type> prev = fromNodes(pred, paths, temp, from_nodes);
      if(prev.size() != curr.size()) { sample_this = true; }
      else
      {
        for(size_type k = 0; k < curr.size(); k++)
        {
          if(curr[k] != prev[k] + 1) { sample_this = true; break; }
        }
      }
    }

    if(sample_this)
    {
      this->sampled_paths[i] = 1;
      for(size_type k = 0; k < curr.size(); k++)
      {
        sample_bits = std::max(sample_bits, bit_length(curr[k]));
        sample_buffer.push_back(curr[k]);
      }
      this->samples[sample_buffer.size() - 1] = 1;
    }
  }
  sdsl::util::clear(from_nodes);

  sdsl::util::init_support(this->sampled_path_rank, &(this->sampled_paths));
  this->samples.resize(sample_buffer.size());
  sdsl::util::init_support(this->sample_select, &(this->samples));

  this->stored_samples = sdsl::int_vector<0>(sample_buffer.size(), 0, sample_bits);
  for(size_type i = 0; i < sample_buffer.size(); i++) { this->stored_samples[i] = sample_buffer[i]; }
  sdsl::util::clear(sample_buffer);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::sample(): " << this->sample_count() << " samples at "
            << this->sampled_positions() << " positions" << std::endl;
#endif
}

//------------------------------------------------------------------------------

void
GCSA::locate(size_type path_node, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { sdsl::util::clear(results); }
  if(path_node >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  this->locateInternal(path_node, results);
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locate(range_type range, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { sdsl::util::clear(results); }
  if(Range::empty(range) || range.second >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  for(size_type i = range.first; i <= range.second; i++)
  {
    this->locateInternal(i, results);
  }
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locateInternal(size_type path_node, std::vector<node_type>& results) const
{
  size_type steps = 0;
  while(this->sampled_paths[path_node] == 0)
  {
    path_node = this->LF(path_node);
    steps++;
  }

  range_type sample_range = this->sampleRange(path_node);
  for(size_type i = sample_range.first; i <= sample_range.second; i++)
  {
    results.push_back(this->stored_samples[i] + steps);
  }
}

//------------------------------------------------------------------------------

} // namespace gcsa
