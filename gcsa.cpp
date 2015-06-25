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

GCSA::GCSA(const std::vector<key_type>& keys, size_type kmer_length, const Alphabet& _alpha)
{
  this->path_node_count = keys.size();
  this->max_query_length = kmer_length;

  size_type total_edges = 0;
  for(size_type i = 0; i < keys.size(); i++) { total_edges += sdsl::bits::lt_cnt[Key::predecessors(keys[i])]; }

  sdsl::int_vector<64> counts(_alpha.sigma, 0);
  sdsl::int_vector<8> buffer(total_edges, 0);
  this->path_nodes = bit_vector(total_edges, 0);
  this->edges = bit_vector(total_edges, 0);
  for(size_type i = 0, bwt_pos = 0, edge_pos = 0; i < keys.size(); i++)
  {
    size_type pred = Key::predecessors(keys[i]);
    for(size_type j = 0; j < _alpha.sigma; j++)
    {
      if(pred & (((size_type)1) << j))
      {
        buffer[bwt_pos] = j; bwt_pos++;
        counts[j]++;
      }
    }
    this->path_nodes[bwt_pos - 1] = 1;
    edge_pos += sdsl::bits::lt_cnt[Key::successors(keys[i])];
    this->edges[edge_pos - 1] = 1;
  }
  directConstruct(this->bwt, buffer);
  this->alpha = Alphabet(counts, _alpha.char2comp, _alpha.comp2char);

  this->initSupport();
}

GCSA::GCSA(std::vector<KMer>& kmers, size_type kmer_length, size_type doubling_steps, const Alphabet& _alpha)
{
  if(kmers.size() == 0) { return; }
  if(doubling_steps > DOUBLING_STEPS)
  {
    std::cerr << "GCSA::GCSA(): The number of doubling steps is too high: " << doubling_steps << std::endl;
    std::cerr << "GCSA::GCSA(): Reverting the number of doubling steps to " << DOUBLING_STEPS << std::endl;
    doubling_steps = DOUBLING_STEPS;
  }

  // Sort the kmers, build the mapper GCSA for generating the edges.
  std::vector<key_type> keys;
  sdsl::int_vector<0> last_char;
  uniqueKeys(kmers, keys, last_char);
  GCSA mapper(keys, kmer_length, _alpha);
  sdsl::util::clear(keys);

  // Transform the kmers into PathNodes.
  // FIXME Later: Save memory by not having both KMers and PathNodes in memory.
  std::vector<PathNode> paths; paths.reserve(kmers.size());
  for(size_type i = 0; i < kmers.size(); i++) { paths.push_back(PathNode(kmers[i])); }
  sdsl::util::clear(kmers);

  // Build the GCSA in PathNodes.
  this->prefixDoubling(paths, kmer_length, doubling_steps);
  std::vector<range_type> from_nodes;
  this->mergeByLabel(paths, from_nodes);
  this->build(paths, mapper, last_char);

  this->sample(paths, from_nodes);
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

  ValueIndex(const std::vector<ValueType>& input)
  {
    std::vector<size_type> buffer;
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

    sdsl::sd_vector<> temp(buffer.begin(), buffer.end());
    this->values.swap(temp);
    sdsl::util::clear(buffer);

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
  Join paths by left.to == right.from.
  FIXME Later: Write the next generation to disk.
  FIXME Later: parallelize
*/
void
joinPaths(std::vector<PathNode>& paths)
{
  PathFromComparator pfc; // Sort the paths by from.
  parallelQuickSort(paths.begin(), paths.end(), pfc);
  size_type old_path_count = paths.size();

  ValueIndex<PathNode, FromGetter> from_index(paths);
  size_type new_path_count = 0;
  for(size_type i = 0; i < paths.size(); i++)
  {
    if(paths[i].sorted()) { new_path_count++; }
    else
    {
      size_type first = from_index.find(paths[i].to);
      for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++) { new_path_count++; }
    }
  }
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "  joinPaths(): Reserving space for " << new_path_count << " paths" << std::endl;
#endif

  std::vector<PathNode> next; next.reserve(new_path_count);
  for(size_type i = 0; i < paths.size(); i++)
  {
    if(paths[i].sorted()) { next.push_back(paths[i]); }
    else
    {
      size_type first = from_index.find(paths[i].to);
      for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++)
      {
        next.push_back(PathNode(paths[i], paths[j]));
      }
    }
  }

  paths.swap(next);
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "  joinPaths(): " << old_path_count << " -> " << paths.size() << " paths" << std::endl;
#endif
}

//------------------------------------------------------------------------------

/*
  Iterates through ranges of paths with the same label. Empty range means start from
  the beginning, while range.first >= paths.size() means the end.
*/
range_type
nextRange(range_type range, const std::vector<PathNode>& paths)
{
  if(Range::empty(range)) { range.first = 0; range.second = 0; }
  else { range.first = range.second = range.second + 1; }

  while(range.second + 1 < paths.size() && !(paths[range.first] < paths[range.second + 1]))
  {
    range.second++;
  }

  return range;
}

range_type
firstRange(const std::vector<PathNode>& paths)
{
  return nextRange(range_type(1, 0), paths);
}

bool
sameFrom(range_type range, const std::vector<PathNode>& paths)
{
  for(size_type i = range.first + 1; i <= range.second; i++)
  {
    if(paths[i].from != paths[range.first].from) { return false; }
  }
  return true;
}

/*
  Extends the range forward into a maximal range of paths starting from the same node
  and sharing a common prefix that no other path has. Assumes that the input range only
  contains paths starting from the same node. Returns the lcp of the range.
*/
size_type
extendRange(range_type& range, const std::vector<PathNode>& paths)
{
  size_type min_lcp = (range.first == 0 ? 0 : paths[range.first].lcp(paths[range.first - 1]) + 1);
  size_type range_lcp = paths[range.first].order();
  size_type candidate_lcp = range_lcp;

  for(size_type i = range.second + 1; i <= paths.size(); i++)
  {
    if(i >= paths.size()) { range.second = paths.size() - 1; range_lcp = candidate_lcp; break; }
    if(paths[i].from != paths[range.first].from) { break; }
    size_type lcp = paths[i].lcp(paths[range.first]);
    if(lcp < candidate_lcp) { range.second = i - 1; range_lcp = candidate_lcp; candidate_lcp = lcp; }
    if(lcp < min_lcp) { break; }
  }

  return range_lcp;
}

PathNode
mergePathNodes(const std::vector<PathNode>& paths, range_type range, size_type lcp)
{
  PathNode res = paths[range.first];

  res.makeSorted(); res.setOrder(lcp); res.pad();
  for(size_type i = range.first + 1; i <= range.second; i++) { res.addPredecessors(paths[i]); }

  return res;
}

/*
  Merges paths with adjacent labels and the same from node. Marks paths with unique
  labels sorted. Returns the number of unsorted path nodes.

  Note: The graph is not reverse deterministic. If a suffix of the path is unique, the
  entire path can be marked sorted, even though its label is non-unique.
*/
size_type
mergePaths(std::vector<PathNode>& paths)
{
  parallelQuickSort(paths.begin(), paths.end());
  size_type old_path_count = paths.size();

  std::vector<PathNode> new_last;
  size_type tail = 0, unique = 0, unsorted = 0, nondeterministic = 0;
  for(range_type range = firstRange(paths); range.first < paths.size(); range = nextRange(range, paths))
  {
    if(sameFrom(range, paths))
    {
      size_type lcp = extendRange(range, paths);
      paths[tail] = mergePathNodes(paths, range, lcp);
      tail++; unique++;
    }
    else
    {
      for(size_type i = range.first; i <= range.second; i++)
      {
        if(paths[i].sorted()) { nondeterministic++; }
        else { unsorted++; }
        paths[tail] = paths[i];
        tail++;
      }
    }
  }

  paths.resize(tail);
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "  mergePaths(): " << old_path_count << " -> " << paths.size() << " paths" << std::endl;
  std::cerr << "  mergePaths(): " << unique << " unique, " << unsorted << " unsorted, "
            << nondeterministic << " nondeterministic paths" << std::endl;
#endif
  return unsorted;
}

//------------------------------------------------------------------------------

void
GCSA::prefixDoubling(std::vector<PathNode>& paths, size_type kmer_length, size_type doubling_steps)
{
  size_type path_order = 1;
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::prefixDoubling(): Initial path length " << kmer_length << std::endl;
#endif
  size_type unsorted = mergePaths(paths);

  for(size_type step = 1; step <= doubling_steps && unsorted > 0; step++)
  {
#ifdef VERBOSE_STATUS_INFO
    std::cerr << "GCSA::prefixDoubling(): Step " << step << " (path length " << (path_order * kmer_length) << " -> "
              << (2 * path_order * kmer_length) << ")" << std::endl;
#endif
    joinPaths(paths); path_order *= 2;
    unsorted = mergePaths(paths);
  }
  this->max_query_length = (unsorted == 0 ? ~(size_type)0 : kmer_length << doubling_steps);
}

//------------------------------------------------------------------------------

void
GCSA::mergeByLabel(std::vector<PathNode>& paths, std::vector<range_type>& from_nodes)
{
  this->path_node_count = 0;
  size_type old_path_count = paths.size();

  for(range_type range = firstRange(paths); range.first < paths.size(); range = nextRange(range, paths))
  {
    paths[this->path_node_count] = mergePathNodes(paths, range, paths[range.first].order());
    for(size_type i = range.first + 1; i <= range.second; i++)
    {
      from_nodes.push_back(range_type(this->path_node_count, paths[i].from));
    }
    this->path_node_count++;
  }

  paths.resize(this->path_node_count);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::mergeByLabel(): " << old_path_count << " -> " << paths.size() << " paths" << std::endl;
#endif
}

//------------------------------------------------------------------------------

PathNode
predecessor(const PathNode& curr, comp_type comp, const GCSA& mapper, const sdsl::int_vector<0>& last_char)
{
  PathNode pred;
  size_type order = curr.order();

  for(size_type i = 0; i < order; i++)
  {
    pred.first_label[i] = pred.last_label[i] = mapper.LF(curr.first_label[i], comp);
    comp = last_char[curr.first_label[i]];
  }
  if(order < PathNode::LABEL_LENGTH)
  {
    pred.first_label[order] = mapper.edge_rank(mapper.alpha.C[comp]);
    pred.last_label[order] = mapper.edge_rank(mapper.alpha.C[comp + 1]) - 1;
    order++;
  }
  pred.setOrder(order); pred.pad();

  return pred;
}

/*
  Path nodes correspond to lexicographic ranges of path labels. The predecessor found using
  mapper.LF() matches all path nodes with ranges intersecting the predecessor range.

  FIXME Later: parallelize
*/
void
GCSA::build(std::vector<PathNode>& paths, GCSA& mapper, sdsl::int_vector<0>& last_char)
{
  for(size_type i = 0; i < paths.size(); i++) { paths[i].initDegree(); }

  // Pointers to the next path nodes with labels starting with the given comp value.
  size_type sigma = mapper.alpha.sigma;
  size_type next[sigma + 1];
  for(size_type comp = 0; comp < sigma; comp++) { next[comp] = mapper.charRange(comp).first; }
  next[sigma] = ~(size_type)0;
  for(size_type i = 0, comp = 0; i < paths.size(); i++)
  {
    while(paths[i].first_label[0] >= next[comp]) { next[comp] = i; comp++; }
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
      PathNode pred = predecessor(paths[i], comp, mapper, last_char);
      while(!(pred.intersect(paths[next[comp]]))) { next[comp]++; }
      paths[i].incrementIndegree(); paths[next[comp]].incrementOutdegree(); total_edges++;
      if(total_edges > bwt_buffer.size()) { bwt_buffer.resize(bwt_buffer.size() + paths.size() / 2); }
      bwt_buffer[total_edges - 1] = comp; counts[comp]++;

      // Create additional edges if the next path nodes also intersect with the predecessor.
      while(next[comp] + 1 < paths.size() && pred.intersect(paths[next[comp] + 1]))
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

  // This is a useful test when something goes wrong.
/*  for(size_type i = 0; i < paths.size(); i++)
  {
    if(paths[i].outdegree() == 0)
    {
      std::cout << "Path " << i << ": " << paths[i] << " has outdegree 0" << std::endl;
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
fromNodes(size_type path, const std::vector<PathNode>& paths,
  size_type& additional, const std::vector<range_type>& from_nodes)
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
GCSA::sample(std::vector<PathNode>& paths, std::vector<range_type>& from_nodes)
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
          if(Node::id(curr[k]) != Node::id(prev[k]) || Node::offset(curr[k]) != Node::offset(prev[k]) + 1)
          {
            sample_this = true; break;
          }
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
}

//------------------------------------------------------------------------------

range_type
GCSA::find(const char_type* pattern, size_type length) const
{
  return this->find(pattern, pattern + length);
}

range_type
GCSA::find(const char* pattern, size_type length) const
{
  return this->find(pattern, pattern + length);
}

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
