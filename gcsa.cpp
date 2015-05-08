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
  this->node_count = 0;
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
  this->node_count = g.node_count;
  this->max_query_length = g.max_query_length;

  this->bwt = g.bwt;
  this->alpha = g.alpha;

  this->nodes = g.nodes;
  this->node_rank = g.node_rank;
  this->node_select = g.node_select;

  this->edges = g.edges;
  this->edge_rank = g.edge_rank;
  this->edge_select = g.edge_select;

  this->sampled_nodes = g.sampled_nodes;
  this->sampled_node_rank = g.sampled_node_rank;

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
    std::swap(this->node_count, g.node_count);
    std::swap(this->max_query_length, g.max_query_length);

    this->bwt.swap(g.bwt);
    this->alpha.swap(g.alpha);

    this->nodes.swap(g.nodes);
    sdsl::util::swap_support(this->node_rank, g.node_rank, &(this->nodes), &(g.nodes));
    sdsl::util::swap_support(this->node_select, g.node_select, &(this->nodes), &(g.nodes));

    this->edges.swap(g.edges);
    sdsl::util::swap_support(this->edge_rank, g.edge_rank, &(this->edges), &(g.edges));
    sdsl::util::swap_support(this->edge_select, g.edge_select, &(this->edges), &(g.edges));

    this->sampled_nodes.swap(g.sampled_nodes);
    sdsl::util::swap_support(this->sampled_node_rank, g.sampled_node_rank, &(this->sampled_nodes), &(g.sampled_nodes));

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
    this->node_count = std::move(g.node_count);
    this->max_query_length = std::move(g.max_query_length);

    this->bwt = std::move(g.bwt);
    this->alpha = std::move(g.alpha);

    this->nodes = std::move(g.nodes);
    this->node_rank = std::move(g.node_rank);
    this->node_select = std::move(g.node_select);

    this->edges = std::move(g.edges);
    this->edge_rank = std::move(g.edge_rank);
    this->edge_select = std::move(g.edge_select);

    this->sampled_nodes = std::move(g.sampled_nodes);
    this->sampled_node_rank = std::move(g.sampled_node_rank);

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
  this->node_rank.set_vector(&(this->nodes));
  this->node_select.set_vector(&(this->nodes));

  this->edge_rank.set_vector(&(this->edges));
  this->edge_select.set_vector(&(this->edges));

  this->sampled_node_rank.set_vector(&(this->sampled_nodes));

  this->sample_select.set_vector(&(this->samples));
}

GCSA::size_type
GCSA::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->node_count, out, child, "node_count");
  written_bytes += sdsl::write_member(this->max_query_length, out, child, "max_query_length");

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->alpha.serialize(out, child, "alpha");

  written_bytes += this->nodes.serialize(out, child, "nodes");
  written_bytes += this->node_rank.serialize(out, child, "node_rank");
  written_bytes += this->node_select.serialize(out, child, "node_select");

  written_bytes += this->edges.serialize(out, child, "edges");
  written_bytes += this->edge_rank.serialize(out, child, "edge_rank");
  written_bytes += this->edge_select.serialize(out, child, "edge_select");

  written_bytes += this->sampled_nodes.serialize(out, child, "sampled_nodes");
  written_bytes += this->sampled_node_rank.serialize(out, child, "sampled_node_rank");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSA::load(std::istream& in)
{
  sdsl::read_member(this->node_count, in);
  sdsl::read_member(this->max_query_length, in);

  this->bwt.load(in);
  this->alpha.load(in);

  this->nodes.load(in);
  this->node_rank.load(in, &(this->nodes));
  this->node_select.load(in, &(this->nodes));

  this->edges.load(in);
  this->edge_rank.load(in, &(this->edges));
  this->edge_select.load(in, &(this->edges));

  this->sampled_nodes.load(in);
  this->sampled_node_rank.load(in, &(this->sampled_nodes));

  this->stored_samples.load(in);
  this->samples.load(in);
  this->sample_select.load(in, &(this->samples));
}

//------------------------------------------------------------------------------

GCSA::GCSA(const std::vector<key_type>& keys, size_type kmer_length, const Alphabet& _alpha)
{
  this->node_count = keys.size();
  this->max_query_length = kmer_length;

  size_type total_edges = 0, alt_edges = 0;
  for(size_type i = 0; i < keys.size(); i++) { total_edges += sdsl::bits::lt_cnt[Key::predecessors(keys[i])];
  alt_edges += sdsl::bits::lt_cnt[Key::successors(keys[i])]; }

  sdsl::int_vector<64> counts(_alpha.sigma, 0);
  sdsl::int_vector<8> buffer(total_edges, 0);
  this->nodes = bit_vector(total_edges, 0);
  this->edges = bit_vector(total_edges, 0);
  for(size_type i = 0, bwt_pos = 0, edge_pos = 0; i < keys.size(); i++)
  {
    byte_type pred = Key::predecessors(keys[i]);
    for(size_type j = 0; j < _alpha.sigma; j++)
    {
      if(pred & (((size_type)1) << j))
      {
        buffer[bwt_pos] = j; bwt_pos++;
        counts[j]++;
      }
    }
    this->nodes[bwt_pos - 1] = 1;
    edge_pos += sdsl::bits::lt_cnt[Key::successors(keys[i])];
    this->edges[edge_pos - 1] = 1;
  }
  directConstruct(this->bwt, buffer);
  this->alpha = Alphabet(counts, _alpha.char2comp, _alpha.comp2char);

  sdsl::util::init_support(this->node_rank, &(this->nodes));
  sdsl::util::init_support(this->node_select, &(this->nodes));
  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->edge_select, &(this->edges));
}

GCSA::GCSA(std::vector<KMer>& kmers, size_type kmer_length, const Alphabet& _alpha)
{
  if(kmers.size() == 0) { return; }

  // Sort the kmers, build the mapper GCSA for generating the edges.
  std::vector<key_type> keys;
  sdsl::int_vector<0> last_chars;
  uniqueKeys(kmers, keys, last_chars);
  GCSA mapper(keys, kmer_length, _alpha);
  sdsl::util::clear(keys);

  // Transform the kmers into PathNodes.
  // FIXME Later: Save memory by not having both KMers and PathNodes in memory.
  std::vector<PathNode<>> paths(kmers.size());
  for(size_type i = 0; i < kmers.size(); i++) { paths[i] = PathNode<>(kmers[i]); }
  sdsl::util::clear(kmers);

  size_type path_order = this->prefixDoubling(paths, kmer_length);
  this->sample(paths, path_order);

  /*
    Edge generation:
    - set the to fields of the last generation to 0,
    - build an index structure to find paths quickly by labels
    - scan the paths in the last generation, writing the number of outgoing edges to the to field
      * use mapper, predecessor field in keys, and last_chars to find the predecessors
      * increment the to field of the predecessor

    FIXME Later: parallelize
  */

  /*
    GCSA construction:
    - count character occurrences
    - set bwt, alpha
    - set nodes, node_rank, node_select, edges, edge_rank, edge_select
  */

  std::cout << std::endl;
}

//------------------------------------------------------------------------------

struct FromIndex
{
  sdsl::bit_vector from_values; // Marks the from values that are present.
  sdsl::bit_vector::rank_1_type from_rank;

  sdsl::bit_vector first_occ;   // Marks the first occurrence of each rank.
  sdsl::bit_vector::select_1_type first_select;

  FromIndex(const std::vector<PathNode<>>& paths)
  {
    this->from_values = sdsl::bit_vector(paths[paths.size() - 1].from + 1, 0);
    this->first_occ = sdsl::bit_vector(paths.size(), 0);

    node_type prev = this->from_values.size();
    for(size_type i = 0; i < paths.size(); i++)
    {
      node_type curr = paths[i].from;
      if(curr != prev)
      {
        this->from_values[curr] = 1;
        this->first_occ[i] = 1;
        prev = curr;
      }
    }

    sdsl::util::init_support(this->from_rank, &(this->from_values));
    sdsl::util::init_support(this->first_select, &(this->first_occ));
  }

  // Finds the first occurrence of the node in the from fields.
  size_type find(node_type node) const
  {
    if(this->from_values[node] == 0) { return this->first_occ.size(); }
    return this->first_select(this->from_rank(node) + 1);
  }
};

/*
  Join paths by left.to == right.from.
  FIXME Later: Predict the size of the next generation or write it to disk.
  FIXME Later: parallelize
*/
void
joinPaths(std::vector<PathNode<>>& paths, size_type& path_order)
{
  FromIndex from_index(paths);

  std::vector<PathNode<>> next;
  for(size_type i = 0; i < paths.size(); i++)
  {
    if(paths[i].sorted()) { next.push_back(paths[i]); continue; }
    size_type first = from_index.find(paths[i].to);
    for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++)
    {
      next.push_back(PathNode<>(paths[i], paths[j], path_order));
    }
  }
  paths.swap(next);

  path_order *= 2;
}

/*
  Iterates through ranges of paths with the same label. Empty range means start from
  the beginning, while range.first >= paths.size() means the end. Returns true if
  the range contains only nodes with the same label. In that case the range may cover
  multiple adjacent labels.
*/
bool
nextRange(range_type& range, const std::vector<PathNode<>>& paths,
  const PathLabelComparator<>& plc)
{
  if(Range::empty(range)) { range.first = 0; range.second = 0; }
  else { range.first = range.second + 1; range.second = range.first; }
  if(range.first >= paths.size()) { return true; }

  bool same_from = true;
  size_type rank_start = range.first, next = range.first + 1;
  while(next < paths.size())
  {
    if(plc(paths[rank_start], paths[next])) // A range ends at next - 1.
    {
      range.second = next - 1; rank_start = next;
      if(!same_from) { return false; }  // Found a range with different from nodes.
    }
    if(paths[next].from != paths[range.first].from)
    {
      if(rank_start != range.first) { return same_from; } // Already found at least one range.
      same_from = false;
    }
    next++;
  }
  range.second = paths.size() - 1;

  return same_from;
}

/*
  Merges paths with adjacent labels and the same from node. Marks paths with unique
  labels sorted.
*/
size_type
mergePaths(std::vector<PathNode<>>& paths, size_type path_order)
{
  PathLabelComparator<> plc(path_order);
  parallelQuickSort(paths.begin(), paths.end(), plc); // Sort paths by labels.

  size_type tail = 0, unique = 0;
  range_type range(1, 0);
  while(true)
  {
    bool same_from = nextRange(range, paths, plc);
    if(range.first >= paths.size()) { break; }
    if(same_from)
    {
      paths[tail] = paths[range.first];
      paths[tail].makeSorted();
      tail++; unique++;
    }
    else
    {
      for(size_type i = range.first; i <= range.second; i++)
      {
        paths[tail] = paths[i]; tail++;
      }
    }
  }
  paths.resize(tail);

  return unique;
}

//------------------------------------------------------------------------------

size_type
GCSA::prefixDoubling(std::vector<PathNode<>>& paths, size_type kmer_length)
{
  bool fully_sorted = false;
  size_type path_order = 1;
  for(size_type step = 1; step <= DOUBLING_STEPS; step++)
  {
    size_type old_paths = paths.size();
    std::cout << "Doubling step " << step << " (path length " << (path_order * kmer_length) << " -> "
              << (2 * path_order * kmer_length) << ")" << std::endl;
    PathFromComparator<> pfc; // Sort the previous generation by from.
    parallelQuickSort(paths.begin(), paths.end(), pfc);

    joinPaths(paths, path_order);
    size_type joined_paths = paths.size();

    size_type unique = mergePaths(paths, path_order);
    std::cout << "  " << old_paths << " -> " << joined_paths << " -> "
              << paths.size() << " paths (" << unique << " with unique labels)" << std::endl;
    if(unique == paths.size()) { fully_sorted = true; break; }
  }
  this->max_query_length = (fully_sorted ? ~(size_type)0 : kmer_length << DOUBLING_STEPS);
  std::cout << "Max query length: " << this->order() << std::endl;

  return path_order;
}

void
GCSA::sample(std::vector<PathNode<>>& paths, size_type path_order)
{
  this->node_count = 0;
  this->sampled_nodes = bit_vector(paths.size(), 0);
  this->samples = bit_vector(paths.size(), 0);

  range_type range(1, 0);
  size_type sample_bits = 0;
  std::vector<node_type> sample_buffer;
  PathLabelComparator<> plc(path_order);
  while(true)
  {
    nextRange(range, paths, plc);
    if(range.first >= paths.size()) { break; }
    paths[this->node_count] = paths[range.first];
    bool sample_this = false;
    for(size_type i = range.first; i <= range.second; i++)
    {
      paths[this->node_count].addPredecessors(paths[i]);
      if(Node::offset(paths[i].from) == 0) { sample_this = true; }
    }
    if(paths[this->node_count].hasPredecessor(ENDMARKER_COMP)) { sample_this = true; }
    if(sdsl::bits::lt_cnt[paths[this->node_count].predecessors] > 1) { sample_this = true; }
    if(sample_this)
    {
      this->sampled_nodes[this->node_count] = 1;
      for(size_type i = range.first; i <= range.second; i++)
      {
        sample_bits = std::max(sample_bits, bit_length(paths[i].from));
        sample_buffer.push_back(paths[i].from);
      }
      this->samples[sample_buffer.size() - 1] = 1;
    }
    this->node_count++;
  }
  paths.resize(this->node_count);

  this->sampled_nodes.resize(this->node_count);
  sdsl::util::init_support(this->sampled_node_rank, &(this->sampled_nodes));

  this->samples.resize(sample_buffer.size());
  sdsl::util::init_support(this->sample_select, &(this->samples));

  this->stored_samples = sdsl::int_vector<0>(sample_buffer.size(), 0, sample_bits);
  for(size_type i = 0; i < sample_buffer.size(); i++) { this->stored_samples[i] = sample_buffer[i]; }
  sdsl::util::clear(sample_buffer);

  std::cout << "Path nodes: " << this->node_count << std::endl;
  std::cout << "Samples: " << this->stored_samples.size() << " values ("
            << (size_type)(this->stored_samples.width()) << " bits each)" << std::endl;
}

//------------------------------------------------------------------------------

range_type
GCSA::find(const char_type* pattern, size_type length) const
{
  return this->find(pattern, pattern + length);
}

void
GCSA::locate(size_type node, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { sdsl::util::clear(results); }
  if(node >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  this->locateInternal(node, results);
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
GCSA::locateInternal(size_type node, std::vector<node_type>& results) const
{
  size_type steps = 0;
  while(this->sampled_nodes[node] == 0)
  {
    node = this->LF(node);
    steps++;
  }

  range_type sample_range = this->sampleRange(node);
  for(size_type i = sample_range.first; i <= sample_range.second; i++)
  {
    results.push_back(this->stored_samples[i] + steps);
  }
}

//------------------------------------------------------------------------------

} // namespace gcsa
