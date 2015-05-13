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

  this->build<key_type, KeyGetter>(keys, _alpha, total_edges);
}

GCSA::GCSA(std::vector<KMer>& kmers, size_type kmer_length, const Alphabet& _alpha)
{
  if(kmers.size() == 0) { return; }

  // Sort the kmers, build the mapper GCSA for generating the edges.
  std::vector<key_type> keys;
  sdsl::int_vector<0> last_char;
  uniqueKeys(kmers, keys, last_char);
  GCSA mapper(keys, kmer_length, _alpha);
  sdsl::util::clear(keys);

  // Transform the kmers into PathNodes.
  // FIXME Later: Save memory by not having both KMers and PathNodes in memory.
  std::vector<PathNode> paths(kmers.size());
  for(size_type i = 0; i < kmers.size(); i++) { paths[i] = PathNode(kmers[i]); }
  sdsl::util::clear(kmers);

  // Build the GCSA in PathNodes.
  std::vector<PathNode> last_labels;
  size_type path_order = this->prefixDoubling(paths, kmer_length, last_labels);
  std::vector<range_type> from_nodes;
  this->mergeByLabel(paths, path_order, from_nodes, last_labels);
  size_type total_edges = this->countEdges(paths, path_order, last_labels, mapper, last_char);
  sdsl::util::clear(last_labels);
  sdsl::util::clear(mapper); sdsl::util::clear(last_char);
  this->build<PathNode, PathNodeGetter>(paths, _alpha, total_edges);

  this->sample(paths, from_nodes);
  sdsl::util::clear(from_nodes);
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
  sdsl::bit_vector                values;     // Marks the values that are present.
  sdsl::bit_vector::rank_1_type   value_rank;

  sdsl::bit_vector                first_occ;  // Marks the first occurrence of each rank.
  sdsl::bit_vector::select_1_type first_select;

  enum Field { field_from, field_label };

  ValueIndex(const std::vector<ValueType>& input)
  {
    this->values = sdsl::bit_vector(Getter::get(input[input.size() - 1]) + 1, 0);
    this->first_occ = sdsl::bit_vector(input.size(), 0);

    size_type prev = ~(size_type)0;
    for(size_type i = 0; i < input.size(); i++)
    {
      size_type curr = Getter::get(input[i]);
      if(curr != prev)
      {
        this->values[curr] = 1;
        this->first_occ[i] = 1;
        prev = curr;
      }
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
  Join paths by left.to == right.from.
  FIXME Later: Predict the size of the next generation or write it to disk.
  FIXME Later: parallelize
*/
void
joinPaths(std::vector<PathNode>& paths, std::vector<PathNode>& last_labels)
{
  PathFromComparator pfc; // Sort the paths by from.
  parallelQuickSort(paths.begin(), paths.end(), pfc);
  size_type old_path_count = paths.size(), old_multilabel = last_labels.size();

  ValueIndex<PathNode, FromGetter> from_index(paths);
  std::vector<PathNode> next, new_last;
  for(size_type i = 0; i < paths.size(); i++)
  {
    if(paths[i].sorted())
    {
      next.push_back(paths[i]);
      if(paths[i].multiLabel())
      {
        new_last.push_back(last_labels[paths[i].lastLabel()]);
        next[next.size() - 1].setLastLabel(new_last.size() - 1);
      }
      continue;
    }
    size_type first = from_index.find(paths[i].to);
    for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++)
    {
      next.push_back(PathNode(paths[i], paths[j]));
      if(paths[j].multiLabel())
      {
        new_last.push_back(PathNode(paths[i], last_labels[paths[j].lastLabel()]));
        next[next.size() - 1].setLastLabel(new_last.size() - 1);
      }
    }
  }

  paths.swap(next); last_labels.swap(new_last);
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "  joinPaths(): " << old_path_count << " -> " << paths.size() << " paths ("
            << old_multilabel << " -> " << last_labels.size() << " multilabel)" << std::endl;
#endif
}

/*
  Iterates through ranges of paths with the same label. Empty range means start from
  the beginning, while range.first >= paths.size() means the end. Returns true if
  the range contains only nodes with the same label. In that case the range may cover
  multiple adjacent labels.
*/
bool
nextRange(range_type& range, const std::vector<PathNode>& paths,
  const PathLabelComparator& plc)
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
  labels sorted. Returns (unique_labels, unsorted_paths).
*/
size_type
mergePaths(std::vector<PathNode>& paths, size_type path_order, std::vector<PathNode>& last_labels)
{
  PathLabelComparator plc(path_order);
  parallelQuickSort(paths.begin(), paths.end(), plc); // Sort paths by labels.
  size_type old_path_count = paths.size(), old_multilabel = last_labels.size();

  std::vector<PathNode> new_last;
  size_type tail = 0, unique = 0, unsorted = 0;
  range_type range(1, 0);
  while(true)
  {
    bool same_from = nextRange(range, paths, plc);
    if(range.first >= paths.size()) { break; }
    if(same_from)
    {
if(paths[range.first].from == 58910 * 256 + 38)
{
  std::cout << "At " << range << std::endl;
}
      paths[tail] = paths[range.first];
      paths[tail].makeSorted();
      PathNode last = paths[range.second];
      for(size_type i = range.first; i <= range.second; i++)
      {
if(paths[range.first].from == 58910 * 256 + 38)
{
  std::cout << "  found " << paths[i] << std::endl;
}
        if(paths[i].multiLabel() && plc(last, last_labels[paths[i].lastLabel()]))
        {
          last = last_labels[paths[i].lastLabel()];
        }
      }
      if(plc(paths[tail], last))
      {
        new_last.push_back(last);
        paths[tail].setLastLabel(new_last.size() - 1);
      }
if(paths[range.first].from == 58910 * 256 + 38)
{
  std::cout << "  produced " << paths[tail] << " at " << tail << std::endl;
  if(paths[tail].multiLabel())
  {
    std::cout << "  last label " << new_last[paths[tail].lastLabel()] << std::endl;
  }
}
      tail++; unique++;
    }
    else
    {
      for(size_type i = range.first; i <= range.second; i++)
      {
        if(!(paths[i].sorted())) { unsorted++; }
        paths[tail] = paths[i];
        if(paths[tail].multiLabel())
        {
          new_last.push_back(last_labels[paths[tail].lastLabel()]);
          paths[tail].setLastLabel(new_last.size() - 1);
std::cout << "This should not happen." << std::endl;
        }
        tail++;
      }
    }
  }

  paths.resize(tail); last_labels.swap(new_last);
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "  mergePaths(): " << old_path_count << " -> " << paths.size() << " paths ("
            << old_multilabel << " -> " << last_labels.size() << " multilabel)" << std::endl;
  std::cerr << "  mergePaths(): " << unique << " unique paths, " << unsorted << " unsorted paths" << std::endl;
#endif
  return unsorted;
}

size_type
GCSA::prefixDoubling(std::vector<PathNode>& paths, size_type kmer_length, std::vector<PathNode>& last_labels)
{
  bool fully_sorted = false;
  size_type path_order = 1;
  for(size_type step = 1; step <= DOUBLING_STEPS; step++)
  {
#ifdef VERBOSE_STATUS_INFO
    std::cerr << "GCSA::prefixDoubling(): Step " << step << " (path length " << (path_order * kmer_length) << " -> "
              << (2 * path_order * kmer_length) << ")" << std::endl;
#endif
    joinPaths(paths, last_labels); path_order *= 2;
    size_type unsorted = mergePaths(paths, path_order, last_labels);
    if(unsorted == 0) { fully_sorted = true; break; }
  }
  this->max_query_length = (fully_sorted ? ~(size_type)0 : kmer_length << DOUBLING_STEPS);

  return path_order;
}

//------------------------------------------------------------------------------

void
GCSA::mergeByLabel(std::vector<PathNode>& paths, size_type path_order,
  std::vector<range_type>& from_nodes, std::vector<PathNode>& last_labels)
{
  this->path_node_count = 0;
  size_type old_path_count = paths.size(), old_multilabel = last_labels.size();

  range_type range(1, 0);
  PathLabelComparator plc(path_order);
  std::vector<PathNode> new_last;
  while(true)
  {
    nextRange(range, paths, plc);
    if(range.first >= paths.size()) { break; }
if(paths[range.first].label[0] == 84464)
{
  std::cout << "Found " << paths[range.first] << " at " << range << std::endl;
}
    paths[this->path_node_count] = paths[range.first];
    PathNode last =
      (paths[range.first].multiLabel() ? last_labels[paths[range.first].lastLabel()] : paths[range.first]);
    for(size_type i = range.first + 1; i <= range.second; i++)
    {
      paths[this->path_node_count].addPredecessors(paths[i]);
      from_nodes.push_back(range_type(this->path_node_count, paths[i].from));
      if(paths[i].multiLabel() && plc(last, last_labels[paths[i].lastLabel()]))
      {
        last = last_labels[paths[i].lastLabel()];
      }
    }
    if(plc(paths[this->path_node_count], last))
    {
      new_last.push_back(last);
      paths[this->path_node_count].setLastLabel(new_last.size() - 1);
    }
if(paths[range.first].label[0] == 84464)
{
  std::cout << "  produced " << paths[this->path_node_count] << " at " << this->path_node_count << std::endl;
  if(paths[this->path_node_count].multiLabel())
  {
    std::cout << "  last label " << new_last[paths[this->path_node_count].lastLabel()] << std::endl;
  }
}
    this->path_node_count++;
  }

  paths.resize(this->path_node_count); last_labels.swap(new_last);
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::mergeByLabel(): " << old_path_count << " -> " << paths.size() << " paths ("
            << old_multilabel << " -> " << last_labels.size() << " multilabel)" << std::endl;
#endif

for(size_type i = 0, next = 0; i < paths.size(); i++)
{
  if(paths[i].label[0] != next)
  {
    std::cout << "Node " << i << ": " << paths[i] << ": expected " << next << std::endl;
  }
  next = paths[i].label[0] + 1;
  if(paths[i].multiLabel()) { next = last_labels[paths[i].lastLabel()].label[0] + 1; }
}
}

//------------------------------------------------------------------------------

/*
  Path node labels are padded with 0s. They are interpreted as lower bounds of path labels starting
  from the path node. Predecessor labels are padded with ~0s and interpreted as upper bounds of the
  path labels starting with the given edge. The predecessor found using mapper.LF() matches
  the last path node with path_node.label <= predecessor.label. If a merged node originally had
  multiple labels, the LF() is based on the lexicographically largest of them.

  FIXME Later: parallelize
*/
size_type
GCSA::countEdges(std::vector<PathNode>& paths, size_type path_order,
  std::vector<PathNode>& last_labels,
  const GCSA& mapper, const sdsl::int_vector<0>& last_char)
{
  for(size_type i = 0; i < paths.size(); i++) { paths[i].to = 0; }
  PathLabelComparator plc(path_order);

  // Pointers to the next path nodes with labels starting with the given comp value.
  size_type sigma = mapper.alpha.sigma;
  size_type next[sigma + 1];
  for(size_type comp = 0; comp < sigma; comp++) { next[comp] = mapper.charRange(comp).first; }
  next[sigma] = ~(size_type)0;
  for(size_type i = 0, comp = 0; i < paths.size(); i++)
  {
    while(paths[i].label[0] >= next[comp]) { next[comp] = i; comp++; }
  }

  size_type total_edges = 0;
  for(size_type i = 0; i < paths.size(); i++)
  {
if(i >= 314891 && i <= 314895)
{
  std::cout << "At " << i << std::endl;
  std::cout << "  Path " << paths[i] << std::endl;
  if(paths[i].multiLabel()) { std::cout << "  Last label " << last_labels[paths[i].lastLabel()] << std::endl; }
  std::cout << "  LF(" << paths[i].label[0] << ", 1) = " << mapper.LF(paths[i].label[0], 1) << std::endl;
}
    size_type order = paths[i].order();
    for(size_type comp = 0; comp < sigma; comp++)
    {
      if(!(paths[i].hasPredecessor(comp))) { continue; }
      PathNode predecessor; predecessor.setOrder(path_order);
      PathNode& curr = (paths[i].multiLabel() ? last_labels[paths[i].lastLabel()] : paths[i]);
      predecessor.label[0] = mapper.LF(curr.label[0], comp);
      for(size_type j = 1; j < order; j++)
      {
        predecessor.label[j] = mapper.LF(curr.label[j], last_char[curr.label[j - 1]]);
      }
      if(order < path_order)
      {
        predecessor.label[order] = mapper.alpha.C[comp + 1] - 1;
        for(size_type j = order + 1; j < PathNode::LABEL_LENGTH; j++)
        {
          predecessor.label[j] = ~(PathNode::rank_type)0;
        }
      }

      // next[comp] should point to the last path node with path_node.label <= predecessor.label.
      while(next[comp] + 1 < paths.size() && !plc(predecessor, paths[next[comp] + 1]))
      {
        next[comp]++;
      }
      paths[next[comp]].to++; total_edges++;

if(predecessor.label[0] >= 84463 && predecessor.label[0] <= 84465)
{
  std::cout << "  Predecessor(" << comp << ") = " << predecessor << std::endl;
  std::cout << "  Maps to " << paths[next[comp]] << " at " << next[comp] << std::endl;
  if(next[comp] > 0) { std::cout << "    previous " << paths[next[comp] - 1] << std::endl; }
  if(next[comp] + 1 < paths.size()) { std::cout << "    next " << paths[next[comp] + 1] << std::endl; }
}
    }
  }

for(size_type i = 0; i < paths.size(); i++)
{
  if(paths[i].outdegree() == 0)
  {
    std::cerr << "Path " << i << " " << paths[i] << " has outdegree 0" << std::endl;
  }
}

  return total_edges;
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
  bit_vector new_sampled_nodes(paths.size(), 0);
  this->sampled_paths = bit_vector(paths.size(), 0);
  this->samples = bit_vector(paths.size() + from_nodes.size(), 0);

  size_type sample_bits = 0;
  std::vector<node_type> sample_buffer;
  ValueIndex<range_type, FirstGetter> from_index(from_nodes);
  for(size_type i = 0, j = 0; i < paths.size(); i++)
  {
    bool sample_this = false;
    std::vector<node_type> curr = fromNodes(i, paths, j, from_nodes);
    if(sdsl::bits::lt_cnt[paths[i].predecessors()] > 1) { sample_this = true; }
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
