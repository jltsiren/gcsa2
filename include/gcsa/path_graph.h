/*
  Copyright (c) 2018, 2019 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#ifndef GCSA_PATH_GRAPH_H
#define GCSA_PATH_GRAPH_H

#include <gcsa/dbg.h>
#include <gcsa/files.h>
#include <gcsa/internal.h>

namespace gcsa
{

/*
  path_graph.h: Internal graph manipulation methods.
*/

//------------------------------------------------------------------------------

struct PathLabel
{
  typedef std::uint32_t rank_type;

  // This should be at least 1 << ConstructionParameters::MAX_STEPS.
  constexpr static size_type LABEL_LENGTH = 16;

  // Labels starting with NO_RANK will be after real labels in lexicographic order.
  // We also use NO_RANK for padding last labels.
  constexpr static rank_type NO_RANK = ~(rank_type)0;

  rank_type label[LABEL_LENGTH];
  bool      first;

  inline bool operator< (const PathLabel& another) const
  {
    for(size_type i = 0; i < LABEL_LENGTH; i++)
    {
      if(this->label[i] != another.label[i]) { return (this->label[i] < another.label[i]); }
    }
    return (this->first && !(another.first));
  }

  inline bool operator<= (const PathLabel& another) const
  {
    for(size_type i = 0; i < LABEL_LENGTH; i++)
    {
      if(this->label[i] != another.label[i]) { return (this->label[i] < another.label[i]); }
    }
    return (this->first || !(another.first));
  }
};

//------------------------------------------------------------------------------

/*
  The node type used during doubling. As in the original GCSA, from and to are nodes
  in the original graph, denoting a path as a semiopen range [from, to). If
  to == -1, the path will not be extended, because it already has a unique label.
  rank_type is the integer type used to store ranks of the original kmers.
  During edge generation, 'to' node will be used to store indegree and the outdegree.

  The rank sequences are stored in an std::vector at position 'pointer()'. The stored
  sequence consists of the first label ('order()' ranks) followed by one rank for the
  diverging last rank of the last label. If the first and the last ranks are identical,
  the last rank is a dummy value.

  There are also alternative versions of most operations using raw pointers to the
  label array.
*/

struct PathNode
{
  typedef PathLabel::rank_type rank_type;

  constexpr static size_type LABEL_LENGTH = PathLabel::LABEL_LENGTH;

  node_type from, to;

  inline bool sorted() const { return (this->to == ~(node_type)0); }
  inline void makeSorted() { this->to = ~(node_type)0; }

//------------------------------------------------------------------------------

  /*
    From low-order to high-order bits:

    8 bits   which predecessor comp values exist
    8 bits   length of the kmer rank sequences representing the path label range
    8 bits   lcp of the above sequences
    40 bits  pointer to the label data
  */
  size_type fields;

  inline byte_type predecessors() const { return (this->fields & 0xFF); }
  inline void setPredecessors(byte_type preds)
  {
    this->fields &= ~(size_type)0xFF;
    this->fields |= (size_type)preds;
  }
  inline bool hasPredecessor(comp_type comp) const
  {
    return (this->fields & (1 << comp));
  }
  inline void addPredecessors(const PathNode& another)
  {
    this->fields |= another.predecessors();
  }

  // Order is the length of the kmer rank sequences representing the path label range.
  inline size_type order() const { return ((this->fields >> 8) & 0xFF); }
  inline void setOrder(size_type new_order)
  {
    this->fields &= ~(size_type)0xFF00;
    this->fields |= new_order << 8;
  }

  // LCP is the length of the common prefix of kmer rank sequences.
  inline size_type lcp() const { return ((this->fields >> 16) & 0xFF); }
  inline void setLCP(size_type new_lcp)
  {
    this->fields &= ~(size_type)0xFF0000;
    this->fields |= new_lcp << 16;
  }

  inline size_type pointer() const { return (this->fields >> 24); }
  inline void setPointer(size_type new_pointer)
  {
    this->fields &= 0xFFFFFF;
    this->fields |= new_pointer << 24;
  }

//------------------------------------------------------------------------------

  inline size_type ranks() const { return this->order() + 1; }
  inline size_type bytes() const { return sizeof(*this) + this->ranks() * sizeof(rank_type); }

  template<class ArrayType>
  inline rank_type firstLabel(size_type i, ArrayType& labels) const
  {
    return labels[this->pointer() + i];
  }

  template<class ArrayType>
  inline rank_type lastLabel(size_type i, ArrayType& labels) const
  {
    if(i < this->lcp()) { return labels[this->pointer() + i]; }
    else { return labels[this->pointer() + this->order()]; }
  }

//------------------------------------------------------------------------------

  PathNode(const KMer& kmer, WriteBuffer<rank_type>& labels);

  PathNode(const PathNode& source,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  PathNode(const PathNode& left, const PathNode& right,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  void print(std::ostream& out, const std::vector<rank_type>& labels) const;
  void print(std::ostream& out, const rank_type* labels) const;

  PathNode();
  explicit PathNode(std::vector<rank_type>& labels);
  PathNode(PathNode&& source);
  ~PathNode();

  inline void swap(PathNode& another)
  {
    if(&another != this)
    {
      std::swap(this->from, another.from); std::swap(this->to, another.to);
      std::swap(this->fields, another.fields);
    }
  }

  PathNode& operator= (PathNode&& source);

  /*
    These are dangerous, because the nodes will share the same label. Changing one will
    change the other as well.
  */
  PathNode(const PathNode& source);
  PathNode& operator= (const PathNode& source);
  void copy(const PathNode& source);
};

// Compares the first labels.
struct PathFirstComparator
{
  const std::vector<PathNode::rank_type>& labels;

  explicit PathFirstComparator(const std::vector<PathNode::rank_type>& _labels) : labels(_labels) { }

  inline bool operator() (const PathNode& a, const PathNode& b) const
  {
    size_type order = std::min(a.order(), b.order());
    for(size_type i = 0, a_ptr = a.pointer(), b_ptr = b.pointer(); i < order; i++, a_ptr++, b_ptr++)
    {
      if(labels[a_ptr] != labels[b_ptr]) { return (labels[a_ptr] < labels[b_ptr]); }
    }
    return (a.order() < b.order());
  }
};

// Compares the 'from' nodes.
struct PathFromComparator
{
  inline bool operator() (const PathNode& a, const PathNode& b) const
  {
    return (a.from < b.from);
  }
};

struct FromGetter
{
  inline static size_type get(const PathNode& path)
  {
    return path.from;
  }
};

//------------------------------------------------------------------------------

struct LCP
{
  typedef PathNode::rank_type rank_type;

  size_type       kmer_length, total_keys;
  sdsl::wt_blcd<> kmer_lcp; // Faster than proper RMQ for small values.

  LCP();
  LCP(const std::vector<key_type>& keys, size_type _kmer_length);

  /*
    Computes the minimal/maximal lcp of the path labels corresponding to path nodes a and b.
    a must be before b in lexicographic order, and the ranges must not overlap.
    The returned lcp value is a pair (x,y), where x is the lcp of the PathNode labels
    and y is the lcp of the first diverging kmers.
  */
  range_type min_lcp(const PathNode& a, const PathNode& b, const std::vector<rank_type>& labels) const;
  range_type max_lcp(const PathNode& a, const PathNode& b, const std::vector<rank_type>& labels) const;

  range_type min_lcp(const PathNode& a, const PathNode& b,
    const rank_type* a_labels, const rank_type* b_labels) const;
  range_type max_lcp(const PathNode& a, const PathNode& b,
    const rank_type* a_labels, const rank_type* b_labels) const;

  // Increments the lcp by 1.
  inline range_type increment(range_type lcp) const
  {
    lcp.second++;
    lcp.first += lcp.second / this->kmer_length;
    lcp.second %= this->kmer_length;
    return lcp;
  }

  void swap(LCP& another);
};

//------------------------------------------------------------------------------

/*
  A path graph is a set of files. Each file in the input graph becomes two temporary
  files: one for the paths and another for the rank sequences corresponding to path
  labels. The PathNodes in each file are sorted by their labels, and the read() member
  functions will also return the PathNodes in sorted order. The labels are stored in
  the same order as the PathNodes.
*/

struct PathGraph
{
  std::vector<std::string> path_names, rank_names;
  std::vector<size_type>   path_counts, rank_counts;

  size_type path_count, rank_count, range_count;
  size_type order, doubling_steps;

  size_type unique, redundant, unsorted, nondeterministic;

  bool delete_files;

  constexpr static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // gcsa

  PathGraph(const InputGraph& source, sdsl::int_vector<0>& distinct_labels);
  PathGraph(size_type file_count, size_type path_order, size_type steps);
  PathGraph(const std::string& path_name, const std::string& rank_name);  // For debugging.
  ~PathGraph();

  void clear();
  void swap(PathGraph& another);

  void open(std::ifstream& path_file, std::ifstream& rank_file, size_type file) const;

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type ranges() const { return this->range_count; } // Only works after prune().
  inline size_type k() const { return this->order; }
  inline size_type step() const { return this->doubling_steps; }
  inline size_type files() const { return this->path_names.size(); }

  inline size_type bytes() const
  {
    return this->size() * sizeof(PathNode) + this->ranks() * sizeof(PathNode::rank_type);
  }

  /*
    The size limit (in bytes) is the space available for the new graph. The current graph is
    not taken into account, so you may want to use something like:

      path_graph.prune(lcp, total_size_limit - path_graph.bytes())
  */
  void prune(const LCP& lcp, size_type size_limit);
  void extend(size_type size_limit);

  void debugExtend();

  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file) const;

  PathGraph(const PathGraph&) = delete;
  PathGraph& operator= (const PathGraph&) = delete;
};

//------------------------------------------------------------------------------

/*
  A merged graph is a path graph with the path nodes in one file, the labels in another
  file, and the additional start nodes in a third file. The PathNodes are sorted by their
  labels and their label pointers are set correctly.
*/

struct MergedGraph
{
  std::string path_name, rank_name, from_name, lcp_name;

  size_type path_count, rank_count, from_count;
  size_type order;

  std::vector<size_type> next;      // paths[next[comp]] is the first path starting with comp.
  std::vector<size_type> next_from; // Where to find the corresponding additional start nodes.

  constexpr static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // gcsa

  /*
    The size limit (in bytes) is the space available for the new graph. The current graph is
    not taken into account, so you may want to use something like:

      MergedGraph merged_graph(source, mapper, kmer_lcp, total_size_limit - source.bytes())
  */
  MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper, const LCP& kmer_lcp, size_type size_limit);
  ~MergedGraph();

  void clear();

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type extra() const { return this->from_count; }
  inline size_type k() const { return this->order; }

  inline size_type path_bytes() const { return this->size() * sizeof(PathNode); }
  inline size_type rank_bytes() const { return this->ranks() * sizeof(PathNode::rank_type); }
  inline size_type from_bytes() const { return this->extra() * sizeof(range_type); }
  inline size_type lcp_bytes() const { return this->size(); }

  inline size_type bytes() const
  {
    return this->path_bytes() + this->rank_bytes() + this->from_bytes() + this->lcp_bytes();
  }

  MergedGraph(const MergedGraph&) = delete;
  MergedGraph& operator= (const MergedGraph&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_UTILS_H
