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

#ifndef _GCSA_PATH_GRAPH_H
#define _GCSA_PATH_GRAPH_H

#include <sdsl/rmq_support.hpp>

#include "dbg.h"
#include "files.h"

namespace gcsa
{

/*
  path_graph.h: Internal graph manipulation methods.
*/

//------------------------------------------------------------------------------

struct PathLabel
{
  typedef std::uint32_t rank_type;

  // This should be at least 1 << GCSA::DOUBLING_STEPS.
  const static size_type LABEL_LENGTH = 8;

  // Labels starting with NO_RANK will be after real labels in lexicographic order.
  // We also use NO_RANK for padding last labels.
  const static rank_type NO_RANK = ~(rank_type)0;

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

  const static size_type LABEL_LENGTH = PathLabel::LABEL_LENGTH;

  node_type from, to;

  inline bool sorted() const { return (this->to == ~(node_type)0); }
  inline void makeSorted() { this->to = ~(node_type)0; }

//------------------------------------------------------------------------------

  /*
    From low-order to high-order bits:

    8 bits   which predecessor comp values exist
    4 bits   length of the kmer rank sequences representing the path label range
    4 bits   lcp of the above sequences
    8 bits   unused
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
  inline size_type order() const { return ((this->fields >> 8) & 0xF); }
  inline void setOrder(size_type new_order)
  {
    this->fields &= ~(size_type)0xF00;
    this->fields |= new_order << 8;
  }

  // LCP is the length of the common prefix of kmer rank sequences.
  inline size_type lcp() const { return ((this->fields >> 12) & 0xF); }
  inline void setLCP(size_type new_lcp)
  {
    this->fields &= ~(size_type)0xF000;
    this->fields |= new_lcp << 12;
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

  inline rank_type firstLabel(size_type i, const rank_type* labels) const
  {
    return labels[this->pointer() + i];
  }

  template<class ArrayType>
  inline rank_type lastLabel(size_type i, ArrayType& labels) const
  {
    if(i < this->lcp()) { return labels[this->pointer() + i]; }
    else { return labels[this->pointer() + this->order()]; }
  }

  inline rank_type lastLabel(size_type i, const rank_type* labels) const
  {
    if(i < this->lcp()) { return labels[this->pointer() + i]; }
    else { return labels[this->pointer() + this->order()]; }
  }

//------------------------------------------------------------------------------

  static std::vector<rank_type> dummyRankVector();

  PathNode(const KMer& kmer, std::vector<rank_type>& labels);

  PathNode(const PathNode& source,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  PathNode(const PathNode& left, const PathNode& right,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  /*
    The constructors set the pointers when loading the path node:
    - the first one sets it to labels.size()
    - the second one sets it to 0

    The third serialize() version sets the pointer to ptr. It can be used to write
    PathNodes into memory mappable files.
  */
  PathNode(std::istream& in, std::vector<rank_type>& labels);
  PathNode(std::istream& in, rank_type* labels);
  void serialize(std::ostream& out, const std::vector<rank_type>& labels) const;
  void serialize(std::ostream& out, const rank_type* labels) const;
  void serialize(std::ostream& node_stream, std::ostream& label_stream, const rank_type* labels, size_type ptr);

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
  typedef sdsl::rmq_succinct_sada<>       rmq_type; // Faster than rmq_support_sct.
  typedef PathNode::rank_type             rank_type;
  typedef std::pair<rank_type, rank_type> rank_range;

  size_type           kmer_length, total_keys;
  sdsl::int_vector<0> kmer_lcp;
  rmq_type            lcp_rmq;

  LCP();
  LCP(const std::vector<key_type>& keys, size_type _kmer_length);

  /*
    Computes the minimal/maximal lcp of the path labels corresponding to path nodes a and b.
    a must be before b in lexicographic order, and the ranges must not overlap.
    The returned lcp value is a pair (x,y), where x is the lcp of the PathNode labels
    and y is the lcp of the first diverging kmers.

    FIXME Later: Do not use the rmq if the kmer ranks are close.
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
    if(lcp.second + 1 < this->kmer_length) { lcp.second++; }
    else { lcp.first++; lcp.second = 0; }
    return lcp;
  }

  /*
    Extends the given rank range into a maximal range having the given lcp.

    FIXME Later: Build an LCP interval tree to do this faster.
  */
  inline rank_range extendRange(rank_range range, size_type lcp) const
  {
    while(range.first > 0 && this->kmer_lcp[range.first] >= lcp) { range.first--; }
    while(range.second + 1 < this->total_keys && this->kmer_lcp[range.second + 1] >= lcp) { range.second++; }
    return range;
  }

  void swap(LCP& another);
};

//------------------------------------------------------------------------------

/*
  A path graph is a set of files, each of them containing the paths derived from one
  of the input files. The PathNodes in each file are sorted by their labels, and the
  read() member functions will also return the PathNodes in sorted order.
*/

struct PathGraph
{
  std::vector<std::string> filenames;
  std::vector<size_type>   sizes, rank_counts;

  size_type path_count, rank_count;
  size_type order;

  size_type unique, unsorted, nondeterministic;

  const static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // .gcsa

  PathGraph(const InputGraph& source, sdsl::sd_vector<>& key_exists);
  PathGraph(size_type file_count, size_type path_order);
  ~PathGraph();

  void clear();
  void swap(PathGraph& another);
  void open(std::ifstream& input, size_type file) const;

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type k() const { return this->order; }
  inline size_type files() const { return this->filenames.size(); }

  inline size_type bytes() const
  {
    return this->size() * sizeof(PathNode) + this->ranks() * sizeof(PathNode::rank_type);
  }

  // Size limits are in bytes.
  void prune(const LCP& lcp, size_type size_limit);
  void extend(size_type size_limit);

  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file) const;

  PathGraph(const PathGraph&) = delete;
  PathGraph& operator= (const PathGraph&) = delete;
};

//------------------------------------------------------------------------------

/*
  A merged graph is a path graph with the path nodes in one file, the labels in another
  file, and the additional from nodes in a third file. The PathNodes are sorted by their
  labels and their label pointers are set correctly.
*/

struct MergedGraph
{
  std::string path_name, rank_name, from_name;

  size_type path_count, rank_count, from_count;
  size_type order;

  std::vector<size_type> next;  // paths[next[comp]] is the first path starting with comp.
  std::vector<size_type> next_from; // Where to find the corresponding additional from nodes.

  const static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // .gcsa

  MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper);
  ~MergedGraph();

  void clear();

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type extra() const { return this->from_count; }
  inline size_type k() const { return this->order; }

  inline size_type path_bytes() const { return this->size() * sizeof(PathNode); }
  inline size_type rank_bytes() const { return this->ranks() * sizeof(PathNode::rank_type); }
  inline size_type from_bytes() const { return this->extra() * sizeof(range_type); }

  inline size_type bytes() const
  {
    return this->path_bytes() + this->rank_bytes() + this->from_bytes();
  }

  MergedGraph(const MergedGraph&) = delete;
  MergedGraph& operator= (const MergedGraph&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_UTILS_H
