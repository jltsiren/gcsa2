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

#ifndef _GCSA_GCSA_H
#define _GCSA_GCSA_H

#include "support.h"

namespace gcsa
{

//------------------------------------------------------------------------------

class GCSA
{
public:
  typedef gcsa::size_type size_type;
  typedef wt_huff<>       bwt_type;

//------------------------------------------------------------------------------

  GCSA();
  GCSA(const GCSA& g);
  GCSA(GCSA&& g);
  ~GCSA();

  void swap(GCSA& g);
  GCSA& operator=(const GCSA& g);
  GCSA& operator=(GCSA&& g);

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  const static std::string EXTENSION;  // .gcsa

//------------------------------------------------------------------------------

  /*
    This interface is intended for indexing 16-mers on an alphabet of size 8 or less. The
    16-mer is encoded as a 64-bit integer (most significant bit first):
      - 16x3 bits for the 16-mer
      - 8 bits for marking which predecessors are present
      - 8 bits for marking which successors are present

    The input vector must be sorted and contain only unique 16-mers. Each 16-mer must have
    at least one predecessor and one successor.
  */

  inline static size_type encode(const Alphabet& alpha, const std::string& label,
    uint8_t pred, uint8_t succ)
  {
    size_type value = 0;
    for(size_type i = 0; i < label.length(); i++) { value = (value << 3) | alpha.char2comp[label[i]]; }
    value = (value << 8) | pred;
    value = (value << 8) | succ;
    return value;
  }

  inline static size_type kmer(size_type code) { return (code >> 16); }
  inline static uint8_t predecessors(size_type code) { return (code >> 8) & 0xFF; }
  inline static uint8_t successors(size_type code) { return code & 0xFF; }

  inline static size_type merge(size_type code1, size_type code2) { return code1 | (code2 & 0xFFFF); }

  explicit GCSA(const std::vector<uint64_t>& kmers, const Alphabet& _alpha = Alphabet());

//------------------------------------------------------------------------------

  /*
    The high-level interface deals with node ranges and actual character values.
    locate() stores the node identifiers in the given vector in sorted order.
    If append is true, the results are appended to the existing vector.
    If sort is true, the results are sorted and the duplicates are removed.

    FIXME find() based on iterators?
  */

  range_type find(const std::string& pattern) const;
  void locate(size_type node, std::vector<size_type>& results, bool append = false, bool sort = true) const;
  void locate(range_type range, std::vector<size_type>& results, bool append = false, bool sort = true) const;

//------------------------------------------------------------------------------

  /*
    The low-level interface deals with nodes, node ranges, and contiguous character values.
  */

  inline size_type size() const { return this->node_count; }
  inline size_type edge_count() const { return this->bwt.size(); }

  inline bool has_samples() const { return (this->stored_samples.size() > 0); }

  inline size_type LF(size_type node, char_type comp) const
  {
    node = this->startPos(node);
    node = gcsa::LF(this->bwt, this->alpha, node, comp);
    return this->edge_rank(node);
  }

  inline range_type LF(range_type range, char_type comp) const
  {
    range = this->bwtRange(range);
    range = gcsa::LF(this->bwt, this->alpha, range, comp);
    return this->nodeRange(range);
  }

  // Follow the first edge backwards.
  inline size_type LF(size_type node) const
  {
    node = this->startPos(node);
    auto temp = this->bwt.inverse_select(node);
    node = this->alpha.C[temp.second] + temp.first;
    return this->edge_rank(node);
  }

//------------------------------------------------------------------------------

  /*
    We could have a bitvector for mapping between the node identifiers returned by locate()
    and the original (node, offset) pairs. Alternatively, if the maximum length of a label
    is reasonable, we can just use (node, offset) = (id / max_length, id % max_length).
  */

  size_type                 node_count;

  bwt_type                  bwt;
  Alphabet                  alpha;

  // The last BWT position in each node is marked with an 1-bit.
  bit_vector                nodes;
  bit_vector::rank_1_type   node_rank;
  bit_vector::select_1_type node_select;

  // The last outgoing edge from each node is marked with an 1-bit.
  bit_vector                edges;
  bit_vector::rank_1_type   edge_rank;
  bit_vector::select_1_type edge_select;

  // Nodes containing samples are marked with an 1-bit.
  bit_vector                sampled_nodes;
  bit_vector::rank_1_type   sampled_node_rank;

  // The last sample belonging to the same node is marked with an 1-bit.
  int_vector<0>             stored_samples;
  bit_vector                samples;
  bit_vector::select_1_type sample_select;

//------------------------------------------------------------------------------

private:
  void copy(const GCSA& g);
  void setVectors();

  void locateInternal(size_type node, std::vector<size_type>& results) const;

//------------------------------------------------------------------------------

  inline size_type startPos(size_type node) const
  {
    return (node > 0 ? this->node_select(node) + 1 : 0);
  }

  inline size_type endPos(size_type node) const
  {
    return this->node_select(node + 1);
  }

  inline range_type bwtRange(range_type node_range) const
  {
    node_range.first = this->startPos(node_range.first);
    node_range.second = this->endPos(node_range.second);
    return node_range;
  }

  inline range_type nodeRange(range_type incoming_range) const
  {
    incoming_range.first = this->edge_rank(incoming_range.first);
    incoming_range.second = this->edge_rank(incoming_range.second);
    return incoming_range;
  }

  inline range_type sampleRange(size_type node) const
  {
    node = this->sampled_node_rank(node);
    range_type sample_range;
    sample_range.first = (node > 0 ? this->sample_select(node) + 1 : 0);
    sample_range.second = this->sample_select(node + 1);
    return sample_range;
  }
};  // class GCSA

//------------------------------------------------------------------------------

/*
  The node type used during doubling. As in the original GCSA, from and to are node
  ids in the original graph, denoting a path as a semiopen range [from, to). If
  from == to, the path will not be extended, because it already has a unique label.
  rank_type is the integer type used to store ranks of the original kmers.
  During edge generation, to will be used to store the number of outgoing edges.

  FIXME Move to .cpp?
*/

template<class rank_type = uint32_t, size_type label_length = 8>
struct DoublingNode
{
  size_type from, to;
  rank_type label[label_length];

  inline bool operator< (const DoublingNode& another) const
  {
    if(this->from != another.from) { return (this->from < another.from); }
    return (this->to < another.to);
  }

  inline bool sorted() const { return (this->from == this->to); }
  inline size_type outdegree() const { return this->to; }

  DoublingNode(size_type _from, size_type _to, size_type rank)
  {
    this->from = _from; this->to = _to;
    this->label[0] = rank;
    for(size_type i = 1; i < label_length; i++) { label[i] = 0; }
  }

  DoublingNode(const DoublingNode& left, const DoublingNode& right, size_type order)
  {
    this->from = left.from; this->to = right.to;
    for(size_type i = 0; i < order; i++) { this->label[i] = left.label[i]; }
    for(size_type i = order; i < 2 * order; i++) { this->label[i] = right.label[i - order]; }
    for(size_type i = 2 * order; i < label_length; i++) { this->label[i] = 0; }
  }

  explicit DoublingNode(std::ifstream& in)
  {
    read_member(&(this->from), in); read_member(&(this->to), in);
    in.read((char*)(this->label), label_length * sizeof(rank_type));
  }

  size_type serialize(std::ostream& out) const
  {
    size_type bytes = 0;
    bytes += write_member(this->from, out); bytes += write_member(this->to, out);
    out.write((char*)(this->label), label_length * sizeof(rank_type));
    bytes += label_length * sizeof(rank_type);
    return bytes;
  }

//------------------------------------------------------------------------------

  DoublingNode() {}
  DoublingNode(const DoublingNode& source) { this->copy(source); }
  DoublingNode(DoublingNode&& source) { *this = std::move(source); }
  ~DoublingNode() {}

  void copy(const DoublingNode& source)
  {
    if(&source != this)
    {
      this->from = source.from; this->to = source.to;
      for(size_type i = 0; i < label_length; i++) { this->label[i] = source.label[i]; }
    }
  }

  void swap(DoublingNode& another)
  {
    if(&another != this)
    {
      std::swap(this->from, another.from); std::swap(this->to, another.to);
      for(size_type i = 0; i < label_length; i++) { std::swap(this->label[i], another.label[i]); }
    }
  }

  DoublingNode& operator= (const DoublingNode& source)
  {
    this->copy(source);
    return *this;
  }

  DoublingNode& operator= (DoublingNode&& source)
  {
    if(&source != this)
    {
      this->from = std::move(source.from);
      this->to = std::move(source.to);
      for(size_type i = 0; i < label_length; i++) { this->label[i] = std::move(source.label[i]); }
    }
    return *this;
  }
};

template<class rank_type = uint32_t, size_type label_length = 8>
struct NodeLabelComparator
{
  typedef DoublingNode<rank_type, label_length> node_type;

  size_type max_length;

  NodeLabelComparator(size_type len = label_length) : max_length(len) {}

  inline bool operator() (const node_type& a, const node_type& b) const
  {
    for(size_t i = 0; i < this->max_length; i++)
    {
      if(a.label[i] != b.label[i]) { return (a.label[i] < b.label[i]); }
    }
    return false;
  }
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_GCSA_H
