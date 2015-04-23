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

  inline static size_type encodeKMer(const Alphabet& alpha, const std::string& kmer,
    uint8_t _predecessors, uint8_t _successors, bool parallel = true)
  {
    size_type value = 0;
    for(size_type i = 0; i < kmer.size(); i++) { value = (value << 3) | alpha.char2comp[kmer[i]]; }
    value = (value << 8) | _predecessors;
    value = (value << 8) | _successors;
    return value;
  }

  inline static size_type kMer(size_type code) { return (code >> 16); }
  inline static uint8_t predecessors(size_type code) { return (code >> 8) & 0xFF; }
  inline static uint8_t successors(size_type code) { return code & 0xFF; }

  GCSA(const std::vector<uint64_t>& kmers, const Alphabet& _alpha);

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
    range.first = this->startPos(range.first);
    range.second = this->endPos(range.second);
    return range;
  }

  inline range_type nodeRange(range_type incoming_range) const
  {
    range.first = this->edge_rank(range.first);
    range.second = this->edge_rank(range.second);
    return range;
  }

  inline size_type sampleRange(size_type node) const
  {
    node = this->sampled_node_rank(node);
    range_type sample_range;
    sample_range.first = (node > 0 ? this->sample_select(node) + 1 : 0);
    sample_range.second = this->sample_select(node + 1);
    return sample_range;
  }
};  // class GCSA

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_GCSA_H
