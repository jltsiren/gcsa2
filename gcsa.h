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
    The high-level interface deals with node ranges and actual character values.
    locate() stores the node identifiers in the given vector in sorted order.
    If append is true, the results are appended to the existing vector.
    If sort is true, the results are sorted and the duplicates are removed.

    FIXME implement
    FIXME find() based on iterators?
  */

  range_type find(const std::string& pattern) const;
  size_type locate(size_type node, std::vector<size_type>& results, bool append = false, bool sort = true) const;
  size_type locate(range_type range, std::vector<size_type>& results, bool append = false, bool sort = true) const;

//------------------------------------------------------------------------------

  /*
    The low-level interface deals with nodes, node ranges, and contiguous character values.
  */

  inline size_type size() const { return this->sampled_nodes.size(); }
  inline size_type edge_count() const { return this->bwt.size(); }

  inline size_type LF(size_type node, char_type comp) const
  {
    if(node > 0) { node = this->node_select(node) + 1; }  // The first incoming edge.
    node = gcsa::LF(this->bwt, this->alpha, node, comp);
    return this->edge_rank(node);
  }

  inline range_type LF(range_type range, char_type comp) const
  {
    if(range_first > 0) { range_first = this->node_select(range_first) + 1; }
    range.second = this->node_select(range_second + 1);
    range = gcsa::LF(this->bwt, this->alpha, range, comp);
    return range_type(this->edge_rank(range.first), this->edge_rank(range.second));
  }

//------------------------------------------------------------------------------

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
};  // class GCSA

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_GCSA_H
