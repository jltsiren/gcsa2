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

#include "files.h"

namespace gcsa
{

/*
  gcsa.h: The main public interface.
*/

//------------------------------------------------------------------------------

class GCSA
{
public:
  typedef gcsa::size_type  size_type;
  typedef sdsl::wt_huff<>  bwt_type;
  typedef sdsl::bit_vector bit_vector;

//------------------------------------------------------------------------------

  GCSA();
  GCSA(const GCSA& g);
  GCSA(GCSA&& g);
  ~GCSA();

  void swap(GCSA& g);
  GCSA& operator=(const GCSA& g);
  GCSA& operator=(GCSA&& g);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  const static std::string EXTENSION;       // .gcsa

  const static size_type DOUBLING_STEPS = 3;
  const static size_type SIZE_LIMIT = 500;  // Gigabytes.
  const static size_type ABSOLUTE_SIZE_LIMIT = 16384;

//------------------------------------------------------------------------------

  /*
    This is the main constructor. We build GCSA from the graph, doubling the path length
    the specified number of times (up to DOUBLING_STEPS). size_limit is the limit for the
    size of individual graph on disk in gigabytes. Actual disk usage may be up to 2x higher.
    If the graph was encoded using non-default Alphabet, an alphabet object must also be
    supplied.

    FIXME Later: We should pass an object containing construction parameters.
  */
  GCSA(const InputGraph& graph,
    size_type doubling_steps = DOUBLING_STEPS, size_type size_limit = SIZE_LIMIT,
    const Alphabet& _alpha = Alphabet());

  /*
    Index verification. The index is queried with all unique kmer labels in the input, and
    the list of occurrences is expected to be the same as the set of start nodes of the
    kmers with that label. To guarantee this, the input should be the same as for the
    constructor, or a subset where all kmers with a given label are either present or
    absent.

    Note: Verification uses multiple threads and sorts the kmer array.

    Returns false if the index fails verification, and true otherwise.
  */
  bool verifyIndex(std::vector<KMer>& kmers, size_type kmer_length) const;
  bool verifyIndex(const InputGraph& graph) const;

//------------------------------------------------------------------------------

  /*
    The high-level interface deals with path ranges and actual character values.
    locate() stores the node identifiers in the given vector in sorted order.
    If append is true, the results are appended to the existing vector.
    If sort is true, the results are sorted and the duplicates are removed.

    The implementation of find() is based on bidirectional iterators.

    If the pattern is longer than the order of the index, there may be false
    positives (but no false negatives). The results of such queries must be
    verified in the original graph.
  */

  template<class Iterator>
  range_type find(Iterator begin, Iterator end) const
  {
    if(begin == end) { return range_type(0, this->size() - 1); }

    --end;
    range_type range = this->charRange(this->alpha.char2comp[*end]);
    while(!Range::empty(range) && end != begin)
    {
      --end;
      range = this->LF(range, this->alpha.char2comp[*end]);
    }

    return range;
  }

  template<class Container>
  range_type find(const Container& pattern) const
  {
    return this->find(pattern.begin(), pattern.end());
  }

  template<class Element>
  range_type find(const Element* pattern, size_type length) const
  {
    return this->find(pattern, pattern + length);
  }

  void locate(size_type path, std::vector<node_type>& results, bool append = false, bool sort = true) const;
  void locate(range_type range, std::vector<node_type>& results, bool append = false, bool sort = true) const;

//------------------------------------------------------------------------------

  /*
    The low-level interface deals with paths, path ranges, and contiguous character values.
  */

  inline size_type size() const { return this->path_node_count; }
  inline size_type edge_count() const { return this->bwt.size(); }
  inline size_type order() const { return this->max_query_length; }

  inline bool has_samples() const { return (this->stored_samples.size() > 0); }
  inline size_type sample_count() const { return this->stored_samples.size(); }
  inline size_type sample_bits() const { return this->stored_samples.width(); }
  inline size_type sampled_positions() const { return this->sampled_path_rank(this->sampled_paths.size()); }

  inline range_type charRange(comp_type comp) const
  {
    return this->pathNodeRange(gcsa::charRange(this->alpha, comp));
  }

  inline range_type LF(range_type range, comp_type comp) const
  {
    range = this->bwtRange(range);
    range = gcsa::LF(this->bwt, this->alpha, range, comp);
    if(Range::empty(range)) { return range; }
    return this->pathNodeRange(range);
  }

  // Follow the first edge backwards.
  inline size_type LF(size_type path_node) const
  {
    path_node = this->startPos(path_node);
    auto temp = this->bwt.inverse_select(path_node);
    path_node = this->alpha.C[temp.second] + temp.first;
    return this->edge_rank(path_node);
  }

//------------------------------------------------------------------------------

  /*
    We could have a bitvector for mapping between the node identifiers returned by locate()
    and the original (node, offset) pairs. Alternatively, if the maximum length of a label
    is reasonable, we can just use (node, offset) = (id / max_length, id % max_length).
  */

  size_type                 path_node_count;
  size_type                 max_query_length;

  bwt_type                  bwt;
  Alphabet                  alpha;

  // The last BWT position in each path is marked with an 1-bit.
  bit_vector                path_nodes;
  bit_vector::rank_1_type   path_rank;
  bit_vector::select_1_type path_select;

  // The last outgoing edge from each path is marked with an 1-bit.
  bit_vector                edges;
  bit_vector::rank_1_type   edge_rank;
  bit_vector::select_1_type edge_select;

  // Paths containing samples are marked with an 1-bit.
  bit_vector                sampled_paths;
  bit_vector::rank_1_type   sampled_path_rank;

  // The last sample belonging to the same path is marked with an 1-bit.
  sdsl::int_vector<0>       stored_samples;
  bit_vector                samples;
  bit_vector::select_1_type sample_select;

//------------------------------------------------------------------------------

private:
  void copy(const GCSA& g);
  void setVectors();
  void initSupport();

  void locateInternal(size_type path, std::vector<node_type>& results) const;

//------------------------------------------------------------------------------

  inline size_type startPos(size_type path_node) const
  {
    return (path_node > 0 ? this->path_select(path_node) + 1 : 0);
  }

  inline size_type endPos(size_type path_node) const
  {
    return this->path_select(path_node + 1);
  }

  inline range_type bwtRange(range_type path_node_range) const
  {
    path_node_range.first = this->startPos(path_node_range.first);
    path_node_range.second = this->endPos(path_node_range.second);
    return path_node_range;
  }

  inline range_type pathNodeRange(range_type outgoing_range) const
  {
    outgoing_range.first = this->edge_rank(outgoing_range.first);
    outgoing_range.second = this->edge_rank(outgoing_range.second);
    return outgoing_range;
  }

  inline range_type sampleRange(size_type path_node) const
  {
    path_node = this->sampled_path_rank(path_node);
    range_type sample_range;
    sample_range.first = (path_node > 0 ? this->sample_select(path_node) + 1 : 0);
    sample_range.second = this->sample_select(path_node + 1);
    return sample_range;
  }
};  // class GCSA

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_GCSA_H
