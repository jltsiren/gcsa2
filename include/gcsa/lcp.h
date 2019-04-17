/*
  Copyright (c) 2019 Jouni Sirén
  Copyright (c) 2016 Genome Research Ltd.

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

#ifndef GCSA_LCP_H
#define GCSA_LCP_H

#include <gcsa/files.h>

namespace gcsa
{

/*
  lcp.h: LCP array that provides additional functionality to a GCSA index.
*/

//------------------------------------------------------------------------------

struct STNode
{
  size_type sp, ep;
  size_type left_lcp, right_lcp; // lcp[sp], lcp[ep + 1]
  size_type node_lcp;            // Can be UNKNOWN.

  constexpr static size_type UNKNOWN = ~(size_type)0;

  STNode() : sp(0), ep(0), left_lcp(0), right_lcp(0), node_lcp(0) {}

  STNode(size_type start, size_type end, size_type left, size_type right, size_type depth) :
    sp(start), ep(end),
    left_lcp(left), right_lcp(right),
    node_lcp(depth)
  {
  }

  inline range_type range() const { return range_type(this->sp, this->ep); }
  inline size_type lcp() const { return this->node_lcp; }

  inline bool operator== (const STNode& node) const
  {
    return (this->sp == node.sp && this->ep == node.ep);
  }

  inline bool operator== (range_type range) const
  {
    return (this->sp == range.first && this->ep == range.second);
  }

  inline bool operator!= (const STNode& node) const
  {
    return (this->sp != node.sp || this->ep != node.ep);
  }

  inline bool operator!= (range_type range) const
  {
    return (this->sp != range.first || this->ep != range.second);
  }
};

std::ostream& operator<< (std::ostream& out, const STNode& node);

//------------------------------------------------------------------------------

/*
  LCP array with support for some suffix tree operations using nsv/psv/rmq queries via
  a range minimum tree.
*/

class LCPArray
{
public:
  typedef gcsa::size_type size_type;
  typedef STNode          node_type;

//------------------------------------------------------------------------------

  LCPArray();
  LCPArray(const LCPArray& source);
  LCPArray(LCPArray&& source);
  ~LCPArray();

  void swap(LCPArray& another);
  LCPArray& operator=(const LCPArray& source);
  LCPArray& operator=(LCPArray&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  const static std::string EXTENSION; // .lcp

//------------------------------------------------------------------------------

  /*
    The InputGraph instance must be the same that has already been used for GCSA construction.
    The graph passes the temporary file containing the LCP array from GCSA construction to
    LCPArray construction.
  */

  explicit LCPArray(const InputGraph& graph, const ConstructionParameters& parameters = ConstructionParameters());

//------------------------------------------------------------------------------

  inline size_type size() const { return this->header.size; }
  inline size_type values() const { return this->data.size(); }
  inline size_type levels() const { return this->offsets.size() - 1; }
  inline size_type branching() const { return this->header.branching; }

  inline size_type operator[] (size_type i) const { return this->data[i]; }

//------------------------------------------------------------------------------

  /*
    High-level interface.
  */

  inline node_type root() const { return node_type(0, this->size() - 1, 0, 0, 0); }

  node_type parent(const node_type& node) const;
  node_type parent(range_type range) const;

  // String depth of the node; the length of the longest pattern matching the range.
  // Not supported for leaf nodes / ranges of length 1.
  size_type depth(const node_type& node) const;
  size_type depth(node_type& node) const;
  size_type depth(range_type range) const;

//------------------------------------------------------------------------------

  /*
    Low-level interface. The return value is a pair (res, LCP[res]) or notFound().
  */

  range_type psv(size_type pos) const;
  range_type psev(size_type pos) const;

  range_type nsv(size_type pos) const;
  range_type nsev(size_type pos) const;

  range_type rmq(size_type sp, size_type ep) const;
  range_type rmq(range_type range) const;

  inline node_type nodeFor(range_type range) const
  {
    if(range.second + 1 < this->size())
    {
      return node_type(range.first, range.second,
                       this->data[range.first], this->data[range.second + 1],
                       node_type::UNKNOWN);
    }
    else
    {
      return node_type(range.first, range.second, this->data[range.first], 0, node_type::UNKNOWN);
    }
  }

  // Returned when a psv/nsv/rmq query cannot find a suitable value.
  inline range_type notFound() const { return range_type(this->values(), this->values()); }

//------------------------------------------------------------------------------

  /*
    We store a k-ary range minimum tree over the LCP array. Each node is identified by
    its position in the data array. Level i occupies range [offsets[i], offsets[i + 1] - 1]
    in the array. Level 0 is the leaves (the LCP array).
  */

  LCPHeader            header;
  sdsl::int_vector<0>  data;
  sdsl::int_vector<64> offsets;

private:
  void copy(const LCPArray& source);
};  // class LCPArray

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_LCP_H
