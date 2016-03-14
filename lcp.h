/*
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

#ifndef _GCSA_LCP_H
#define _GCSA_LCP_H

#include "files.h"

namespace gcsa
{

/*
  dbg.h: Specialization of GCSA for de Bruijn graphs.
*/

//------------------------------------------------------------------------------

struct STNode
{
  size_type sp, ep;
  size_type left_lcp, right_lcp; // lcp[sp], lcp[ep + 1]

  STNode() : sp(0), ep(0), left_lcp(0), right_lcp(0) {}

  STNode(size_type start, size_type end, size_type left, size_type right) :
    sp(start), ep(end),
    left_lcp(left), right_lcp(right)
  {
  }

  inline range_type range() const { return range_type(this->sp, this->ep); }
  inline size_type lcp() const { return std::max(this->left_lcp, this->right_lcp); }

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

  inline size_type size() const { return this->lcp_size; }
  inline size_type values() const { return this->data.size(); }
  inline size_type levels() const { return this->offsets.size() - 1; }
  inline size_type branching() const { return this->branching_factor; }

  inline size_type operator[] (size_type i) const { return this->data[i]; }

  // High-level interface.
  inline STNode root() const { return STNode(0, this->size() - 1, 0, 0); }
  STNode parent(const STNode& node) const;
  STNode parent(range_type range) const;

  // Low-level interface. The return value is a pair (res, LCP[res]) or notFound().
  range_type psv(size_type pos) const;
  range_type psev(size_type pos) const;
  range_type nsv(size_type pos) const;
  range_type nsev(size_type pos) const;

  inline STNode nodeFor(range_type range) const
  {
    if(range.second + 1 < this->size())
    {
      return STNode(range.first, range.second, this->data[range.first], this->data[range.second + 1]);
    }
    else
    {
      return STNode(range.first, range.second, this->data[range.first], 0);
    }
  }

  // Returned when a psv/nsv query cannot find a smaller value.
  inline range_type notFound() const { return range_type(this->values(), 0); }

//------------------------------------------------------------------------------

  /*
    We store a k-ary range minimum tree over the LCP array. Each node is identified by
    its position in the data array. Level i occupies range [offsets[i], offsets[i + 1] - 1]
    in the array. Level 0 is the leaves (the LCP array).
  */

  sdsl::int_vector<0>  data;
  sdsl::int_vector<64> offsets;
  size_type            lcp_size, branching_factor;

private:
  void copy(const LCPArray& source);

//------------------------------------------------------------------------------

  /*
    Tree operations on the range minimum tree. If level is required, it refers to the
    level of the parameter node.
  */

public:
  inline size_type rmtRoot() const
  {
    return this->values() - 1;
  }

  inline size_type rmtParent(size_type node, size_type level) const
  {
    return this->offsets[level + 1] + (node - this->offsets[level]) / this->branching();
  }

  inline bool rmtIsFirst(size_type node, size_type level) const
  {
    return ((node - this->offsets[level]) % this->branching() == 0);
  }

  inline size_type rmtFirstSibling(size_type node, size_type level) const
  {
    return node - (node - this->offsets[level]) % this->branching();
  }

  inline size_type rmtLastSibling(size_type first_child, size_type level) const
  {
    return std::min(this->offsets[level + 1], first_child + this->branching()) - 1;
  }

  inline size_type rmtFirstChild(size_type node, size_type level) const
  {
    return this->offsets[level - 1] + (node - this->offsets[level]) * this->branching();
  }

  inline size_type rmtLastChild(size_type node, size_type level) const
  {
    return this->rmtLastSibling(this->rmtFirstChild(node, level), level - 1);
  }
};  // class LCPArray

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_LCP_H
