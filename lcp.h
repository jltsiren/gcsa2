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

/*
  This is a specialization of GCSA for de Bruijn graphs. Because the graph is
  reverse deterministic, we can use indicator bitvectors for encoding the BWT.
  This simplifies LF() to two rank() operations.
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
  inline size_type levels() const { return this->offsets.size() - 1; }

  // FIXME queries

//------------------------------------------------------------------------------

  /*
    We store a k-ary range minimum tree over the LCP array. Each node is identified by
    itse position in the data array. Level i occupies range [offsets[i], offsets[i + 1] - 1]
    in the array. Level 0 is the leaves (the LCP array).
  */

  sdsl::int_vector<8>  data;
  sdsl::int_vector<64> offsets;
  size_type            lcp_size, branching_factor;

private:
  void copy(const LCPArray& source);

//------------------------------------------------------------------------------

  /*
    Tree operations on the range minimum tree. If level is required, it refers to the
    level of the parameter node.
  */

  inline size_type rmtRoot() const
  {
    return this->data.size() - 1;
  }

  inline size_type rmtLevel(size_type node) const
  {
    size_type level = 0;
    while(this->offsets[level + 1] <= node) { level++; }
    return level;
  }

  inline size_type rmtParent(size_type node, size_type level) const
  {
    return this->offsets[level + 1] + (node - this->offsets[level]) / this->branching_factor;
  }

  inline size_type rmtFirstSibling(size_type last_child, size_type level) const
  {
    return last_child - (last_child - this->offsets[level]) % this->branching_factor;
  }

  inline size_type rmtLastSibling(size_type first_child, size_type level) const
  {
    return std::min(this->offsets[level + 1] - 1, first_child + this->branching_factor - 1);
  }

  inline size_type rmtFirstChild(size_type node, size_type level) const
  {
    return this->offsets[level - 1] + (node - this->offsets[level]) * this->branching_factor;
  }

  inline size_type rmtLastChild(size_type node, size_type level) const
  {
    return this->rmtLastSibling(this->rmtFirstChild(node, level), level - 1);
  }

};  // class LCPArray

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_LCP_H
