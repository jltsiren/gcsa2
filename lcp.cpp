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

#include <stack>

#include "internal.h"
#include "lcp.h"

namespace gcsa
{
//------------------------------------------------------------------------------

std::ostream&
operator<< (std::ostream& out, const STNode& node)
{
  out << "(" << node.left_lcp << ", " << node.range() << " at depth " << node.lcp()
      << ", " << node.right_lcp << ")";
  return out;
}

//------------------------------------------------------------------------------

const std::string LCPArray::EXTENSION = ".lcp";

LCPArray::LCPArray() :
  lcp_size(0), branching_factor(2)
{
}

LCPArray::LCPArray(const LCPArray& source)
{
  this->copy(source);
}

LCPArray::LCPArray(LCPArray&& source)
{
  *this = std::move(source);
}

LCPArray::~LCPArray()
{
}

void
LCPArray::copy(const LCPArray& source)
{
  this->data = source.data;
  this->offsets = source.offsets;

  this->lcp_size = source.lcp_size;
  this->branching_factor = source.branching_factor;
}

void
LCPArray::swap(LCPArray& another)
{
  if(this != &another)
  {
    this->data.swap(another.data);
    this->offsets.swap(another.offsets);

    std::swap(this->lcp_size, another.lcp_size);
    std::swap(this->branching_factor, another.branching_factor);
  }
}

LCPArray&
LCPArray::operator=(const LCPArray& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

LCPArray&
LCPArray::operator=(LCPArray&& source)
{
  if(this != &source)
  {
    this->data = std::move(source.data);
    this->offsets = std::move(source.offsets);

    this->lcp_size = std::move(source.lcp_size);
    this->branching_factor = std::move(source.branching_factor);
  }
  return *this;
}

LCPArray::size_type
LCPArray::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->data.serialize(out, child, "data");
  written_bytes += this->offsets.serialize(out, child, "offsets");

  written_bytes += sdsl::write_member(this->lcp_size, out, child, "lcp_size");
  written_bytes += sdsl::write_member(this->branching_factor, out, child, "branching_factor");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
LCPArray::load(std::istream& in)
{
  this->data.load(in);
  this->offsets.load(in);

  sdsl::read_member(this->lcp_size, in);
  sdsl::read_member(this->branching_factor, in);
}

//------------------------------------------------------------------------------

/*
  Tree operations on the range minimum tree. If level is required, it refers to the
  level of the parameter node.
*/

inline size_type
rmtRoot(const LCPArray& lcp)
{
  return lcp.values() - 1;
}

inline size_type
rmtParent(const LCPArray& lcp, size_type node, size_type level)
{
  return lcp.offsets[level + 1] + (node - lcp.offsets[level]) / lcp.branching();
}

inline bool
rmtIsFirst(const LCPArray& lcp, size_type node, size_type level)
{
  return ((node - lcp.offsets[level]) % lcp.branching() == 0);
}

inline size_type
rmtFirstSibling(const LCPArray& lcp, size_type node, size_type level)
{
  return node - (node - lcp.offsets[level]) % lcp.branching();
}

inline size_type
rmtLastSibling(const LCPArray& lcp, size_type first_child, size_type level)
{
  return std::min(lcp.offsets[level + 1], first_child + lcp.branching()) - 1;
}

inline size_type
rmtFirstChild(const LCPArray& lcp, size_type node, size_type level)
{
  return lcp.offsets[level - 1] + (node - lcp.offsets[level]) * lcp.branching();
}

inline size_type
rmtLastChild(const LCPArray& lcp, size_type node, size_type level)
{
  return rmtLastSibling(lcp, rmtFirstChild(lcp, node, level), level - 1);
}

inline size_type
rmtLevel(const LCPArray& lcp, size_type node)
{
  size_type level = 0;
  while(lcp.offsets[level + 1] <= node) { level++; }
  return level;
}

//------------------------------------------------------------------------------

LCPArray::LCPArray(const InputGraph& graph, const ConstructionParameters& parameters) :
  branching_factor(parameters.lcp_branching)
{
  double start = readTimer();

  if(graph.lcp_name.empty())
  {
    std::cerr << "LCPArray::LCPArray: The input graph does not contain the LCP file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::ifstream in(graph.lcp_name.c_str(), std::ios_base::binary);
  if(!in)
  {
    std::cerr << "LCPArray::LCPArray: Cannot open LCP file " << graph.lcp_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Determine the number of levels.
  this->lcp_size = fileSize(in);
  size_type level_count = 1, level_size = this->lcp_size;
  while(level_size > 1)
  {
    level_count++; level_size = (level_size + this->branching() - 1) / this->branching();
  }

  // Initialize offsets.
  this->offsets = sdsl::int_vector<64>(level_count + 1, 0);
  level_size = this->lcp_size;
  size_type total_size = 0;
  for(size_type level = 0; level < this->levels(); level++)
  {
    total_size += level_size;
    this->offsets[level + 1] = total_size;
    level_size = (level_size + this->branching() - 1) / this->branching();
  }

  // Initialize data.
  this->data = sdsl::int_vector<0>(total_size, ~(uint8_t)0, 8);
  DiskIO::read(in, (const uint8_t*)(this->data.data()), this->lcp_size);
  in.close();
  for(size_type level = 0; level < this->levels(); level++)
  {
    for(size_type i = this->offsets[level]; i < this->offsets[level + 1]; i++)
    {
      size_type par = rmtParent(*this, i, level);
      if(this->data[i] < this->data[par]) { this->data[par] = this->data[i]; }
    }
  }
  sdsl::util::bit_compress(this->data);

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double seconds = readTimer() - start;
    std::cerr << "LCPArray::LCPArray(): Construction: " << seconds << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
  }
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "LCPArray::LCPArray(): " << this->values() << " values at " << this->levels()
              << " levels (branching factor " << this->branching() << ")" << std::endl;
  }
}

//------------------------------------------------------------------------------

LCPArray::node_type
LCPArray::parent(const LCPArray::node_type& node) const
{
  if(node == this->root()) { return this->root(); }

  size_type node_lcp = std::max(node.left_lcp, node.right_lcp);
  range_type left(node.sp, node.left_lcp), right(node.ep + 1, node.right_lcp);
  if(node.left_lcp == node_lcp)
  {
    left = this->psv(node.sp);
    if(left == this->notFound()) { left = range_type(0, 0); }
  }
  if(node.right_lcp == node_lcp)
  {
    right = this->nsv(node.ep + 1);
    if(right == this->notFound()) { right = range_type(this->size(), 0); }
  }

  return node_type(left.first, right.first - 1, left.second, right.second, node_lcp);
}

LCPArray::node_type
LCPArray::parent(range_type range) const
{
  return this->parent(this->nodeFor(range));
}

//------------------------------------------------------------------------------

size_type
LCPArray::depth(const LCPArray::node_type& node) const
{
  if(node.lcp() != node_type::UNKNOWN) { return node.lcp(); }
  return this->depth(node.range());
}

size_type
LCPArray::depth(LCPArray::node_type& node) const
{
  if(node.lcp() == node_type::UNKNOWN) { node.node_lcp = this->depth(node.range()); }
  return node.lcp();
}

size_type
LCPArray::depth(range_type range) const
{
  if(Range::length(range) <= 1) { return node_type::UNKNOWN; }
  range_type res = this->rmq(range.first + 1, range.second);
  return (res == this->notFound() ? node_type::UNKNOWN : res.second);
}

//------------------------------------------------------------------------------

/*
  Find the last value less than 'val' between 'from' (inclusive) and 'to' (exclusive).
  Return value will be lcp.notFound() if no such value exists.
*/
template<class Comparator>
range_type
psv(const LCPArray& lcp, size_type from, size_type to, size_type val, const Comparator& comp)
{
  while(to > from)
  {
    to--;
    if(comp(lcp[to], val)) { return range_type(to, lcp[to]); }
  }
  return lcp.notFound();
}

template<class Comparator>
range_type
psv(const LCPArray& lcp, size_type to, const Comparator& comp)
{
  if(to == 0 || to >= lcp.size()) { return lcp.notFound(); }

  // Find the children of the lowest common ancestor of psv(to) and 'to'.
  size_type level = 0, val = lcp[to];
  range_type res = lcp.notFound();
  while(to != rmtRoot(lcp))
  {
    res = gcsa::psv(lcp, rmtFirstSibling(lcp, to, level), to, val, comp);
    if(res.first < lcp.values()) { break; }
    to = rmtParent(lcp, to, level); level++;
  }
  if(res.first >= lcp.values()) { return res; } // Not found.

  // Go to the leaf containing psv(to).
  while(level > 0)
  {
    size_type from = rmtFirstChild(lcp, res.first, level); level--;
    res = gcsa::psv(lcp, from, rmtLastSibling(lcp, from, level) + 1, val, comp);
  }

  return res;
}

range_type
LCPArray::psv(size_type pos) const
{
  return gcsa::psv(*this, pos, std::less<size_type>());
}

range_type
LCPArray::psev(size_type pos) const
{
  return gcsa::psv(*this, pos, std::less_equal<size_type>());
}

//------------------------------------------------------------------------------

/*
  Find the first value less than 'val' between 'from' and 'to' (inclusive).
  Return value will be lcp.notFound() if no such value exists.
*/
template<class Comparator>
range_type
nsv(const LCPArray& lcp, size_type from, size_type to, size_type val, const Comparator& comp)
{
  for(size_type i = from; i <= to; i++)
  {
    if(comp(lcp[i], val)) { return range_type(i, lcp[i]); }
  }
  return lcp.notFound();
}

template<class Comparator>
range_type
nsv(const LCPArray& lcp, size_type from, const Comparator& comp)
{
  if(from + 1 >= lcp.size()) { return lcp.notFound(); }

  // Find the children of the lowest common ancestor for 'from' and nsv(from).
  size_type level = 0, val = lcp[from];
  range_type res = lcp.notFound();
  while(from != rmtRoot(lcp))
  {
    res = gcsa::nsv(lcp, from + 1, rmtLastSibling(lcp, from, level), val, comp);
    if(res.first < lcp.values()) { break; }
    from = rmtParent(lcp, from, level); level++;
  }
  if(res.first >= lcp.values()) { return res; }

  // Go to the leaf containing nsv(to).
  while(level > 0)
  {
    from = rmtFirstChild(lcp, res.first, level); level--;
    res = gcsa::nsv(lcp, from, rmtLastSibling(lcp, from, level), val, comp);
  }

  return res;
}

range_type
LCPArray::nsv(size_type pos) const
{
  return gcsa::nsv(*this, pos, std::less<size_type>());
}

range_type
LCPArray::nsev(size_type pos) const
{
  return gcsa::nsv(*this, pos, std::less_equal<size_type>());
}

//------------------------------------------------------------------------------

inline void
updateRes(const LCPArray& lcp, range_type& res, size_type i)
{
  if(lcp.data[i] < res.second) { res.first = i; res.second = lcp.data[i]; }
}

range_type
LCPArray::rmq(size_type sp, size_type ep) const
{
  if(sp > ep || ep >= this->size()) { return this->notFound(); }
  if(sp == ep) { return range_type(sp, this->data[sp]); }

  /*
    Search for a subtree containing the rmq, maintaining the following invariants:
      - left < right
      - nodes before subtree(left) are processed
      - nodes after subtree(right) are in tail
      - res contains (i, lcp[i]) for the rmq in the processed range
  */
  range_type res(this->values(), this->size());
  size_type level = 0, left = sp, right = ep;
  std::stack<range_type> tail;
  while(true)
  {
    size_type left_par = rmtParent(*this, left, level), right_par = rmtParent(*this, right, level);
    if(left_par == right_par)
    {
      for(size_type i = left; i <= right; i++) { updateRes(*this, res, i); }
      break;
    }

    size_type left_child = rmtFirstChild(*this, left_par, level + 1);
    if(left != left_child)
    {
      size_type last_child = rmtLastSibling(*this, left_child, level);
      for(size_type i = left; i <= last_child; i++) { updateRes(*this, res, i); }
      left_par++;
    }

    size_type right_child = rmtLastChild(*this, right_par, level + 1);
    if(right != right_child)
    {
      size_type first_child = rmtFirstSibling(*this, right_child, level);
      for(size_type i = right; i >= first_child; i--) { tail.push(range_type(i, this->data[i])); }
      right_par--;
    }

    if(left_par >= right_par)
    {
      if(left_par == right_par) { updateRes(*this, res, left_par); }
      break;
    }
    left = left_par; right = right_par; level++;
  }

  // Check the tail.
  while(!(tail.empty()))
  {
    range_type temp = tail.top(); tail.pop();
    if(temp.second < res.second) { res = temp; }
  }

  // Find the leftmost leaf in subtree(res.first) containing LCP value res.second.
  level = rmtLevel(*this, res.first);
  while(level > 0)
  {
    res.first = rmtFirstChild(*this, res.first, level); level--;
    while(this->data[res.first] != res.second) { res.first++; }
  }

  return res;
}

range_type
LCPArray::rmq(range_type range) const
{
  return this->rmq(range.first, range.second);
}

//------------------------------------------------------------------------------

} // namespace gcsa
