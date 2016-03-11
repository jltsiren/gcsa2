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

#include "internal.h"
#include "lcp.h"

namespace gcsa
{
//------------------------------------------------------------------------------

std::ostream&
operator<< (std::ostream& out, const STNode& node)
{
  out << "(" << node.left_lcp << ", " << node.range() << ", " << node.right_lcp << ")";
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

LCPArray::LCPArray(const InputGraph& graph, const ConstructionParameters& parameters) :
  branching_factor(parameters.lcp_branching)
{
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
  this->data = sdsl::int_vector<8>(total_size, ~(uint8_t)0);
  DiskIO::read(in, (const uint8_t*)(this->data.data()), this->lcp_size);
  in.close();
  for(size_type level = 0; level < this->levels(); level++)
  {
    for(size_type i = this->offsets[level]; i < this->offsets[level + 1]; i++)
    {
      size_type par = this->rmtParent(i, level);
      if(this->data[i] < this->data[par]) { this->data[par] = this->data[i]; }
    }
  }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "LCPArray::LCPArray(): " << this->values() << " values at " << this->levels()
            << " levels (branching factor " << this->branching() << ")" << std::endl;
#endif
}

//------------------------------------------------------------------------------

STNode
LCPArray::parent(const STNode& node) const
{
  if(node == this->root()) { return this->root(); }

  range_type left(node.sp, node.left_lcp), right(node.ep + 1, node.right_lcp);
  if(node.left_lcp >= node.right_lcp)
  {
    left = this->psv(node.sp);
    if(left == this->notFound()) { left.first = 0; }
  }
  if(node.right_lcp >= node.left_lcp)
  {
    right = this->nsv(node.ep + 1);
    if(right == this->notFound()) { right.first = this->size(); }
  }

  return STNode(left.first, right.first - 1, left.second, right.second);
}

STNode
LCPArray::parent(range_type range) const
{
  return this->parent(this->nodeFor(range));
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
  while(to != lcp.rmtRoot())
  {
    res = gcsa::psv(lcp, lcp.rmtFirstSibling(to, level), to, val, comp);
    if(res.first < lcp.values()) { break; }
    to = lcp.rmtParent(to, level); level++;
  }
  if(res.first >= lcp.values()) { return res; } // Not found.

  // Go to the leaf containing psv(to).
  while(level > 0)
  {
    size_type from = lcp.rmtFirstChild(res.first, level); level--;
    res = gcsa::psv(lcp, from, lcp.rmtLastSibling(from, level) + 1, val, comp);
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
  while(from != lcp.rmtRoot())
  {
    res = gcsa::nsv(lcp, from + 1, lcp.rmtLastSibling(from, level), val, comp);
    if(res.first < lcp.values()) { break; }
    from = lcp.rmtParent(from, level); level++;
  }
  if(res.first >= lcp.values()) { return res; }

  // Go to the leaf containing nsv(to).
  while(level > 0)
  {
    from = lcp.rmtFirstChild(res.first, level); level--;
    res = gcsa::nsv(lcp, from, lcp.rmtLastSibling(from, level), val, comp);
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

} // namespace gcsa
