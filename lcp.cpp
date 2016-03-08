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
    level_count++; level_size = (level_size + this->branching_factor - 1) / this->branching_factor;
  }

  // Initialize offsets.
  this->offsets = sdsl::int_vector<64>(level_count + 1, 0);
  level_size = this->lcp_size;
  size_type total_size = this->lcp_size;
  for(size_type level = 0; level < this->levels(); level++)
  {
    this->offsets[level + 1] = total_size;
    level_size = (level_size + this->branching_factor - 1) / this->branching_factor;
    total_size += level_size;
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
}

//------------------------------------------------------------------------------

} // namespace gcsa
