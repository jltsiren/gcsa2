/*
  Copyright (c) 2019 Jouni Sirén
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

#include <gcsa/dbg.h>

namespace gcsa
{

//------------------------------------------------------------------------------

// Other class variables.

const std::string DeBruijnGraph::EXTENSION = ".dbg";

//------------------------------------------------------------------------------

DeBruijnGraph::DeBruijnGraph()
{
  this->node_count = 0;
  this->graph_order = 0;
}

DeBruijnGraph::DeBruijnGraph(const DeBruijnGraph& g)
{
  this->copy(g);
}

DeBruijnGraph::DeBruijnGraph(DeBruijnGraph&& g)
{
  *this = std::move(g);
}

DeBruijnGraph::~DeBruijnGraph()
{
}

void
DeBruijnGraph::copy(const DeBruijnGraph& g)
{
  this->node_count = g.node_count;
  this->graph_order = g.graph_order;

  this->alpha = g.alpha;

  this->bwt = g.bwt;
  this->bwt_rank = g.bwt_rank;

  this->nodes = g.nodes;
  this->node_rank = g.node_rank;

  this->setVectors();
}

void
DeBruijnGraph::swap(DeBruijnGraph& g)
{
  if(this != &g)
  {
    std::swap(this->node_count, g.node_count);
    std::swap(this->graph_order, g.graph_order);

    this->alpha.swap(g.alpha);

    this->bwt.swap(g.bwt);
    sdsl::util::swap_support(this->bwt_rank, g.bwt_rank, &(this->bwt), &(g.bwt));

    this->nodes.swap(g.nodes);
    sdsl::util::swap_support(this->node_rank, g.node_rank, &(this->nodes), &(g.nodes));
  }
}

DeBruijnGraph&
DeBruijnGraph::operator=(const DeBruijnGraph& g)
{
  if(this != &g) { this->copy(g); }
  return *this;
}

DeBruijnGraph&
DeBruijnGraph::operator=(DeBruijnGraph&& g)
{
  if(this != &g)
  {
    this->node_count = std::move(g.node_count);
    this->graph_order = std::move(g.graph_order);

    this->alpha = std::move(g.alpha);

    this->bwt = std::move(g.bwt);
    this->bwt_rank = std::move(g.bwt_rank);

    this->nodes = std::move(g.nodes);
    this->node_rank = std::move(g.node_rank);

    this->setVectors();
  }
  return *this;
}

void
DeBruijnGraph::setVectors()
{
  this->bwt_rank.set_vector(&(this->bwt));
  this->node_rank.set_vector(&(this->nodes));
}

DeBruijnGraph::size_type
DeBruijnGraph::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->node_count, out, child, "node_count");
  written_bytes += sdsl::write_member(this->graph_order, out, child, "graph_order");

  written_bytes += this->alpha.serialize(out, child, "alpha");

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->bwt_rank.serialize(out, child, "bwt_rank");

  written_bytes += this->nodes.serialize(out, child, "nodes");
  written_bytes += this->node_rank.serialize(out, child, "node_rank");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
DeBruijnGraph::load(std::istream& in)
{
  sdsl::read_member(this->node_count, in);
  sdsl::read_member(this->graph_order, in);

  this->alpha.load(in);

  this->bwt.load(in);
  this->bwt_rank.load(in, &(this->bwt));

  this->nodes.load(in);
  this->node_rank.load(in, &(this->nodes));
}

//------------------------------------------------------------------------------

DeBruijnGraph::DeBruijnGraph(const std::vector<key_type>& keys, size_type kmer_length, const Alphabet& alphabet)
{
  this->node_count = keys.size();
  this->graph_order = kmer_length;

  size_type total_edges = 0;
  for(size_type i = 0; i < keys.size(); i++) { total_edges += sdsl::bits::lt_cnt[Key::predecessors(keys[i])]; }

  sdsl::int_vector<64> counts(alphabet.sigma, 0);
  bit_vector bwt_buffer(alphabet.sigma * total_edges, 0);
  bit_vector node_buffer(total_edges, 0);
  for(size_type i = 0, edge_pos = 0; i < keys.size(); i++)
  {
    size_type pred = Key::predecessors(keys[i]);
    for(size_type j = 0; j < alphabet.sigma; j++)
    {
      if(pred & (((size_type)1) << j))
      {
        bwt_buffer[j * this->size() + i] = 1;
        counts[j]++;
      }
    }
    edge_pos += sdsl::bits::lt_cnt[Key::successors(keys[i])];
    node_buffer[edge_pos - 1] = 1;
  }
  this->alpha = Alphabet(counts, alphabet.char2comp, alphabet.comp2char);
  this->bwt = bwt_buffer; sdsl::util::clear(bwt_buffer);
  this->nodes = node_buffer; sdsl::util::clear(node_buffer);

  this->initSupport();
}

//------------------------------------------------------------------------------

void
DeBruijnGraph::initSupport()
{
  sdsl::util::init_support(this->bwt_rank, &(this->bwt));
  sdsl::util::init_support(this->node_rank, &(this->nodes));
}

//------------------------------------------------------------------------------

} // namespace gcsa
