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

#include "gcsa.h"

namespace gcsa
{

//------------------------------------------------------------------------------

const std::string GCSA::EXTENSION = ".gcsa";

GCSA::GCSA()
{
  this->node_count = 0;
}

GCSA::GCSA(const GCSA& g)
{
  this->copy(g);
}

GCSA::GCSA(GCSA&& g)
{
  *this = std::move(g);
}

GCSA::~GCSA()
{
}

void
GCSA::copy(const GCSA& g)
{
  this->node_count = g.node_count;

  this->bwt = g.bwt;
  this->alpha = g.alpha;

  this->nodes = g.nodes;
  this->node_rank = g.node_rank;
  this->node_select = g.node_select;

  this->edges = g.edges;
  this->edge_rank = g.edge_rank;
  this->edge_select = g.edge_select;

  this->sampled_nodes = g.sampled_nodes;
  this->sampled_node_rank = g.sampled_node_rank;

  this->stored_samples = g.stored_samples;
  this->samples = g.samples;
  this->sample_select = g.sample_select;

  this->setVectors();
}

void
GCSA::swap(GCSA& g)
{
  if(this != &g)
  {
    std::swap(this->node_count, g.node_count);

    this->bwt.swap(g.bwt);
    this->alpha.swap(g.alpha);

    this->nodes.swap(g.nodes);
    util::swap_support(this->node_rank, g.node_rank, &(this->nodes), &(g.nodes));
    util::swap_support(this->node_select, g.node_select, &(this->nodes), &(g.nodes));

    this->edges.swap(g.edges);
    util::swap_support(this->edge_rank, g.edge_rank, &(this->edges), &(g.edges));
    util::swap_support(this->edge_select, g.edge_select, &(this->edges), &(g.edges));

    this->sampled_nodes.swap(g.sampled_nodes);
    util::swap_support(this->sampled_node_rank, g.sampled_node_rank, &(this->sampled_nodes), &(g.sampled_nodes));

    this->stored_samples.swap(g.stored_samples);
    this->samples.swap(g.samples);
    util::swap_support(this->sample_select, g.sample_select, &(this->samples), &(g.samples));
  }
}

GCSA&
GCSA::operator=(const GCSA& g)
{
  if(this != &g) { this->copy(g); }
  return *this;
}

GCSA&
GCSA::operator=(GCSA&& g)
{
  if(this != &g)
  {
    this->node_count = std::move(g.node_count);

    this->bwt = std::move(g.bwt);
    this->alpha = std::move(g.alpha);

    this->nodes = std::move(g.nodes);
    this->node_rank = std::move(g.node_rank);
    this->node_select = std::move(g.node_select);

    this->edges = std::move(g.edges);
    this->edge_rank = std::move(g.edge_rank);
    this->edge_select = std::move(g.edge_select);

    this->sampled_nodes = std::move(g.sampled_nodes);
    this->sampled_node_rank = std::move(g.sampled_node_rank);

    this->stored_samples = std::move(g.stored_samples);
    this->samples = std::move(g.samples);
    this->sample_select = std::move(g.sample_select);

    this->setVectors();
  }
  return *this;
}

void
GCSA::setVectors()
{
  this->node_rank.set_vector(&(this->nodes));
  this->node_select.set_vector(&(this->nodes));

  this->edge_rank.set_vector(&(this->edges));
  this->edge_select.set_vector(&(this->edges));

  this->sampled_node_rank.set_vector(&(this->sampled_nodes));

  this->sample_select.set_vector(&(this->samples));
}

GCSA::size_type
GCSA::serialize(std::ostream& out, structure_tree_node* s, std::string name) const
{
  structure_tree_node* child = structure_tree::add_child(s, name, util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += write_member(this->node_count, out, child, "node_count");

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->alpha.serialize(out, child, "alpha");

  written_bytes += this->nodes.serialize(out, child, "nodes");
  written_bytes += this->node_rank.serialize(out, child, "node_rank");
  written_bytes += this->node_select.serialize(out, child, "node_select");

  written_bytes += this->edges.serialize(out, child, "edges");
  written_bytes += this->edge_rank.serialize(out, child, "edge_rank");
  written_bytes += this->edge_select.serialize(out, child, "edge_select");

  written_bytes += this->sampled_nodes.serialize(out, child, "sampled_nodes");
  written_bytes += this->sampled_node_rank.serialize(out, child, "sampled_node_rank");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSA::load(std::istream& in)
{
  read_member(this->node_count, in);

  this->bwt.load(in);
  this->alpha.load(in);

  this->nodes.load(in);
  this->node_rank.load(in, &(this->nodes));
  this->node_select.load(in, &(this->nodes));

  this->edges.load(in);
  this->edge_rank.load(in, &(this->edges));
  this->edge_select.load(in, &(this->edges));

  this->sampled_nodes.load(in);
  this->sampled_node_rank.load(in, &(this->sampled_nodes));

  this->stored_samples.load(in);
  this->samples.load(in);
  this->sample_select.load(in, &(this->samples));
}

//------------------------------------------------------------------------------

GCSA::GCSA(const std::vector<uint64_t>& kmers, const Alphabet& _alpha)
{
  this->node_count = kmers.size();

  size_type total_edges = 0;
  for(size_type i = 0; i < kmers.size(); i++) { total_edges += bits::lt_cnt[successors(kmers[i])]; }

  int_vector<64> counts(_alpha.sigma, 0);
  int_vector<8> buffer(total_edges, 0);
  util::assign(this->nodes, bit_vector(total_edges, 0));
  util::assign(this->edges, bit_vector(total_edges, 0));
  for(size_type i = 0, bwt_pos = 0, edge_pos = 0; i < kmers.size(); i++)
  {
    uint8_t pred = predecessors(kmers[i]);
    for(size_type j = 0; j < _alpha.sigma; j++)
    {
      if(pred & (((size_type)1) << j))
      {
        buffer[bwt_pos] = j; bwt_pos++;
        counts[j]++;
      }
    }
    this->nodes[bwt_pos - 1] = 1;
    edge_pos += bits::lt_cnt[successors(kmers[i])];
    this->edges[edge_pos - 1] = 1;
  }
  directConstruct(this->bwt, buffer);
  this->alpha = Alphabet(counts, _alpha.char2comp, _alpha.comp2char);

  util::init_support(this->node_rank, &(this->nodes));
  util::init_support(this->node_select, &(this->nodes));
  util::init_support(this->edge_rank, &(this->edges));
  util::init_support(this->edge_select, &(this->edges));
}

//------------------------------------------------------------------------------

range_type
GCSA::find(const std::string& pattern) const
{
  if(pattern.length() == 0) { return range_type(0, this->size() - 1); }

  auto iter = pattern.rbegin();
  range_type range = this->nodeRange(gcsa::charRange(this->alpha, this->alpha.char2comp[*iter]));
  --iter;

  while(!isEmpty(range) && iter != pattern.rend())
  {
    range = this->LF(range, *iter);
    --iter;
  }

  return range;
}

void
GCSA::locate(size_type node, std::vector<size_type>& results, bool append, bool sort) const
{
  if(!append) { util::clear(results); }
  if(node >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  this->locateInternal(node, results);
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locate(range_type range, std::vector<size_type>& results, bool append, bool sort) const
{
  if(!append) { util::clear(results); }
  if(isEmpty(range) || range.second >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  for(size_type i = range.first; i <= range.second; i++)
  {
    this->locateInternal(i, results);
  }
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locateInternal(size_type node, std::vector<size_type>& results) const
{
  size_type steps = 0;
  while(this->sampled_nodes[node] == 0)
  {
    node = this->LF(node);
    steps++;
  }

  range_type sample_range = this->sampleRange(node);
  for(size_type i = sample_range.first; i <= sample_range.second; i++)
  {
    results.push_back(this->stored_samples[i] + steps);
  }
}

//------------------------------------------------------------------------------

} // namespace gcsa
