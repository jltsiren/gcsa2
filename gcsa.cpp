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
  this->max_query_length = 0;
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
  this->max_query_length = g.max_query_length;

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
    std::swap(this->max_query_length, g.max_query_length);

    this->bwt.swap(g.bwt);
    this->alpha.swap(g.alpha);

    this->nodes.swap(g.nodes);
    sdsl::util::swap_support(this->node_rank, g.node_rank, &(this->nodes), &(g.nodes));
    sdsl::util::swap_support(this->node_select, g.node_select, &(this->nodes), &(g.nodes));

    this->edges.swap(g.edges);
    sdsl::util::swap_support(this->edge_rank, g.edge_rank, &(this->edges), &(g.edges));
    sdsl::util::swap_support(this->edge_select, g.edge_select, &(this->edges), &(g.edges));

    this->sampled_nodes.swap(g.sampled_nodes);
    sdsl::util::swap_support(this->sampled_node_rank, g.sampled_node_rank, &(this->sampled_nodes), &(g.sampled_nodes));

    this->stored_samples.swap(g.stored_samples);
    this->samples.swap(g.samples);
    sdsl::util::swap_support(this->sample_select, g.sample_select, &(this->samples), &(g.samples));
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
    this->max_query_length = std::move(g.max_query_length);

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
GCSA::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->node_count, out, child, "node_count");
  written_bytes += sdsl::write_member(this->max_query_length, out, child, "max_query_length");

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

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSA::load(std::istream& in)
{
  sdsl::read_member(this->node_count, in);
  sdsl::read_member(this->max_query_length, in);

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

GCSA::GCSA(const std::vector<key_type>& keys, size_type kmer_length, const Alphabet& _alpha)
{
  this->node_count = keys.size();
  this->max_query_length = kmer_length;

  size_type total_edges = 0, alt_edges = 0;
  for(size_type i = 0; i < keys.size(); i++) { total_edges += sdsl::bits::lt_cnt[Key::predecessors(keys[i])];
  alt_edges += sdsl::bits::lt_cnt[Key::successors(keys[i])]; }

  sdsl::int_vector<64> counts(_alpha.sigma, 0);
  sdsl::int_vector<8> buffer(total_edges, 0);
  this->nodes = bit_vector(total_edges, 0);
  this->edges = bit_vector(total_edges, 0);
  for(size_type i = 0, bwt_pos = 0, edge_pos = 0; i < keys.size(); i++)
  {
    byte_type pred = Key::predecessors(keys[i]);
    for(size_type j = 0; j < _alpha.sigma; j++)
    {
      if(pred & (((size_type)1) << j))
      {
        buffer[bwt_pos] = j; bwt_pos++;
        counts[j]++;
      }
    }
    this->nodes[bwt_pos - 1] = 1;
    edge_pos += sdsl::bits::lt_cnt[Key::successors(keys[i])];
    this->edges[edge_pos - 1] = 1;
  }
  directConstruct(this->bwt, buffer);
  this->alpha = Alphabet(counts, _alpha.char2comp, _alpha.comp2char);

  sdsl::util::init_support(this->node_rank, &(this->nodes));
  sdsl::util::init_support(this->node_select, &(this->nodes));
  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->edge_select, &(this->edges));
}
//------------------------------------------------------------------------------

GCSA::GCSA(std::vector<KMer>& kmers, size_type kmer_length, const Alphabet& _alpha)
{
  std::vector<key_type> keys;
  sdsl::int_vector<0> last_chars;
  uniqueKeys(kmers, keys, last_chars);
  GCSA mapper(keys, kmer_length, _alpha);

  // FIXME implement

  /*
    Transform the KMers into DoublingNodes.
    FIXME Later: save memory by not having both KMers and DoublingNodes  in memory

    A single doubling step (out of three):
    - sort the previous generation by from
    - build an index structure to find paths quickly by from value
    - scan the previous generation, build the next generation
      * next = DoublingNode(left, right) if left to == right.from
      * if left.sorted(), output it directly instead
    - delete the previous generation
    - sort the next generation by labels
    - merge sorted ranges of size > 1
      * a range is either a single path or multiple adjacent label values with all paths
        having the same from node
      * set to = from
    - stop if fully sorted
    - after all steps set max_query_length

    FIXME Later: parallelize, save memory by writing the next generation to disk
  */

  /*
    Merging and sampling:
      - scan the paths
        * sample if multiple predecessors, is source, or if offset == 0 for at least one node
      - set node count
      - set sampled_nodes, sampled_node_rank, stored_samples, samples, sample_rank

    FIXME Later: alternate sampling scheme not based on (id, offset) pairs. Sample a node,
    if it has multiple predecessors, it is a source node, or if
    { from(node) } != { from(predecessor(node)) } + 1. Also some samples if distance to
    the next sample is too large.
  */

  /*
    Edge generation:
    - set the to fields of the last generation to 0,
    - build an index structure to find paths quickly by labels
    - scan the paths in the last generation, writing the number of outgoing edges to the to field
      * use mapper, predecessor field in keys, and last_chars to find the predecessors
      * increment the to field of the predecessor

    FIXME Later: parallelize; load KMers from disk before proceeding
  */

  /*
    GCSA construction:
    - count character occurrences
    - set bwt, alpha
    - set nodes, node_rank, node_select, edges, edge_rank, edge_select
  */
}

//------------------------------------------------------------------------------

range_type
GCSA::find(const std::string& pattern) const
{
  return this->find(pattern.begin(), pattern.end());
}

range_type
GCSA::find(const char* pattern, size_type length) const
{
  return this->find(pattern, pattern + length);
}

void
GCSA::locate(size_type node, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { sdsl::util::clear(results); }
  if(node >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  this->locateInternal(node, results);
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locate(range_type range, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { sdsl::util::clear(results); }
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
GCSA::locateInternal(size_type node, std::vector<node_type>& results) const
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
