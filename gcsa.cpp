#include "gcsa.h"

namespace gcsa
{

//------------------------------------------------------------------------------

GCSA::GCSA()
{
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

} // namespace gcsa
