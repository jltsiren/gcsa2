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

#include <cstdio>
#include <cstdlib>

#include "gcsa.h"
#include "path_graph.h"

namespace gcsa
{

//------------------------------------------------------------------------------

const std::string GCSA::EXTENSION = ".gcsa";

GCSA::GCSA()
{
  this->path_node_count = 0;
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
  this->path_node_count = g.path_node_count;
  this->max_query_length = g.max_query_length;

  this->bwt = g.bwt;
  this->alpha = g.alpha;

  this->path_nodes = g.path_nodes;
  this->path_rank = g.path_rank;
  this->path_select = g.path_select;

  this->edges = g.edges;
  this->edge_rank = g.edge_rank;
  this->edge_select = g.edge_select;

  this->sampled_paths = g.sampled_paths;
  this->sampled_path_rank = g.sampled_path_rank;

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
    std::swap(this->path_node_count, g.path_node_count);
    std::swap(this->max_query_length, g.max_query_length);

    this->bwt.swap(g.bwt);
    this->alpha.swap(g.alpha);

    this->path_nodes.swap(g.path_nodes);
    sdsl::util::swap_support(this->path_rank, g.path_rank, &(this->path_nodes), &(g.path_nodes));
    sdsl::util::swap_support(this->path_select, g.path_select, &(this->path_nodes), &(g.path_nodes));

    this->edges.swap(g.edges);
    sdsl::util::swap_support(this->edge_rank, g.edge_rank, &(this->edges), &(g.edges));
    sdsl::util::swap_support(this->edge_select, g.edge_select, &(this->edges), &(g.edges));

    this->sampled_paths.swap(g.sampled_paths);
    sdsl::util::swap_support(this->sampled_path_rank, g.sampled_path_rank, &(this->sampled_paths), &(g.sampled_paths));

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
    this->path_node_count = std::move(g.path_node_count);
    this->max_query_length = std::move(g.max_query_length);

    this->bwt = std::move(g.bwt);
    this->alpha = std::move(g.alpha);

    this->path_nodes = std::move(g.path_nodes);
    this->path_rank = std::move(g.path_rank);
    this->path_select = std::move(g.path_select);

    this->edges = std::move(g.edges);
    this->edge_rank = std::move(g.edge_rank);
    this->edge_select = std::move(g.edge_select);

    this->sampled_paths = std::move(g.sampled_paths);
    this->sampled_path_rank = std::move(g.sampled_path_rank);

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
  this->path_rank.set_vector(&(this->path_nodes));
  this->path_select.set_vector(&(this->path_nodes));

  this->edge_rank.set_vector(&(this->edges));
  this->edge_select.set_vector(&(this->edges));

  this->sampled_path_rank.set_vector(&(this->sampled_paths));

  this->sample_select.set_vector(&(this->samples));
}

GCSA::size_type
GCSA::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->path_node_count, out, child, "path_node_count");
  written_bytes += sdsl::write_member(this->max_query_length, out, child, "max_query_length");

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->alpha.serialize(out, child, "alpha");

  written_bytes += this->path_nodes.serialize(out, child, "path_nodes");
  written_bytes += this->path_rank.serialize(out, child, "path_rank");
  written_bytes += this->path_select.serialize(out, child, "path_select");

  written_bytes += this->edges.serialize(out, child, "edges");
  written_bytes += this->edge_rank.serialize(out, child, "edge_rank");
  written_bytes += this->edge_select.serialize(out, child, "edge_select");

  written_bytes += this->sampled_paths.serialize(out, child, "sampled_paths");
  written_bytes += this->sampled_path_rank.serialize(out, child, "sampled_path_rank");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSA::load(std::istream& in)
{
  sdsl::read_member(this->path_node_count, in);
  sdsl::read_member(this->max_query_length, in);

  this->bwt.load(in);
  this->alpha.load(in);

  this->path_nodes.load(in);
  this->path_rank.load(in, &(this->path_nodes));
  this->path_select.load(in, &(this->path_nodes));

  this->edges.load(in);
  this->edge_rank.load(in, &(this->edges));
  this->edge_select.load(in, &(this->edges));

  this->sampled_paths.load(in);
  this->sampled_path_rank.load(in, &(this->sampled_paths));

  this->stored_samples.load(in);
  this->samples.load(in);
  this->sample_select.load(in, &(this->samples));
}

//------------------------------------------------------------------------------

void
predecessor(const PathNode& curr, const PathNode::rank_type* labels,
  comp_type comp, const DeBruijnGraph& mapper, const sdsl::int_vector<0>& last_char,
  PathLabel& first, PathLabel& last)
{
  size_type i = 0, j = curr.pointer();
  first.first = true; last.first = false;

  // Handle the common prefix of the labels.
  while(i < curr.lcp())
  {
    first.label[i] = last.label[i] = mapper.LF(labels[j], comp);
    comp = last_char[labels[j]];
    i++; j++;
  }

  // Handle the diverging suffixes of the labels.
  comp_type first_comp = comp, last_comp = comp;
  if(i < curr.order())
  {
    first.label[i] = mapper.LF(labels[j], first_comp);
    first_comp = last_char[labels[j]];
    last.label[i] = mapper.LF(labels[j + 1], last_comp);
    last_comp = last_char[labels[j + 1]];
    i++;
  }
  if(i < PathLabel::LABEL_LENGTH)
  {
    first.label[i] = mapper.node_rank(mapper.alpha.C[first_comp]);
    last.label[i] = mapper.node_rank(mapper.alpha.C[last_comp + 1]) - 1;
    i++;
  }

  // Pad the labels.
  while(i < PathLabel::LABEL_LENGTH)
  {
    first.label[i] = 0; last.label[i] = PathLabel::NO_RANK; i++;
  }
}

void
fromNodes(const MergedGraph& graph, size_type path, size_type& pointer, std::vector<node_type>& results, bool revert)
{
  results.clear();
  results.push_back(graph.paths[path].from);

  while(pointer < graph.size() && graph.from_nodes[pointer].first < path) { pointer++; }
  size_type old_pointer = pointer;
  while(pointer < graph.size() && graph.from_nodes[pointer].first == path)
  {
    results.push_back(graph.from_nodes[pointer].second); pointer++;
  }
  if(revert) { pointer = old_pointer; } // We may have the same predecessor for multiple nodes.

  removeDuplicates(results, false);
}

GCSA::GCSA(const InputGraph& graph,
  size_type doubling_steps, size_type size_limit, const Alphabet& _alpha)
{
  if(graph.size() == 0) { return; }
  if(doubling_steps > DOUBLING_STEPS)
  {
    std::cerr << "GCSA::GCSA(): The number of doubling steps is too high: " << doubling_steps << std::endl;
    std::cerr << "GCSA::GCSA(): Reverting the number of doubling steps to " << DOUBLING_STEPS << std::endl;
    doubling_steps = DOUBLING_STEPS;
  }
  if(size_limit > ABSOLUTE_SIZE_LIMIT)
  {
    std::cerr << "GCSA::GCSA(): The size limit is suspiciously high: " << size_limit << " GB" << std::endl;
    std::cerr << "GCSA::GCSA(): Reverting the size limit to " << ABSOLUTE_SIZE_LIMIT << " GB" << std::endl;
    size_limit = ABSOLUTE_SIZE_LIMIT;
  }
  size_type bytes_required = graph.size() * (sizeof(PathNode) + 2 * sizeof(PathNode::rank_type));
  if(bytes_required > size_limit * GIGABYTE)
  {
    std::cerr << "GCSA::GCSA(): The input is too large: " << (bytes_required / GIGABYTE_DOUBLE) << " GB" << std::endl;
    std::cerr << "GCSA::GCSA(): Construction aborted" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Extract the keys and build the necessary support structures.
  // FIXME Later: Write the structures to disk until needed?
  std::vector<key_type> keys;
  graph.read(keys);
  DeBruijnGraph mapper(keys, graph.k(), _alpha);
  LCP lcp(keys, graph.k());
  sdsl::int_vector<0> last_char;
  Key::lastChars(keys, last_char);
  sdsl::sd_vector_builder builder(Key::label(keys[keys.size() - 1]) + 1, keys.size());
  for(size_type i = 0; i < keys.size(); i++) { builder.set(Key::label(keys[i])); }
  sdsl::sd_vector<> key_exists(builder);
  sdsl::util::clear(keys);

  // Create the initial PathGraph.
  PathGraph path_graph(graph, key_exists);
  sdsl::util::clear(key_exists);

  // Prefix-doubling.
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::GCSA(): Initial path length: " << path_graph.k() << std::endl;
#endif
  path_graph.prune(lcp, size_limit * GIGABYTE);
  for(size_type step = 1; step <= doubling_steps && path_graph.unsorted > 0; step++)
  {
#ifdef VERBOSE_STATUS_INFO
    std::cerr << "GCSA::GCSA(): Step " << step << " (path length " << path_graph.k() << " -> "
              << (2 * path_graph.k()) << ")" << std::endl;
#endif
    path_graph.extend(size_limit * GIGABYTE);
    path_graph.prune(lcp, size_limit * GIGABYTE);
  }
  this->max_query_length = (path_graph.unsorted == 0 ? ~(size_type)0 : path_graph.k());
  sdsl::util::clear(lcp);

  // Merge the paths into the nodes of a pruned de Bruijn graph.
  MergedGraph merged_graph(path_graph, mapper);
  this->path_node_count = merged_graph.size();
  path_graph.clear();

  // Structures used to build GCSA.
  sdsl::int_vector<64> counts(mapper.alpha.sigma, 0); // alpha
  sdsl::int_vector<8> bwt_buffer(merged_graph.size() + merged_graph.size() / 2, 0); // bwt
  this->path_nodes = bit_vector(bwt_buffer.size(), 0);
  SLArray outdegrees(merged_graph.size());  // edges
  this->sampled_paths = bit_vector(merged_graph.size(), 0);
  std::vector<node_type> sample_buffer; // stored_samples
  this->samples = bit_vector(merged_graph.size() + merged_graph.extra(), 0);

  // The actual construction.
  // FIXME use MemoryMapManagers
  PathLabel first, last;
  size_type from_pointer = 0, total_edges = 0, sample_bits = 0;
  std::vector<node_type> pred_from, curr_from;
  for(size_type i = 0; i < merged_graph.size(); i++)
  {
    // Find the predecessors.
    size_type indegree = 0, pred_id = 0, pred_comp = 0;
    bool sample_this = false;
    for(size_type comp = 0; comp < mapper.alpha.sigma; comp++)
    {
      if(!(merged_graph.paths[i].hasPredecessor(comp))) { continue; }

      // Find the predecessor of paths[i] with comp and the first path intersecting it.
      predecessor(merged_graph.paths[i], merged_graph.labels, comp, mapper, last_char, first, last);
      while(!(merged_graph.paths[merged_graph.next[comp]].intersect(first, last, merged_graph.labels))) {
          merged_graph.next[comp]++;
      }

      // Add an edge.
      indegree++; outdegrees.increment(merged_graph.next[comp]); total_edges++;
      if(total_edges > bwt_buffer.size()) { bwt_buffer.resize(bwt_buffer.size() + merged_graph.size() / 2); }
      bwt_buffer[total_edges - 1] = comp; counts[comp]++;
      pred_id = merged_graph.next[comp]; pred_comp = comp; // For sampling.

      // Create additional edges if the next path nodes also intersect with the predecessor.
      while(merged_graph.next[comp] + 1 < merged_graph.size() &&
            merged_graph.paths[merged_graph.next[comp] + 1].intersect(first, last, merged_graph.labels))
      {
        merged_graph.next[comp]++;
        indegree++; outdegrees.increment(merged_graph.next[comp]); total_edges++;
        if(total_edges > bwt_buffer.size()) { bwt_buffer.resize(bwt_buffer.size() + merged_graph.size() / 2); }
        bwt_buffer[total_edges - 1] = comp; counts[comp]++;
      }
    }
    if(total_edges > this->path_nodes.size()) { this->path_nodes.resize(bwt_buffer.size()); }
    this->path_nodes[total_edges - 1] = 1;

    /*
      Simple cases for sampling the node:
      - multiple predecessors
      - at the beginning of the source node with no real predecessors
      - at the beginning of a node in the original graph (makes the previous case redundant)
    */
    fromNodes(merged_graph, i, from_pointer, curr_from, false);
    if(indegree > 1) { sample_this = true; }
    if(merged_graph.paths[i].hasPredecessor(Alphabet::SINK_COMP)) { sample_this = true; }
    for(size_type k = 0; k < curr_from.size(); k++)
    {
      if(Node::offset(curr_from[k]) == 0) { sample_this = true; break; }
    }

    // Sample if the from nodes cannot be derived from the only predecessor.
    if(!sample_this)
    {
      fromNodes(merged_graph, pred_id, merged_graph.next_from[pred_comp], pred_from, true);
      if(pred_from.size() != curr_from.size()) { sample_this = true; }
      else
      {
        for(size_type k = 0; k < curr_from.size(); k++)
        {
          if(curr_from[k] != pred_from[k] + 1) { sample_this = true; break; }
        }
      }
    }

    // Store the samples.
    if(sample_this)
    {
      this->sampled_paths[i] = 1;
      for(size_type k = 0; k < curr_from.size(); k++)
      {
        sample_bits = std::max(sample_bits, bit_length(curr_from[k]));
        sample_buffer.push_back(curr_from[k]);
      }
      this->samples[sample_buffer.size() - 1] = 1;
    }
  }
  sdsl::util::clear(last_char);

  // Initialize alpha.
  this->alpha = Alphabet(counts, mapper.alpha.char2comp, mapper.alpha.comp2char);
  sdsl::util::clear(mapper);

  // Initialize bwt.
  bwt_buffer.resize(total_edges);
  directConstruct(this->bwt, bwt_buffer);
  sdsl::util::clear(bwt_buffer);

  // Initialize bitvectors (path_nodes, edges, sampled_paths, samples).
  this->path_nodes.resize(total_edges);
  this->edges = bit_vector(total_edges, 0); total_edges = 0;
  for(size_type i = 0; i < merged_graph.size(); i++)
  {
    total_edges += outdegrees[i];
    this->edges[total_edges - 1] = 1;
  }
  outdegrees.clear();
  this->samples.resize(sample_buffer.size());
  this->initSupport();

  // Initialize stored_samples.
  this->stored_samples = sdsl::int_vector<0>(sample_buffer.size(), 0, sample_bits);
  for(size_type i = 0; i < sample_buffer.size(); i++) { this->stored_samples[i] = sample_buffer[i]; }
  sdsl::util::clear(sample_buffer);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::GCSA(): " << total_edges << " edges" << std::endl;
  std::cerr << "GCSA::GCSA(): " << this->sample_count() << " samples at "
            << this->sampled_positions() << " positions" << std::endl;
#endif
}

//------------------------------------------------------------------------------

std::ostream&
printOccs(const std::vector<node_type>& occs, std::ostream& out)
{
  out << "{";
  for(size_type i = 0; i < occs.size(); i++)
  {
    out << (i == 0 ? " " : ", ") << Node::decode(occs[i]);
  }
  out << " }";
  return out;
}

void
printFailure(const std::string& kmer,
             const std::vector<node_type>& expected, const std::vector<node_type>& occs)
{
  std::cerr << "GCSA::verifyIndex(): locate(" << kmer << ") failed" << std::endl;
  std::cerr << "GCSA::verifyIndex(): Expected ";
  printOccs(expected, std::cerr) << std::endl;
  std::cerr << "GCSA::verifyIndex(): Got ";
  printOccs(occs, std::cerr) << std::endl;
}

struct KMerSplitComparator
{
  inline bool operator() (const KMer& left, const KMer& right) const
    {
      return (Key::label(left.key) != Key::label(right.key));
    }
};

bool
GCSA::verifyIndex(std::vector<KMer>& kmers, size_type kmer_length) const
{
  size_type threads = omp_get_max_threads();
  parallelQuickSort(kmers.begin(), kmers.end());
  KMerSplitComparator k_comp;
  std::vector<range_type> bounds = getBounds(kmers, threads, k_comp);
  assert(bounds.size() == threads);

  size_type fails = 0;
#pragma omp parallel for schedule(static)
  for(size_type thread = 0; thread < threads; thread++)
  {
    size_type i = bounds[thread].first;
    while(i <= bounds[thread].second)
    {
      size_type next = i + 1;
      while(next <= bounds[thread].second && Key::label(kmers[next].key) == Key::label(kmers[i].key)) {
          next++;
      }

      std::string kmer = Key::decode(kmers[i].key, kmer_length, alpha);
      size_type endmarker_pos = kmer.find('$'); // The actual kmer ends at the first endmarker.
      if(endmarker_pos != std::string::npos) { kmer = kmer.substr(0, endmarker_pos + 1); }

      range_type range = this->find(kmer);
      if(Range::empty(range))
      {
        #pragma omp critical
        {
          std::cerr << "GCSA::verifyIndex(): find(" << kmer << ") returned empty range" << std::endl;
          fails++;
        }
        i = next; continue;
      }

      std::vector<node_type> expected;
      for(size_type j = i; j < next; j++) { expected.push_back(kmers[j].from); }
      removeDuplicates(expected, false);
      std::vector<node_type> occs;
      this->locate(range, occs);

      if(occs.size() != expected.size())
      {
        #pragma omp critical
        {
          std::cerr << "GCSA::verifyIndex(): Expected " << expected.size()
                    << " occurrences, got " << occs.size() << std::endl;
          printFailure(kmer, expected, occs); fails++;
        }
      }
      else
      {
        for(size_type j = 0; j < occs.size(); j++)
        {
          if(occs[j] != expected[j])
          {
            #pragma omp critical
            {
              std::cerr << "GCSA::verifyIndex(): Failure at " << j << ": "
                        << "expected " << Node::decode(expected[j])
                        << ", got " << Node::decode(occs[j]) << std::endl;
              printFailure(kmer, expected, occs); fails++;
            }
            break;
          }
        }
      }

      i = next;
    }
  }

  if(fails == 0)
  {
    std::cout << "Index verification complete." << std::endl;
  }
  else
  {
    std::cout << "Index verification failed for " << fails << " patterns." << std::endl;
  }
  std::cout << std::endl;

  return fails == 0;
}

//------------------------------------------------------------------------------

void
GCSA::initSupport()
{
  sdsl::util::init_support(this->path_rank, &(this->path_nodes));
  sdsl::util::init_support(this->path_select, &(this->path_nodes));
  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->edge_select, &(this->edges));
  sdsl::util::init_support(this->sampled_path_rank, &(this->sampled_paths));
  sdsl::util::init_support(this->sample_select, &(this->samples));
}

//------------------------------------------------------------------------------

void
GCSA::locate(size_type path_node, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { sdsl::util::clear(results); }
  if(path_node >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  this->locateInternal(path_node, results);
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locate(range_type range, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { sdsl::util::clear(results); }
  if(Range::empty(range) || range.second >= this->size())
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
GCSA::locateInternal(size_type path_node, std::vector<node_type>& results) const
{
  size_type steps = 0;
  while(this->sampled_paths[path_node] == 0)
  {
    path_node = this->LF(path_node);
    steps++;
  }

  range_type sample_range = this->sampleRange(path_node);
  for(size_type i = sample_range.first; i <= sample_range.second; i++)
  {
    results.push_back(this->stored_samples[i] + steps);
  }
}

//------------------------------------------------------------------------------

} // namespace gcsa
