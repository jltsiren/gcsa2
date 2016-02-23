/*
  Copyright (c) 2015, 2016 Genome Research Ltd.

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
#include <stack>

#include "gcsa.h"
#include "internal.h"
#include "path_graph.h"

namespace gcsa
{

//------------------------------------------------------------------------------

ConstructionParameters::ConstructionParameters() :
  doubling_steps(DOUBLING_STEPS), size_limit(SIZE_LIMIT * GIGABYTE)
{
}

void
ConstructionParameters::setSteps(size_type steps)
{
  this->doubling_steps = Range::bound(steps, 1, DOUBLING_STEPS);
}

void
ConstructionParameters::setLimit(size_type gigabytes)
{
  this->size_limit = Range::bound(gigabytes, 1, ABSOLUTE_LIMIT) * GIGABYTE;
}

//------------------------------------------------------------------------------

const std::string GCSA::EXTENSION = ".gcsa";

GCSA::GCSA()
{
}

GCSA::GCSA(const GCSA& source)
{
  this->copy(source);
}

GCSA::GCSA(GCSA&& source)
{
  *this = std::move(source);
}

GCSA::~GCSA()
{
}

void
GCSA::swap(GCSA& another)
{
  if(this != &another)
  {
    this->header.swap(another.header);

    this->bwt.swap(another.bwt);
    this->alpha.swap(another.alpha);

    this->path_nodes.swap(another.path_nodes);
    sdsl::util::swap_support(this->path_rank, another.path_rank, &(this->path_nodes), &(another.path_nodes));
    sdsl::util::swap_support(this->path_select, another.path_select, &(this->path_nodes), &(another.path_nodes));

    this->edges.swap(another.edges);
    sdsl::util::swap_support(this->edge_rank, another.edge_rank, &(this->edges), &(another.edges));
    sdsl::util::swap_support(this->edge_select, another.edge_select, &(this->edges), &(another.edges));

    this->sampled_paths.swap(another.sampled_paths);
    sdsl::util::swap_support(this->sampled_path_rank, another.sampled_path_rank, &(this->sampled_paths), &(another.sampled_paths));

    this->stored_samples.swap(another.stored_samples);
    this->samples.swap(another.samples);
    sdsl::util::swap_support(this->sample_select, another.sample_select, &(this->samples), &(another.samples));

    this->counter.swap(another.counter);
    this->total_pointers.swap(another.total_pointers);
    this->redundant_pointers.swap(another.redundant_pointers);
  }
}

GCSA&
GCSA::operator=(const GCSA& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

GCSA&
GCSA::operator=(GCSA&& source)
{
  if(this != &source)
  {
    this->header = std::move(source.header);

    this->bwt = std::move(source.bwt);
    this->alpha = std::move(source.alpha);

    this->path_nodes = std::move(source.path_nodes);
    this->path_rank = std::move(source.path_rank);
    this->path_select = std::move(source.path_select);

    this->edges = std::move(source.edges);
    this->edge_rank = std::move(source.edge_rank);
    this->edge_select = std::move(source.edge_select);

    this->sampled_paths = std::move(source.sampled_paths);
    this->sampled_path_rank = std::move(source.sampled_path_rank);

    this->stored_samples = std::move(source.stored_samples);
    this->samples = std::move(source.samples);
    this->sample_select = std::move(source.sample_select);

    this->counter = std::move(source.counter);
    this->total_pointers = std::move(source.total_pointers);
    this->redundant_pointers = std::move(source.redundant_pointers);

    this->setVectors();
  }
  return *this;
}

GCSA::size_type
GCSA::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");

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

  written_bytes += this->counter.serialize(out, child, "counter");
  written_bytes += this->total_pointers.serialize(out, child, "total_pointers");
  written_bytes += this->redundant_pointers.serialize(out, child, "redundant_pointers");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSA::load(std::istream& in)
{
  this->header.load(in);
  if(!(this->header.check()))
  {
    std::cerr << "GCSA::load(): Invalid header: " << this->header << std::endl;
  }

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

  this->counter.load(in);
  this->total_pointers.load(in);
  this->redundant_pointers.load(in);
}

//------------------------------------------------------------------------------

struct MergedGraphReader
{
  ReadBuffer<PathNode>            paths;
  ReadBuffer<PathNode::rank_type> labels;
  ReadBuffer<range_type>          from_nodes;

  size_type path, rank, from;

  const DeBruijnGraph*            mapper;
  const sdsl::int_vector<0>*      last_char;

  void init(const MergedGraph& graph, const DeBruijnGraph* _mapper, const sdsl::int_vector<0>* _last_char);
  void init(const MergedGraph& graph, size_type comp);
  void close();

  void seek();

  inline void advance()
  {
    if(this->path + 1 >= this->paths.size()) { return; }
    this->path++;
    this->seek();
  }

  void predecessor(comp_type comp, PathLabel& first, PathLabel& last);

  /*
    Does paths[path + offset] intersect with the given range of labels?
  */
  bool intersect(const PathLabel& first, const PathLabel& last, size_type offset);

  void fromNodes(std::vector<node_type>& results);
};

void
MergedGraphReader::init(const MergedGraph& graph,
  const DeBruijnGraph* _mapper, const sdsl::int_vector<0>* _last_char)
{
  this->paths.init(graph.path_name);
  this->labels.init(graph.rank_name);
  this->from_nodes.init(graph.from_name);

  this->path = this->rank = this->from = 0;
  this->seek();

  this->mapper = _mapper;
  this->last_char = _last_char;
}

void
MergedGraphReader::init(const MergedGraph& graph, size_type comp)
{
  this->paths.init(graph.path_name);
  this->labels.init(graph.rank_name);
  this->from_nodes.init(graph.from_name);

  this->path = graph.next[comp];
  this->from = graph.next_from[comp];
  this->seek();

  this->mapper = 0;
  this->last_char = 0;
}

void
MergedGraphReader::close()
{
  this->paths.close(),
  this->labels.close();
  this->from_nodes.close();

  this->path = this->rank = this->from = 0;
}

void
MergedGraphReader::seek()
{
  this->paths.seek(this->path);
  this->rank = this->paths[this->path].pointer();
  this->labels.seek(this->rank);
  while(this->from < this->from_nodes.size() && this->from_nodes[this->from].first < this->path)
  {
    this->from++;
  }
  this->from_nodes.seek(this->from);
}

void
MergedGraphReader::predecessor(comp_type comp, PathLabel& first, PathLabel& last)
{
  const PathNode& curr = this->paths[this->path];
  size_type i = 0, j = curr.pointer();
  first.first = true; last.first = false;

  // Handle the common prefix of the labels.
  while(i < curr.lcp())
  {
    first.label[i] = last.label[i] = this->mapper->LF(this->labels[j], comp);
    comp = (*(this->last_char))[labels[j]];
    i++; j++;
  }

  // Handle the diverging suffixes of the labels.
  comp_type first_comp = comp, last_comp = comp;
  if(i < curr.order())
  {
    first.label[i] = this->mapper->LF(this->labels[j], first_comp);
    first_comp = (*(this->last_char))[this->labels[j]];
    last.label[i] = this->mapper->LF(this->labels[j + 1], last_comp);
    last_comp = (*(this->last_char))[this->labels[j + 1]];
    i++;
  }
  if(i < PathLabel::LABEL_LENGTH)
  {
    first.label[i] = this->mapper->node_rank(this->mapper->alpha.C[first_comp]);
    last.label[i] = this->mapper->node_rank(this->mapper->alpha.C[last_comp + 1]) - 1;
    i++;
  }

  // Pad the labels.
  while(i < PathLabel::LABEL_LENGTH)
  {
    first.label[i] = 0; last.label[i] = PathLabel::NO_RANK; i++;
  }
}

inline PathLabel
firstLabel(const PathNode& path, ReadBuffer<PathNode::rank_type>& labels)
{
  PathLabel res; res.first = true;
  size_type limit = std::min(path.order(), PathLabel::LABEL_LENGTH);
  for(size_type i = 0; i < limit; i++) { res.label[i] = path.firstLabel(i, labels); }
  for(size_type i = limit; i < PathLabel::LABEL_LENGTH; i++) { res.label[i] = 0; }
  return res;
}

inline PathLabel
lastLabel(const PathNode& path, ReadBuffer<PathNode::rank_type>& labels)
{
  PathLabel res; res.first = false;
  size_type limit = std::min(path.order(), PathLabel::LABEL_LENGTH);
  for(size_type i = 0; i < limit; i++) { res.label[i] = path.lastLabel(i, labels); }
  for(size_type i = limit; i < PathLabel::LABEL_LENGTH; i++) { res.label[i] = PathLabel::NO_RANK; }
  return res;
}

/*
  Does the path node intersect with the given range of labels?
*/
bool
MergedGraphReader::intersect(const PathLabel& first, const PathLabel& last, size_type offset)
{
  PathLabel my_first = firstLabel(this->paths[this->path + offset], this->labels);
  if(my_first <= first)
  {
    PathLabel my_last = lastLabel(this->paths[this->path + offset], this->labels);
    return (first <= my_last);
  }
  else
  {
    return my_first <= last;
  }
}

void
MergedGraphReader::fromNodes(std::vector<node_type>& results)
{
  results.clear();
  results.push_back(this->paths[this->path].from);

  size_type old_pointer = this->from;
  while(this->from < this->from_nodes.size() && this->from_nodes[this->from].first == this->path)
  {
    results.push_back(this->from_nodes[this->from].second); this->from++;
  }
  this->from = old_pointer;

  removeDuplicates(results, false);
}

//------------------------------------------------------------------------------

void
updateRedundant(CounterArray& redundant, const MergedGraph& merged_graph, std::vector<range_type>& ranges)
{
  #pragma omp parallel for schedule(static)
  for(size_type i = 0; i < ranges.size(); i++)
  {
    ranges[i].first = merged_graph.lcp_rmq(ranges[i].first, ranges[i].second);
  }
  for(size_type i = 0; i < ranges.size(); i++) { redundant.increment(ranges[i].first - 1); }
  ranges.clear();
}

GCSA::GCSA(const InputGraph& graph, const ConstructionParameters& parameters, const Alphabet& _alpha)
{
  if(graph.size() == 0) { return; }
  size_type bytes_required = graph.size() * (sizeof(PathNode) + 2 * sizeof(PathNode::rank_type));
  if(bytes_required > parameters.size_limit)
  {
    std::cerr << "GCSA::GCSA(): The input is too large: " << (bytes_required / GIGABYTE_DOUBLE) << " GB" << std::endl;
    std::cerr << "GCSA::GCSA(): Construction aborted" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Extract the keys and build the necessary support structures.
  // FIXME Later: Write the structures to disk until needed?
  std::vector<key_type> keys;
  graph.readKeys(keys);
  DeBruijnGraph mapper(keys, graph.k(), _alpha);
  LCP lcp(keys, graph.k());
  sdsl::int_vector<0> last_char;
  Key::lastChars(keys, last_char);
  sdsl::sd_vector_builder builder(Key::label(keys[keys.size() - 1]) + 1, keys.size());
  for(size_type i = 0; i < keys.size(); i++) { builder.set(Key::label(keys[i])); }
  sdsl::sd_vector<> key_exists(builder);
  sdsl::util::clear(keys);

  // Determine the existing from nodes.
  std::vector<node_type> from_node_buffer;;
  graph.readFrom(from_node_buffer);
  sdsl::sd_vector<> from_nodes(from_node_buffer.begin(), from_node_buffer.end());
  sdsl::sd_vector<>::rank_1_type from_rank;
  sdsl::util::init_support(from_rank, &(from_nodes));
  size_type unique_from_nodes = from_node_buffer.size();
  sdsl::util::clear(from_node_buffer);

  // Create the initial PathGraph.
  PathGraph path_graph(graph, key_exists);
  sdsl::util::clear(key_exists);

  // Prefix-doubling.
#ifdef VERBOSE_STATUS_INFO
  std::cerr << "GCSA::GCSA(): Initial path length: " << path_graph.k() << std::endl;
#endif
  path_graph.prune(lcp, parameters.size_limit);
  for(size_type step = 1; step <= parameters.doubling_steps && path_graph.unsorted > 0; step++)
  {
#ifdef VERBOSE_STATUS_INFO
    std::cerr << "GCSA::GCSA(): Step " << step << " (path length " << path_graph.k() << " -> "
              << (2 * path_graph.k()) << ")" << std::endl;
#endif
    path_graph.extend(parameters.size_limit);
    path_graph.prune(lcp, parameters.size_limit);
  }
  this->header.order = (path_graph.unsorted == 0 ? ~(size_type)0 : path_graph.k());

  // Merge the paths into the nodes of a pruned de Bruijn graph.
  MergedGraph merged_graph(path_graph, mapper, lcp);
  this->header.path_nodes = merged_graph.size();
  path_graph.clear();
  sdsl::util::clear(lcp);

  // Structures used for building GCSA.
  sdsl::int_vector<64> counts(mapper.alpha.sigma, 0); // alpha
  sdsl::int_vector<8> bwt_buffer(merged_graph.size() + merged_graph.size() / 2, 0); // bwt
  this->path_nodes = bit_vector(bwt_buffer.size(), 0);
  CounterArray outdegrees(merged_graph.size()); // edges
  this->sampled_paths = bit_vector(merged_graph.size(), 0);
  std::vector<node_type> sample_buffer; // stored_samples
  this->samples = bit_vector(merged_graph.size() + merged_graph.extra(), 0);

  // Structures used for building OccurrenceCounter.
  // Invariant: The previous occurrence of from node x was at path prev_occ[from_rank(x)] - 1.
  // The RMQs used for updating the redundant array are buffered and done in parallel.
  CounterArray occurrences(merged_graph.size()), redundant(merged_graph.size() - 1);
  sdsl::int_vector<0> prev_occ(unique_from_nodes, 0, bit_length(merged_graph.size()));
  std::vector<range_type> rmq_ranges; rmq_ranges.reserve(RMQ_BUFFER);

  // Read pointers to the MergedGraph files.
  std::vector<MergedGraphReader> reader(mapper.alpha.sigma + 1);
  reader[0].init(merged_graph, &mapper, &last_char);
  for(size_type comp = 0; comp < mapper.alpha.sigma; comp++)
  {
    reader[comp + 1].init(merged_graph, comp);
  }

  // The actual construction.
  PathLabel first, last;
  size_type total_edges = 0, sample_bits = 0;
  std::vector<node_type> pred_from, curr_from;
  for(size_type i = 0; i < merged_graph.size(); i++, reader[0].advance())
  {
    // Find the predecessors.
    size_type indegree = 0, pred_comp = 0;
    bool sample_this = false;
    for(size_type comp = 0; comp < mapper.alpha.sigma; comp++)
    {
      if(!(reader[0].paths[reader[0].path].hasPredecessor(comp))) { continue; }

      // Find the predecessor of paths[i] with comp and the first path intersecting it.
      reader[0].predecessor(comp, first, last);
      if(!(reader[comp + 1].intersect(first, last, 0)))
      {
        reader[comp + 1].advance();
      }

      // Add an edge.
      indegree++; outdegrees.increment(reader[comp + 1].path); total_edges++;
      if(total_edges > bwt_buffer.size()) { bwt_buffer.resize(bwt_buffer.size() + merged_graph.size() / 2); }
      bwt_buffer[total_edges - 1] = comp; counts[comp]++;
      pred_comp = comp; // For sampling.

      // Create additional edges if the next path nodes also intersect with the predecessor.
      while(reader[comp + 1].path + 1 < merged_graph.size() && reader[comp + 1].intersect(first, last, 1))
      {
        reader[comp + 1].advance();
        indegree++; outdegrees.increment(reader[comp + 1].path); total_edges++;
        if(total_edges > bwt_buffer.size()) { bwt_buffer.resize(bwt_buffer.size() + merged_graph.size() / 2); }
        bwt_buffer[total_edges - 1] = comp; counts[comp]++;
      }
    }
    if(total_edges > this->path_nodes.size()) { this->path_nodes.resize(bwt_buffer.size()); }
    this->path_nodes[total_edges - 1] = 1;

    // Get the from nodes and update the occurrences array.
    // Buffer the RMQ ranges needed for updating the redundant array.
    reader[0].fromNodes(curr_from);
    occurrences.increment(i, curr_from.size());
    for(size_type j = 0; j < curr_from.size(); j++)
    {
      size_type temp = from_rank(curr_from[j]);
      if(prev_occ[temp] > 0)
      {
        rmq_ranges.push_back(range_type(prev_occ[temp], i));
        if(rmq_ranges.size() >= RMQ_BUFFER) { updateRedundant(redundant, merged_graph, rmq_ranges); }
      }
      prev_occ[temp] = i + 1;
    }

    /*
      Simple cases for sampling the node:
      - multiple predecessors
      - at the beginning of the source node with no real predecessors
      - at the beginning of a node in the original graph (makes the previous case redundant)
    */
    if(indegree > 1) { sample_this = true; }
    if(reader[0].paths[reader[0].path].hasPredecessor(Alphabet::SINK_COMP)) { sample_this = true; }
    for(size_type k = 0; k < curr_from.size(); k++)
    {
      if(Node::offset(curr_from[k]) == 0) { sample_this = true; break; }
    }

    // Sample if the from nodes cannot be derived from the only predecessor.
    if(!sample_this)
    {
      reader[pred_comp + 1].fromNodes(pred_from);
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
  for(size_type i = 0; i < reader.size(); i++) { reader[i].close(); }
  sdsl::util::clear(last_char); sdsl::util::clear(from_nodes); sdsl::util::clear(prev_occ);
  if(rmq_ranges.size() > 0) { updateRedundant(redundant, merged_graph, rmq_ranges); }
  sdsl::util::clear(rmq_ranges);
  this->header.edges = total_edges;

  // Initialize alpha.
  this->alpha = Alphabet(counts, mapper.alpha.char2comp, mapper.alpha.comp2char);
  sdsl::util::clear(mapper);

  // Initialize counter.
#ifdef VERBOSE_STATUS_INFO
  size_type occ_count = occurrences.sum(), red_count = redundant.sum();
#endif
  this->counter = OccurrenceCounter(occurrences, redundant);
std::cout << "OccurrenceCounter: " << inMegabytes(sdsl::size_in_bytes(this->counter)) << " MB" << std::endl;
  for(size_type i = 0; i < occurrences.size(); i++) { occurrences.decrement(i); }
  this->total_pointers = SadaSparse(occurrences);
std::cout << "SadaSparse / total: " << inMegabytes(sdsl::size_in_bytes(this->total_pointers)) << " MB" << std::endl;
  this->redundant_pointers = SadaSparse(redundant);
std::cout << "SadaSparse / redundant: " << inMegabytes(sdsl::size_in_bytes(this->redundant_pointers)) << " MB" << std::endl;
  sdsl::util::clear(occurrences); sdsl::util::clear(redundant);

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
  std::cerr << "GCSA::GCSA(): " << this->size() << " paths, " << this->edgeCount() << " edges" << std::endl;
  std::cerr << "GCSA::GCSA(): " << occ_count << " pointers (" << red_count << " redundant)" << std::endl;
  std::cerr << "GCSA::GCSA(): " << this->sampleCount() << " samples at "
            << this->sampledPositions() << " positions" << std::endl;
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
locateFailure(const std::vector<node_type>& expected, const std::vector<node_type>& occs)
{
  std::cerr << "GCSA::verifyIndex(): Expected ";
  printOccs(expected, std::cerr) << std::endl;
  std::cerr << "GCSA::verifyIndex(): Got ";
  printOccs(occs, std::cerr) << std::endl;
}

bool
printFailure(size_type& failure_count)
{
  if(failure_count == GCSA::MAX_ERRORS)
  {
    std::cerr << "GCSA::verifyIndex(): There were further errors" << std::endl;
  }
  failure_count++;
  return (failure_count <= GCSA::MAX_ERRORS);
}

struct KMerSplitComparator
{
  inline bool operator() (const KMer& left, const KMer& right) const
    {
      return (Key::label(left.key) != Key::label(right.key));
    }
};

bool
GCSA::verifyIndex(const InputGraph& graph) const
{
  std::vector<KMer> kmers; graph.read(kmers);
  return this->verifyIndex(kmers, graph.k());
}

bool
GCSA::verifyIndex(std::vector<KMer>& kmers, size_type kmer_length) const
{
  double start = readTimer();

  size_type threads = omp_get_max_threads();
  parallelQuickSort(kmers.begin(), kmers.end());
  KMerSplitComparator k_comp;
  std::vector<range_type> bounds = getBounds(kmers, threads, k_comp);
  assert(bounds.size() == threads);

  size_type fails = 0, unique = 0;
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
      #pragma omp atomic
      unique++;

      std::string kmer = Key::decode(kmers[i].key, kmer_length, alpha);
      size_type endmarker_pos = kmer.find('$'); // The actual kmer ends at the first endmarker.
      if(endmarker_pos != std::string::npos) { kmer = kmer.substr(0, endmarker_pos + 1); }

      range_type range = this->find(kmer);
      if(Range::empty(range))
      {
        #pragma omp critical
        {
          if(printFailure(fails))
          {
            std::cerr << "GCSA::verifyIndex(): find(" << kmer << ") returned empty range" << std::endl;
          }
        }
        i = next; continue;
      }

      std::vector<node_type> expected;
      for(size_type j = i; j < next; j++) { expected.push_back(kmers[j].from); }
      removeDuplicates(expected, false);
      size_type unique_count = this->count(range);
      if(unique_count != expected.size())
      {
        #pragma omp critical
        {
          if(printFailure(fails))
          {
            std::cerr << "GCSA::verifyIndex(): count" << range << " failed: Expected "
                      << expected.size() << " occurrences, got " << unique_count << std::endl;
          }
        }
        i = next; continue;
      }

      std::vector<node_type> occs;
      this->locate(range, occs);
      if(occs.size() != expected.size())
      {
        #pragma omp critical
        {
          if(printFailure(fails))
          {
            std::cerr << "GCSA::verifyIndex(): locate(" << kmer << ") failed: Expected "
                      << expected.size() << " occurrences, got " << occs.size() << std::endl;
            locateFailure(expected, occs);
          }
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
              if(printFailure(fails))
              {
                std::cerr << "GCSA::verifyIndex(): locate(" << kmer << ") failed: Expected "
                          << Node::decode(expected[j]) << ", got " << Node::decode(occs[j]) << std::endl;
                locateFailure(expected, occs);
              }
            }
            break;
          }
        }
      }

      i = next;
    }
  }

  double seconds = readTimer() - start;
  std::cout << "Queried the index with " << unique << " patterns in " << seconds << " seconds ("
            << (unique / seconds) << " patterns / second)" << std::endl;
  if(fails == 0)
  {
    std::cout << "Index verification complete" << std::endl;
  }
  else
  {
    std::cout << "Index verification failed for " << fails << " patterns" << std::endl;
  }
  std::cout << std::endl;

  return fails == 0;
}

//------------------------------------------------------------------------------

struct ReverseTrieNode
{
  range_type range;
  size_type  depth;
  bool       ends_at_sink;

  const static size_type SEED_LENGTH = 5; // Create all patterns of that length before parallelizing.

  ReverseTrieNode() : range(0, 0), depth(0), ends_at_sink(false) {}
  ReverseTrieNode(range_type rng, size_type d, bool sink) : range(rng), depth(d), ends_at_sink(sink) {}
};

struct KMerCounter
{
  size_type count, limit;

  KMerCounter(size_type _limit) : count(0), limit(_limit) {}

  inline void add(const ReverseTrieNode& node)
  {
    if(node.ends_at_sink || node.depth >= this->limit) { this->count++; }
  }
};

struct SeedCollector
{
  size_type count, limit;
  std::vector<ReverseTrieNode> seeds;

  SeedCollector(size_type _limit) : count(0), limit(_limit), seeds() {}

  inline void add(const ReverseTrieNode& node)
  {
    if(node.depth >= this->limit) { this->seeds.push_back(node); }
    else if(node.ends_at_sink) { this->count++; } // Seeds ending at the sink are counted later.
  }
};

template<class Handler>
void
processSubtree(const GCSA& index, std::stack<ReverseTrieNode>& node_stack, Handler& handler)
{
  std::vector<range_type> predecessors(index.alpha.sigma);
  while(!(node_stack.empty()))
  {
    ReverseTrieNode curr = node_stack.top(); node_stack.pop();
    if(Range::empty(curr.range)) { continue; }
    handler.add(curr);
    if(curr.depth < handler.limit)
    {
      index.LF(curr.range, predecessors);
      for(size_type comp = 1; comp + 1 < index.alpha.sigma; comp++)
      {
        node_stack.push(ReverseTrieNode(predecessors[comp], curr.depth + 1, curr.ends_at_sink));
      }
    }
  }
}

size_type
GCSA::countKMers(size_type k, bool force) const
{
  if(k == 0) { return 1; }
  if(k > this->order())
  {
    if(force)
    {
      std::cerr << "GCSA::countKMers(): Warning: The value of k (" << k
                << ") is greater than the order of the graph (" << this->order() << ")" << std::endl;
    }
    else
    {
      std::cerr << "GCSA::countKMers(): The value of k (" << k
                << ") is greater than the order of the graph (" << this->order() << ")" << std::endl;
      return 0;
    }
  }

  // Create an array of seed kmers of length ReverseTrieNode::SEED_LENGTH.
  size_type result = 0;
  std::vector<ReverseTrieNode> seeds;
  {
    std::stack<ReverseTrieNode> node_stack;
    node_stack.push(ReverseTrieNode(this->charRange(0), 1, true));
    for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++)
    {
      node_stack.push(ReverseTrieNode(this->charRange(comp), 1, false));
    }
    SeedCollector collector(std::min(k, ReverseTrieNode::SEED_LENGTH));
    processSubtree(*this, node_stack, collector);
    result = collector.count;
    seeds = collector.seeds;
  }

  // Extend the seeds in parallel.
  #pragma omp parallel for schedule (dynamic, 1)
  for(size_type i = 0; i < seeds.size(); i++)
  {
    std::stack<ReverseTrieNode> node_stack;
    node_stack.push(seeds[i]);
    KMerCounter counter(k);
    processSubtree(*this, node_stack, counter);
    #pragma omp atomic
    result += counter.count;
  }

  return result;
}

//------------------------------------------------------------------------------

void
GCSA::copy(const GCSA& g)
{
  this->header = g.header;

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
GCSA::setVectors()
{
  this->path_rank.set_vector(&(this->path_nodes));
  this->path_select.set_vector(&(this->path_nodes));

  this->edge_rank.set_vector(&(this->edges));
  this->edge_select.set_vector(&(this->edges));

  this->sampled_path_rank.set_vector(&(this->sampled_paths));

  this->sample_select.set_vector(&(this->samples));
}

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
GCSA::LF(range_type range, std::vector<range_type>& results) const
{
  for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++) { results[comp] = Range::empty_range(); }
  range = this->bwtRange(range);

  if(range.first == range.second) // Single edge.
  {
    auto temp = this->bwt.inverse_select(range.first);
    size_type path_node = this->edge_rank(this->alpha.C[temp.second] + temp.first);
    results[temp.second] = range_type(path_node, path_node);
  }
  else if(Range::length(range) <= SHORT_RANGE) // Use brute force for a few edges.
  {
    for(size_type i = range.first; i <= range.second; i++)
    {
      auto temp = this->bwt.inverse_select(i);
      if(Range::empty(results[temp.second]))
      {
        size_type edge_pos = this->alpha.C[temp.second] + temp.first;
        results[temp.second] = range_type(edge_pos, edge_pos);
      }
      else { results[temp.second].second++; }
    }
    for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++)
    {
      if(!Range::empty(results[comp])) { results[comp] = this->pathNodeRange(results[comp]); }
    }
  }
  else  // General case.
  {
    for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++)
    {
      results[comp] = gcsa::LF(this->bwt, this->alpha, range, comp);
      if(!Range::empty(results[comp])) { results[comp] = this->pathNodeRange(results[comp]); }
    }
  }
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
  while(!(this->sampled(path_node)))
  {
    path_node = this->LF(path_node);
    steps++;
  }

  range_type sample_range = this->sampleRange(path_node);
  for(size_type i = sample_range.first; i <= sample_range.second; i++)
  {
    results.push_back(this->sample(i) + steps);
  }
}

//------------------------------------------------------------------------------

} // namespace gcsa
