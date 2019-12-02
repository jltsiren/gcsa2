/*
  Copyright (c) 2018, 2019 Jouni Siren
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

#include <gcsa/algorithms.h>
#include <gcsa/internal.h>
#include <gcsa/path_graph.h>

#include <random>
#include <unordered_set>

namespace gcsa
{

//------------------------------------------------------------------------------

// Other class variables

const std::string GCSA::EXTENSION = ".gcsa";

//------------------------------------------------------------------------------

GCSA::GCSA() :
  header(), alpha(),
  fast_bwt(this->alpha.sigma), fast_rank(this->alpha.sigma),
  sparse_bwt(this->alpha.sigma), sparse_rank(this->alpha.sigma),
  edges(), edge_rank(),
  sampled_paths(), sampled_path_rank(),
  stored_samples(), samples(), sample_select(),
  extra_pointers(), redundant_pointers()
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
    this->alpha.swap(another.alpha);

    this->fast_bwt.swap(another.fast_bwt);
    this->fast_rank.swap(another.fast_rank);

    this->sparse_bwt.swap(another.sparse_bwt);
    this->sparse_rank.swap(another.sparse_rank);

    this->edges.swap(another.edges);
    sdsl::util::swap_support(this->edge_rank, another.edge_rank, &(this->edges), &(another.edges));

    this->sampled_paths.swap(another.sampled_paths);
    sdsl::util::swap_support(this->sampled_path_rank, another.sampled_path_rank, &(this->sampled_paths), &(another.sampled_paths));

    this->stored_samples.swap(another.stored_samples);
    this->samples.swap(another.samples);
    sdsl::util::swap_support(this->sample_select, another.sample_select, &(this->samples), &(another.samples));

    this->extra_pointers.swap(another.extra_pointers);
    this->redundant_pointers.swap(another.redundant_pointers);

    this->setVectors();
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
    this->alpha = std::move(source.alpha);

    this->fast_bwt = std::move(source.fast_bwt);
    this->fast_rank = std::move(source.fast_rank);

    this->sparse_bwt = std::move(source.sparse_bwt);
    this->sparse_rank = std::move(source.sparse_rank);

    this->edges = std::move(source.edges);
    this->edge_rank = std::move(source.edge_rank);

    this->sampled_paths = std::move(source.sampled_paths);
    this->sampled_path_rank = std::move(source.sampled_path_rank);

    this->stored_samples = std::move(source.stored_samples);
    this->samples = std::move(source.samples);
    this->sample_select = std::move(source.sample_select);

    this->extra_pointers = std::move(source.extra_pointers);
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
  written_bytes += this->alpha.serialize(out, child, "alpha");

  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->fast_bwt[comp].serialize(out, child, "fast_bwt");
  }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->fast_rank[comp].serialize(out, child, "fast_rank");
  }

  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->sparse_bwt[comp].serialize(out, child, "sparse_bwt");
  }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->sparse_rank[comp].serialize(out, child, "sparse_rank");
  }

  written_bytes += this->edges.serialize(out, child, "edges");
  written_bytes += this->edge_rank.serialize(out, child, "edge_rank");

  written_bytes += this->sampled_paths.serialize(out, child, "sampled_paths");
  written_bytes += this->sampled_path_rank.serialize(out, child, "sampled_path_rank");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  written_bytes += this->extra_pointers.serialize(out, child, "extra_pointers");
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
  this->alpha.load(in);

  this->fast_bwt.resize(this->alpha.sigma); this->fast_rank.resize(this->alpha.sigma);
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->fast_bwt[comp].load(in); }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->fast_rank[comp].load(in, &(this->fast_bwt[comp])); }

  this->sparse_bwt.resize(this->alpha.sigma); this->sparse_rank.resize(this->alpha.sigma);
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->sparse_bwt[comp].load(in); }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->sparse_rank[comp].load(in, &(this->sparse_bwt[comp])); }

  this->edges.load(in);
  this->edge_rank.load(in, &(this->edges));

  this->sampled_paths.load(in);
  this->sampled_path_rank.load(in, &(this->sampled_paths));

  this->stored_samples.load(in);
  this->samples.load(in);
  this->sample_select.load(in, &(this->samples));

  this->extra_pointers.load(in);
  this->redundant_pointers.load(in);
}

void
GCSA::copy(const GCSA& source)
{
  this->header = source.header;
  this->alpha = source.alpha;

  this->fast_bwt = source.fast_bwt;
  this->fast_rank = source.fast_rank;

  this->sparse_bwt = source.sparse_bwt;
  this->sparse_rank = source.sparse_rank;

  this->edges = source.edges;
  this->edge_rank = source.edge_rank;

  this->sampled_paths = source.sampled_paths;
  this->sampled_path_rank = source.sampled_path_rank;

  this->stored_samples = source.stored_samples;
  this->samples = source.samples;
  this->sample_select = source.sample_select;

  this->extra_pointers = source.extra_pointers;
  this->redundant_pointers = source.redundant_pointers;

  this->setVectors();
}

void
GCSA::setVectors()
{
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    this->fast_rank[comp].set_vector(&(this->fast_bwt[comp]));
    this->sparse_rank[comp].set_vector(&(this->sparse_bwt[comp]));
  }

  this->edge_rank.set_vector(&(this->edges));

  this->sampled_path_rank.set_vector(&(this->sampled_paths));

  this->sample_select.set_vector(&(this->samples));
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

  void fromNodes(std::vector<node_type>& results, const NodeMapping& mapping);
};

void
MergedGraphReader::init(const MergedGraph& graph,
  const DeBruijnGraph* _mapper, const sdsl::int_vector<0>* _last_char)
{
  this->paths.open(graph.path_name);
  this->labels.open(graph.rank_name);
  this->from_nodes.open(graph.from_name);

  this->path = this->rank = this->from = 0;
  this->seek();

  this->mapper = _mapper;
  this->last_char = _last_char;
}

void
MergedGraphReader::init(const MergedGraph& graph, size_type comp)
{
  this->paths.open(graph.path_name);
  this->labels.open(graph.rank_name);
  this->from_nodes.open(graph.from_name);

  this->path = graph.next[comp];
  this->from = graph.next_from[comp];
  this->seek();

  this->mapper = nullptr;
  this->last_char = nullptr;
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
MergedGraphReader::fromNodes(std::vector<node_type>& results, const NodeMapping& mapping)
{
  results.clear();
  results.push_back(this->paths[this->path].from);

  size_type old_pointer = this->from;
  while(this->from < this->from_nodes.size() && this->from_nodes[this->from].first == this->path)
  {
    results.push_back(this->from_nodes[this->from].second); this->from++;
  }
  this->from = old_pointer;

  Node::map(results, mapping);
  removeDuplicates(results, false);
}

//------------------------------------------------------------------------------

GCSA::GCSA(InputGraph& graph, const ConstructionParameters& parameters) :
  GCSA()
{
  double start = readTimer();

  if(graph.size() == 0) { return; }
  size_type bytes_required = graph.size() * (sizeof(PathNode) + 2 * sizeof(PathNode::rank_type));
  if(bytes_required > parameters.getLimitBytes())
  {
    std::cerr << "GCSA::GCSA(): The input is too large: " << (bytes_required / GIGABYTE_DOUBLE) << " GB" << std::endl;
    std::cerr << "GCSA::GCSA(): Construction aborted" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Extract the keys and build the necessary support structures.
  // FIXME Later: Write the structures to disk until needed?
  std::vector<key_type> keys;
  graph.readKeys(keys);
  DeBruijnGraph mapper(keys, graph.k(), graph.alpha);
  LCP lcp(keys, graph.k());
  sdsl::int_vector<0> last_char;
  Key::lastChars(keys, last_char);
  sdsl::int_vector<0> distinct_labels(keys.size(), 0, bit_length(Key::label(keys.back())));
  for(size_type i = 0; i < keys.size(); i++) { distinct_labels[i] = Key::label(keys[i]); }
  sdsl::util::clear(keys);

  // Determine the existing start nodes. Because the information is only used for index
  // construction, we use NodeMapping to map the node ids used for construction to the
  // node ids reported by locate().
  std::vector<node_type> from_node_buffer;
  graph.readFrom(from_node_buffer, true);
  sdsl::sd_vector<> from_nodes(from_node_buffer.begin(), from_node_buffer.end());
  sdsl::sd_vector<>::rank_1_type from_rank;
  sdsl::util::init_support(from_rank, &(from_nodes));
  size_type unique_from_nodes = from_node_buffer.size();
  sdsl::util::clear(from_node_buffer);

  // Create the initial PathGraph.
  PathGraph path_graph(graph, distinct_labels);
  sdsl::util::clear(distinct_labels);
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Preprocessing: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
    start = stop;
  }

  // Prefix-doubling.
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): Prefix-doubling from path length " << path_graph.k() << std::endl;
  }
  for(size_type step = 1; step <= parameters.getSteps(); step++)
  {
    if(Verbosity::level >= Verbosity::BASIC)
    {
      std::cerr << "GCSA::GCSA(): Step " << step << " (path length " << path_graph.k() << " -> "
                << (2 * path_graph.k()) << ")" << std::endl;
    }
    path_graph.prune(lcp, parameters.getLimitBytes() - path_graph.bytes());
    path_graph.extend(parameters.getLimitBytes() - path_graph.bytes());
  }
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Prefix-doubling: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
    start = stop;
  }

  // Merge the paths into the nodes of a maximally pruned de Bruijn graph.
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): Merging the paths" << std::endl;
  }
  MergedGraph merged_graph(path_graph, mapper, lcp, parameters.getLimitBytes() - path_graph.bytes());
  this->header.path_nodes = merged_graph.size();
  this->header.order = merged_graph.k();
  path_graph.clear();
  sdsl::util::clear(lcp);
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Merging: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
    start = stop;
  }

  // Structures used for building GCSA.
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): Building the index" << std::endl;
  }
  sdsl::int_vector<64> counts(graph.alpha.sigma, 0); // alpha
  std::vector<bit_vector> bwt(graph.alpha.sigma); // fast_bwt, sparse_bwt
  for(size_type comp = 0; comp < bwt.size(); comp++) { bwt[comp] = bit_vector(merged_graph.size(), 0); }
  CounterArray outdegrees(merged_graph.size(), 4); // edges
  bit_vector sampled_positions(merged_graph.size(), 0); // sampled_paths
  std::vector<node_type> sample_buffer; // stored_samples
  this->samples = bit_vector(merged_graph.size() + merged_graph.extra(), 0);

  // Structures used for building counting support.
  // Invariant: The previous occurrence of from node x was at path prev_occ[from_rank(x)] - 1.
  CounterArray occurrences(merged_graph.size(), 4), redundant(merged_graph.size() - 1, 4);
  sdsl::int_vector<0> prev_occ(unique_from_nodes, 0, bit_length(merged_graph.size()));
  std::vector<size_type> node_lcp, first_time, last_time;

  // Read pointers to the MergedGraph files.
  std::vector<MergedGraphReader> reader(graph.alpha.sigma + 1);
  reader[0].init(merged_graph, &mapper, &last_char);
  for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
  {
    reader[comp + 1].init(merged_graph, comp);
  }
  ReadBuffer<uint8_t> lcp_array; lcp_array.open(merged_graph.lcp_name);

  // The actual construction.
  PathLabel first, last;
  size_type total_edges = 0, sample_bits = 0;
  std::vector<node_type> pred_from, curr_from;
  for(size_type i = 0; i < merged_graph.size(); i++, reader[0].advance())
  {
    // Find the predecessors.
    size_type indegree = 0, pred_comp = 0;
    bool sample_this = false;
    for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
    {
      if(!(reader[0].paths[reader[0].path].hasPredecessor(comp))) { continue; }

      // Find the predecessor of paths[i] with comp and the path intersecting it.
      reader[0].predecessor(comp, first, last);
      if(!(reader[comp + 1].intersect(first, last, 0)))
      {
        reader[comp + 1].advance();
      }

      // Add the edge.
      bwt[comp][i] = 1; counts[comp]++;
      indegree++; outdegrees.increment(reader[comp + 1].path); total_edges++;
      pred_comp = comp; // For sampling.
    }

    /*
      Get the start nodes and update the occurrences/redundant arrays.

      We traverse the ST in inorder using the LCP array. For each internal node, we
      record the LCP value and the first and the last times (positions) we have
      encountered that value within the subtree. If we have encountered the current
      from node before, the LCA of the previous and current occurrences is the
      highest ST node we have encountered after the previous occurrence. We then
      increment the redundant array at the first encounter with that node.
    */
    reader[0].fromNodes(curr_from, graph.mapping);
    occurrences.increment(i, curr_from.size() - 1);
    lcp_array.seek(i);
    size_type curr_lcp = lcp_array[i] + (i > 0 ? 1 : 0); // Handle LCP[0] as -1.
    while(!(node_lcp.empty()) && node_lcp.back() > curr_lcp)
    {
      node_lcp.pop_back(); first_time.pop_back(); last_time.pop_back();
    }
    if(!(node_lcp.empty()) && node_lcp.back() == curr_lcp) { last_time.back() = i; }
    else { node_lcp.push_back(curr_lcp); first_time.push_back(i); last_time.push_back(i); }
    for(size_type j = 0; j < curr_from.size(); j++)
    {
      size_type temp = from_rank(curr_from[j]);
      if(prev_occ[temp] > 0)
      {
        size_type pos = std::lower_bound(last_time.begin(), last_time.end(), prev_occ[temp]) - last_time.begin();
        redundant.increment(first_time[pos] - 1);
      }
      prev_occ[temp] = i + 1;
    }

    /*
      Simple cases for sampling the node:
      - multiple predecessors
      - at the beginning of the source node with no real predecessors
      - if a from node is divisible by the sample period
    */
    if(indegree > 1) { sample_this = true; }
    if(reader[0].paths[reader[0].path].hasPredecessor(Alphabet::SINK_COMP)) { sample_this = true; }
    for(size_type k = 0; k < curr_from.size(); k++)
    {
      if(curr_from[k] % parameters.getSamplePeriod() == 0) { sample_this = true; break; }
    }

    // Sample if the start nodes cannot be derived from the only predecessor.
    if(!sample_this)
    {
      reader[pred_comp + 1].fromNodes(pred_from, graph.mapping);
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
      sampled_positions[i] = 1;
      for(size_type k = 0; k < curr_from.size(); k++)
      {
        sample_bits = std::max(sample_bits, bit_length(curr_from[k]));
        sample_buffer.push_back(curr_from[k]);
      }
      this->samples[sample_buffer.size() - 1] = 1;
    }
  }
  for(size_type i = 0; i < reader.size(); i++) { reader[i].close(); }
  lcp_array.close();
  sdsl::util::clear(last_char); sdsl::util::clear(from_nodes); sdsl::util::clear(prev_occ);
  this->header.edges = total_edges;

  // Initialize alpha.
  this->alpha = Alphabet(counts, graph.alpha.char2comp, graph.alpha.comp2char);
  sdsl::util::clear(mapper);

  // Initialize extra_pointers and redundant_pointers.
  size_type occ_count = occurrences.sum() + occurrences.size(), red_count = redundant.sum();
  this->extra_pointers = SadaSparse(occurrences);
  this->redundant_pointers = SadaCount(redundant);
  sdsl::util::clear(occurrences); sdsl::util::clear(redundant);

  // Initialize bwt.
  this->fast_bwt.resize(this->alpha.sigma); this->fast_rank.resize(this->alpha.sigma);
  this->sparse_bwt.resize(this->alpha.sigma); this->sparse_rank.resize(this->alpha.sigma);
  this->sparse_bwt[0] = bwt[0]; sdsl::util::clear(bwt[0]);
  for(size_type comp = 1; comp <= graph.alpha.fast_chars; comp++)
  {
    this->fast_bwt[comp] = bwt[comp]; sdsl::util::clear(bwt[comp]);
  }
  for(size_type comp = graph.alpha.fast_chars + 1; comp < graph.alpha.sigma; comp++)
  {
    this->sparse_bwt[comp] = bwt[comp]; sdsl::util::clear(bwt[comp]);
  }

  // Initialize bitvectors (edges, sampled_positions, samples).
  bit_vector edge_buffer(total_edges, 0); total_edges = 0;
  for(size_type i = 0; i < merged_graph.size(); i++)
  {
    total_edges += outdegrees[i];
    edge_buffer[total_edges - 1] = 1;
  }
  outdegrees.clear();
  this->edges = edge_buffer; sdsl::util::clear(edge_buffer);
  this->sampled_paths = sampled_positions; sdsl::util::clear(sampled_positions);
  this->samples.resize(sample_buffer.size());
  this->initSupport();

  // Initialize stored_samples.
  this->stored_samples = sdsl::int_vector<0>(sample_buffer.size(), 0, sample_bits);
  for(size_type i = 0; i < sample_buffer.size(); i++) { this->stored_samples[i] = sample_buffer[i]; }
  sdsl::util::clear(sample_buffer);

  // Transfer the LCP array from MergedGraph to InputGraph.
  TempFile::remove(graph.lcp_name);
  graph.lcp_name = merged_graph.lcp_name;
  merged_graph.lcp_name.clear();

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Construction: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
  }
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): " << this->size() << " paths, " << this->edgeCount() << " edges" << std::endl;
    std::cerr << "GCSA::GCSA(): " << occ_count << " pointers (" << red_count << " redundant)" << std::endl;
    std::cerr << "GCSA::GCSA(): " << this->sampleCount() << " samples at "
              << this->sampledPositions() << " positions" << std::endl;
  }
}

void
GCSA::initSupport()
{
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    sdsl::util::init_support(this->fast_rank[comp], &(this->fast_bwt[comp]));
    sdsl::util::init_support(this->sparse_rank[comp], &(this->sparse_bwt[comp]));
  }

  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->sampled_path_rank, &(this->sampled_paths));
  sdsl::util::init_support(this->sample_select, &(this->samples));
}

//------------------------------------------------------------------------------

void
GCSA::LF_fast(range_type range, std::vector<range_type>& results) const
{
  for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++) { results[comp] = Range::empty_range(); }
  if(Range::empty(range)) { return; }

  if(range.first == range.second) // Single path node.
  {
    for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++)
    {
      if(this->fast_bwt[comp][range.first])
      {
        results[comp].first = results[comp].second = this->edge_rank(this->LF(this->fast_rank, range.first, comp));
      }
    }
  }
  else  // General case.
  {
    for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++)
    {
      results[comp] = this->LF(this->fast_rank, range, comp);
      if(!Range::empty(results[comp])) { results[comp] = this->pathNodeRange(results[comp]); }
    }
  }
}

void
GCSA::LF_all(range_type range, std::vector<range_type>& results) const
{
  for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++) { results[comp] = Range::empty_range(); }
  if(Range::empty(range)) { return; }

  if(range.first == range.second) // Single path node.
  {
    for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++)
    {
      if(this->fast_bwt[comp][range.first])
      {
        results[comp].first = results[comp].second = this->edge_rank(this->LF(this->fast_rank, range.first, comp));
      }
    }
    for(size_type comp = this->alpha.fast_chars + 1; comp + 1 < this->alpha.sigma; comp++)
    {
      if(this->sparse_bwt[comp][range.first])
      {
        results[comp].first = results[comp].second = this->edge_rank(this->LF(this->sparse_rank, range.first, comp));
      }
    }
  }
  else  // General case.
  {
    for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++)
    {
      results[comp] = this->LF(range, comp);
    }
  }
}

//------------------------------------------------------------------------------

size_type
GCSA::count(range_type range) const
{
  if(Range::empty(range) || range.second >= this->size()) { return 0; }
  size_type res = this->extra_pointers.count(range.first, range.second) + Range::length(range);
  if(range.second > range.first) { res -= this->redundant_pointers.count(range.first, range.second - 1); }
  return res;
}

//------------------------------------------------------------------------------

void
GCSA::locate(size_type path_node, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { results.clear(); }
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
  if(!append) { results.clear(); }
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
GCSA::locate(range_type range, size_type max_positions, std::vector<node_type>& results) const
{
  results.clear();

  size_type total_positions = this->count(range);
  if(total_positions <= 0) { return; }
  max_positions = std::min(max_positions, total_positions);

  std::mt19937_64 rng(range.first ^ range.second);
  if(max_positions >= total_positions / 2)  // Just locate everything.
  {
    this->locate(range, results);
  }
  else  // Locate at random positions until we have enough distinct occurrences.
  {
    std::unordered_set<node_type, size_type(*)(size_type)> found(16, wang_hash_64);
    while(found.size() < max_positions)
    {
      size_type pos = range.first + rng() % Range::length(range);
      this->locateInternal(pos, results);
      for(node_type node : results) { found.insert(node); }
      results.clear();
    }
    for(node_type node : found) { results.push_back(node); }
  }

  // If there are too many results, select a random subset of them.
  if(results.size() > max_positions)
  {
    deterministicShuffle(results, rng, false);
    results.resize(max_positions);
  }
  sequentialSort(results.begin(), results.end());
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

  size_type sample = this->firstSample(path_node);
  do
  {
    results.push_back(this->sample(sample) + steps); sample++;
  }
  while(!(this->lastSample(sample - 1)));
}

//------------------------------------------------------------------------------

} // namespace gcsa
