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

#include <cstring>
#include <deque>

// For memory mapped files.
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include "path_graph.h"

namespace gcsa
{

//------------------------------------------------------------------------------

struct PriorityNode
{
  typedef PathLabel::rank_type rank_type;

  const static size_type LABEL_LENGTH = PathLabel::LABEL_LENGTH;
  const static rank_type NO_RANK = PathLabel::NO_RANK;

  rank_type file;
  rank_type label[LABEL_LENGTH + 1];
  PathNode  node;

  inline bool eof() const { return (this->label[0] == NO_RANK); }

  // Compare by first labels.
  inline bool operator< (const PriorityNode& another) const
  {
    size_type order = std::min(this->node.order(), another.node.order());
    for(size_type i = 0; i < order; i++)
    {
      if(this->label[i] != another.label[i]) { return (this->label[i] < another.label[i]); }
    }
    return (this->node.order() < another.node.order());
  }

  inline void load(std::vector<std::ifstream>& files)
  {
    this->node = PathNode(files[this->file], this->label);
    if(files[this->file].eof())
    {
      this->node.setOrder(1);
      this->node.setLCP(1);
      this->label[0] = NO_RANK;
    }
  }

  inline rank_type firstLabel(size_type i) const { return this->label[i]; }

  inline void serialize(std::ostream& out) const
  {
    this->node.serialize(out, this->label);
  }

  inline void serialize(std::ostream& path_file, std::ostream& label_file, size_type ptr)
  {
    this->node.serialize(path_file, label_file, this->label, ptr);
  }

  inline void serialize(std::vector<std::ofstream>& files) const
  {
    this->node.serialize(files[this->file], this->label);
  }

  inline size_type bytes() const { return this->node.bytes(); }
};

//------------------------------------------------------------------------------

struct PathGraphBuilder
{
  PathGraph graph;
  std::vector<std::ofstream> files;
  size_type limit;  // Bytes of disk space.

  const static size_type WRITE_BUFFER_SIZE = MEGABYTE; // PathNodes per thread.

  PathGraphBuilder(size_type file_count, size_type path_order, size_type size_limit);
  void close();

  /*
    The single write is not thread safe, while the batch write is. The file number is
    assumed to be valid.
  */
  void write(const PriorityNode& path);
  void write(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file);

  void sort(size_type file);
};

PathGraphBuilder::PathGraphBuilder(size_type file_count, size_type path_order, size_type size_limit) :
  graph(file_count, path_order), files(file_count), limit(size_limit)
{
  for(size_type file = 0; file < this->files.size(); file++)
  {
    this->files[file].open(this->graph.filenames[file].c_str(), std::ios_base::binary);
    if(!(this->files[file]))
    {
      std::cerr << "PathGraphBuilder::PathGraphBuilder(): Cannot open output file "
                << this->graph.filenames[file] << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

void
PathGraphBuilder::close()
{
  for(size_type file = 0; file < this->files.size(); file++)
  {
    this->files[file].close();
  }
}

void
PathGraphBuilder::write(const PriorityNode& path)
{
  if(this->graph.bytes() + path.bytes() > limit)
  {
    std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  path.serialize(this->files);
  this->graph.sizes[path.file]++; this->graph.path_count++;
  this->graph.rank_counts[path.file] += path.node.ranks();
  this->graph.rank_count += path.node.ranks();
}

void
PathGraphBuilder::write(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file)
{
  size_type bytes_required = paths.size() * sizeof(PathNode) + labels.size() * sizeof(PathNode::rank_type);
  #pragma omp critical
  {
    if(bytes_required + this->graph.bytes() > limit)
    {
      std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for(auto& path : paths) { path.serialize(this->files[file], labels); }
    this->graph.sizes[file] += paths.size(); this->graph.path_count += paths.size();
    this->graph.rank_counts[file] += labels.size(); this->graph.rank_count += labels.size();
  }
  paths.clear(); labels.clear();
}

void
PathGraphBuilder::sort(size_type file)
{
  this->files[file].close();

  std::vector<PathNode> paths;
  std::vector<PathNode::rank_type> labels;
  this->graph.read(paths, labels, file);

  PathFirstComparator first_c(labels);
  parallelQuickSort(paths.begin(), paths.end(), first_c);

  this->files[file].open(this->graph.filenames[file].c_str(), std::ios_base::binary);
  for(auto& path : paths) { path.serialize(this->files[file], labels); }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "PathGraphBuilder::sort(): File " << file << ": Sorted "
            << this->graph.sizes[file] << " paths" << std::endl;
#endif
}

//------------------------------------------------------------------------------

struct PriorityQueue
{
  std::vector<PriorityNode> data;

  explicit PriorityQueue(size_type n);

  inline size_type size() const { return this->data.size(); }
  inline static size_type parent(size_type i) { return (i - 1) / 2; }
  inline static size_type left(size_type i) { return 2 * i + 1; }
  inline static size_type right(size_type i) { return 2 * i + 2; }

  inline size_type smaller(size_type i, size_type j) const
  {
    return (this->data[j] < this->data[i] ? j : i);
  }

  inline void down(size_type i)
  {
    while(left(i) < this->size())
    {
      size_type next = this->smaller(i, left(i));
      if(right(i) < this->size()) { next = this->smaller(next, right(i)); }
      if(next == i) { return; }
      std::swap(this->data[i], this->data[next]);
      i = next;
    }
  }

  inline PriorityNode& operator[] (size_type i) { return this->data[i]; }
  inline const PriorityNode& operator[] (size_type i) const { return this->data[i]; }

  void heapify();
};

PriorityQueue::PriorityQueue(size_type n) :
  data(n)
{
}

void
PriorityQueue::heapify()
{
  if(this->size() <= 1) { return; }

  size_type i = parent(this->size() - 1);
  while(true)
  {
    this->down(i);
    if(i == 0) { break; }
    i--;
  }
}

//------------------------------------------------------------------------------

struct PathGraphMerger
{
  const PathGraph& graph;
  std::vector<std::ifstream> files;

  PriorityQueue inputs;
  std::deque<PriorityNode> buffer;

  const static size_type BUFFER_SIZE = MEGABYTE;  // Minimum size in nodes.

  explicit PathGraphMerger(const PathGraph& path_graph);

  void fill(size_type limit = BUFFER_SIZE); // Resize buffer to the new limit if currently smaller.
  inline void readItem(size_type i) { this->fill(std::max(i + 1, BUFFER_SIZE)); }

  void close();

  inline range_type firstRange() { return this->nextRange(Range::empty_range()); }

  // Call readItem(i) first.
  inline bool atEnd(size_type i) const { return (i >= this->buffer.size()); }
  inline bool atEnd(range_type range) const { return this->atEnd(range.first); }

  /*
    Iterates through ranges of paths with the same label. Empty range means start from
    the beginning, and a range past buffer size means the end.
  */
  range_type nextRange(range_type range);

  /*
    Each file must have its own source/sink nodes, even though they might share their
    labels. Hence we check for the file number here.
  */
  inline bool sameFrom(size_type i, size_type j) const
  {
    return (this->buffer[i].node.from == this->buffer[j].node.from &&
            this->buffer[i].file == this->buffer[j].file);
  }

  bool sameFrom(range_type range) const;

  inline range_type min_lcp(const LCP& lcp, size_type i, size_type j) const // i < j
  {
    return lcp.min_lcp(this->buffer[i].node, this->buffer[j].node,
      this->buffer[i].label, this->buffer[j].label);
  }

  inline range_type max_lcp(const LCP& lcp, size_type i, size_type j) const // i < j
  {
    return lcp.max_lcp(this->buffer[i].node, this->buffer[j].node,
      this->buffer[i].label, this->buffer[j].label);
  }

  /*
    Clears elements 0 to range.second - 1 from the buffer and updates the range to (0, 0).
    The remaining element at buffer[0] is assumed to contain the previous merged path.
  */
  void clearUntil(range_type& range);

  /*
    Extends the range forward into a maximal range of paths starting from the same node
    and sharing a common prefix that no other path has. Assumes that the input range only
    contains paths starting from the same node. Returns the lcp of the range as a pair
    (a,b), where a is the lcp of the labels and b is the lcp of the first diverging kmers.
    If the range cannot be extended, the returned lcp value may be incorrect.

    If range.first > 0, the previous merged path must be at range.first - 1.
  */
  range_type extendRange(range_type& range, const LCP& lcp);

  /*
    Merges the path nodes into paths[range.second]. The second version also produces
    the list of all unique from nodes different from the one at paths[range.second]
    as (path_id, from_node).
  */
  void mergePathNodes(range_type range, range_type range_lcp, const LCP& lcp);
  void mergePathNodes(range_type range, std::vector<range_type>& from_nodes, size_type path_id);
};

PathGraphMerger::PathGraphMerger(const PathGraph& path_graph) :
  graph(path_graph), files(path_graph.files()), inputs(path_graph.files())
{
  for(size_type file = 0; file < this->files.size(); file++)
  {
    this->graph.open(this->files[file], file);
    this->inputs[file].file = file;
    this->inputs[file].load(this->files);
  }
  this->inputs.heapify();
  this->fill();
}

void
PathGraphMerger::fill(size_type limit)
{
  while(this->buffer.size() < limit)
  {
    if(this->inputs[0].eof()) { return; }
    this->buffer.push_back(this->inputs[0]);
    this->inputs[0].load(this->files);
    this->inputs.down(0);
  }
}

void
PathGraphMerger::close()
{
  for(size_type file = 0; file < this->files.size(); file++)
  {
    this->files[file].close();
  }
}

range_type
PathGraphMerger::nextRange(range_type range)
{
  if(Range::empty(range)) { range.first = range.second = 0; }
  else { range.first = range.second = range.second + 1; }

  this->readItem(range.second + 1);
  while(range.second + 1 < this->buffer.size()
    && !(this->buffer[range.first] < this->buffer[range.second + 1]))
  {
    range.second++;
    this->readItem(range.second + 1);
  }

  return range;
}

bool
PathGraphMerger::sameFrom(range_type range) const
{
  for(size_type i = range.first + 1; i <= range.second; i++)
  {
    if(!(this->sameFrom(range.first, i))) { return false; }
  }
  return true;
}

void
PathGraphMerger::clearUntil(range_type& range)
{
  for(size_type i = 0; i < range.second; i++) { this->buffer.pop_front(); }
  range.first = range.second = 0;
}

// FIXME Disabled because kmer_lcp does not work
range_type
PathGraphMerger::extendRange(range_type& range, const LCP&)
{
/*  range_type lower_bound(0, 0); // Minimum acceptable LCP.
  if(range.first > 0) { lower_bound = lcp.increment(this->max_lcp(lcp, range.first - 1, range.first)); }*/
  range_type range_lcp(this->buffer[range.first].node.order(), 0);

  /*
    Iterate over one range at a time, doing the following tests:
    1. Is the from node still the same? Stop if not.
    2. Is the LCP still high enough? Stop if not.
    3. Is [range.first, next_range.second] an ancestor of range and next_range? Extend if true.
  */
/*  for(range_type next_range = this->nextRange(range);
    !(this->atEnd(next_range)); next_range = this->nextRange(next_range))
  {
    if(!(this->sameFrom(range_type(next_range.first - 1, next_range.first))) || !(this->sameFrom(next_range)))
    {
      break;
    }
    range_type parent_lcp = this->min_lcp(lcp, range.first, next_range.second);
    if(parent_lcp < lower_bound) { break; }
    this->readItem(next_range.second + 1);
    if(!(this->atEnd(next_range.second + 1)))
    {
      range_type border_lcp = this->max_lcp(lcp, next_range.second, next_range.second + 1);
      if(border_lcp >= parent_lcp) { continue; }
    }
    range.second = next_range.second; range_lcp = parent_lcp;
  }*/

  return range_lcp;
}

void
PathGraphMerger::mergePathNodes(range_type range, range_type range_lcp, const LCP&)
{
  PathNode& merged = this->buffer[range.second].node;
  if(Range::length(range) == 1)
  {
    merged.makeSorted();
    return;
  }

  size_type order = range_lcp.first;
  if(range_lcp.second > 0)
  {
    LCP::rank_range ranks(this->buffer[range.first].node.firstLabel(order, this->buffer[range.first].label),
      this->buffer[range.second].node.lastLabel(order, this->buffer[range.second].label));
// FIXME disabled because kmer lcps don't work
//    ranks = lcp.extendRange(ranks, range_lcp.second);
    this->buffer[range.second].label[order] = ranks.first;
    this->buffer[range.second].label[order + 1] = ranks.second;
    order++;
  }

  merged.makeSorted(); merged.setOrder(order); merged.setLCP(range_lcp.first);
  for(size_type i = range.first; i < range.second; i++)
  {
    merged.addPredecessors(this->buffer[i].node);
  }
}

void
PathGraphMerger::mergePathNodes(range_type range, std::vector<range_type>& from_nodes, size_type path_id)
{
  from_nodes.clear();

  PathNode& merged = this->buffer[range.second].node;
  for(size_type i = range.first; i < range.second; i++)
  {
    merged.addPredecessors(this->buffer[i].node);
    node_type from = this->buffer[i].node.from;
    if(from != merged.from) { from_nodes.push_back(range_type(path_id, from)); }
  }

  removeDuplicates(from_nodes, false);
}

//------------------------------------------------------------------------------

const std::string PathGraph::PREFIX = ".gcsa";

PathGraph::PathGraph(const InputGraph& source, sdsl::sd_vector<>& key_exists)
{
  this->path_count = 0; this->rank_count = 0;
  this->order = source.k();
  this->unique = UNKNOWN; this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;

  sdsl::sd_vector<>::rank_1_type key_rank(&key_exists);
  for(size_type file = 0; file < source.files(); file++)
  {
    std::string temp_file = tempFile(PREFIX);
    this->filenames.push_back(temp_file);
    this->sizes.push_back(source.sizes[file]); this->path_count += source.sizes[file];
    this->rank_counts.push_back(2 * source.sizes[file]); this->rank_count += 2 * source.sizes[file];

    std::ofstream out(temp_file.c_str(), std::ios_base::binary);
    if(!out)
    {
      std::cerr << "PathGraph::PathGraph(): Cannot open output file " << temp_file << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // Read KMers, sort them, and convert them to PathNodes.
    std::vector<KMer> kmers;
    source.read(kmers, file);
    parallelQuickSort(kmers.begin(), kmers.end());
    std::vector<PathNode::rank_type> temp_labels = PathNode::dummyRankVector();
    for(size_type i = 0; i < kmers.size(); i++)
    {
      kmers[i].key = Key::replace(kmers[i].key, key_rank(Key::label(kmers[i].key)));
      PathNode temp(kmers[i], temp_labels);
      temp.serialize(out, temp_labels);
      temp_labels.resize(0);
    }
    out.close();
  }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "PathGraph::PathGraph(): " << this->size() << " paths with "
            << this->ranks() << " ranks" << std::endl;
  std::cerr << "PathGraph::PathGraph(): " << inGigabytes(this->bytes()) << " GB in "
            << this->files() << " file(s)" << std::endl;
#endif
}

PathGraph::PathGraph(size_type file_count, size_type path_order) :
  filenames(file_count), sizes(file_count, 0), rank_counts(file_count, 0),
  path_count(0), rank_count(0), order(path_order),
  unique(0), unsorted(0), nondeterministic(0)
{
  for(size_type file = 0; file < this->files(); file++)
  {
    this->filenames[file] = tempFile(PREFIX);
  }
}

PathGraph::~PathGraph()
{
  this->clear();
}

void
PathGraph::clear()
{
  for(size_type file = 0; file < this->files(); file++)
  {
    remove(this->filenames[file].c_str());
  }
  this->filenames.clear();
  this->sizes.clear();
  this->rank_counts.clear();

  this->path_count = 0; this->rank_count = 0;
  this->order = 0;
  this->unique = UNKNOWN; this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;
}

void
PathGraph::swap(PathGraph& another)
{
  this->filenames.swap(another.filenames);
  this->sizes.swap(another.sizes);
  this->rank_counts.swap(another.rank_counts);

  std::swap(this->path_count, another.path_count);
  std::swap(this->rank_count, another.rank_count);
  std::swap(this->order, another.order);

  std::swap(this->unique, another.unique);
  std::swap(this->unsorted, another.unsorted);
  std::swap(this->nondeterministic, another.nondeterministic);
}

void
PathGraph::open(std::ifstream& input, size_type file) const
{
  if(file >= this->files())
  {
    std::cerr << "PathGraph::open(): Invalid file number: " << file << std::endl;
    std::exit(EXIT_FAILURE);
  }

  input.open(this->filenames[file].c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "PathGraph::open(): Cannot open graph file " << this->filenames[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//------------------------------------------------------------------------------

void
PathGraph::prune(const LCP& lcp, size_type size_limit)
{
#ifdef VERBOSE_STATUS_INFO
  size_type old_path_count = this->size();
#endif

  PathGraphMerger merger(*this);
  PathGraphBuilder builder(this->files(), this->k(), size_limit);
  for(range_type range = merger.firstRange(); !(merger.atEnd(range)); range = merger.nextRange(range))
  {
    if(merger.sameFrom(range))
    {
      range_type range_lcp = merger.extendRange(range, lcp);
      merger.mergePathNodes(range, range_lcp, lcp);
      builder.write(merger.buffer[range.second]);
      builder.graph.unique++;
    }
    else
    {
      for(size_type i = range.first; i <= range.second; i++)
      {
        if(merger.buffer[i].node.sorted()) { builder.graph.nondeterministic++; }
        else { builder.graph.unsorted++; }
        builder.write(merger.buffer[i]);
      }
    }
    merger.clearUntil(range);
  }
  merger.close(); builder.close();
  this->clear(); this->swap(builder.graph);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "PathGraph::prune(): " << old_path_count << " -> " << this->size() << " paths" << std::endl;
  std::cerr << "PathGraph::prune(): " << this->unique << " unique, " << this->unsorted << " unsorted, "
            << this->nondeterministic << " nondeterministic paths" << std::endl;
  std::cerr << "PathGraph::prune(): " << inGigabytes(this->bytes()) << " GB in "
            << this->files() << " file(s)" << std::endl;
#endif
}

//------------------------------------------------------------------------------

void
PathGraph::extend(size_type size_limit)
{
#ifdef VERBOSE_STATUS_INFO
  size_type old_path_count = this->size();
#endif

  PathGraphBuilder builder(this->files(), 2 * this->k(), size_limit);
  for(size_type file = 0; file < this->files(); file++)
  {
    // Read the current file.
    std::vector<PathNode> paths;
    std::vector<PathNode::rank_type> labels;
    this->read(paths, labels, file);

    // Initialization.
    PathFromComparator from_c;  // Sort the paths by from.
    parallelQuickSort(paths.begin(), paths.end(), from_c);
    ValueIndex<PathNode, FromGetter> from_index(paths);
    size_type threads = omp_get_max_threads();

    // Create the next generation.
    std::vector<PathNode> temp_nodes[threads];
    std::vector<PathNode::rank_type> temp_labels[threads];
    #pragma omp parallel for schedule(static)
    for(size_type i = 0; i < paths.size(); i++)
    {
      size_type thread = omp_get_thread_num();
      if(paths[i].sorted())
      {
        temp_nodes[thread].push_back(PathNode(paths[i], labels, temp_labels[thread]));
      }
      else
      {
        size_type first = from_index.find(paths[i].to);
        for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++)
        {
          temp_nodes[thread].push_back(PathNode(paths[i], paths[j], labels, temp_labels[thread]));
        }
      }
      if(temp_nodes[thread].size() >= PathGraphBuilder::WRITE_BUFFER_SIZE)
      {
        builder.write(temp_nodes[thread], temp_labels[thread], file);
      }
    }
    for(size_type thread = 0; thread < threads; thread++)
    {
      builder.write(temp_nodes[thread], temp_labels[thread], file);
    }
    sdsl::util::clear(paths); sdsl::util::clear(labels);

#ifdef VERBOSE_STATUS_INFO
    std::cerr << "PathGraph::extend(): File " << file << ": Created " << builder.graph.sizes[file]
              << " order-" << builder.graph.k() << " paths" << std::endl;
#endif

    // Sort the next generation.
    builder.sort(file);
  }
  builder.close();
  this->clear(); this->swap(builder.graph);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "PathGraph::extend(): " << old_path_count << " -> " << this->size() << " paths ("
            << this->ranks() << " ranks)" << std::endl;
  std::cerr << "PathGraph::extend(): " << inGigabytes(this->bytes()) << " GB in "
            << this->files() << " file(s)" << std::endl;
#endif
}

//------------------------------------------------------------------------------

void
PathGraph::read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file) const
{
  sdsl::util::clear(paths); sdsl::util::clear(labels);

  std::ifstream input; this->open(input, file);
  paths.reserve(this->sizes[file]); labels.reserve(this->rank_counts[file]);
  for(size_type i = 0; i < this->sizes[file]; i++) { paths.push_back(PathNode(input, labels)); }
  input.close();

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "PathGraph::read(): File " << file << ": Read " << paths.size()
            << " order-" << this->k()<< " paths" << std::endl;
#endif
}

//------------------------------------------------------------------------------

const std::string MergedGraph::PREFIX = ".gcsa";

MergedGraph::MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper) :
  path_name(tempFile(PREFIX)), rank_name(tempFile(PREFIX)), from_name(tempFile(PREFIX)),
  path_d(NO_FILE), rank_d(NO_FILE), from_d(NO_FILE),
  paths(0), labels(0), from_nodes(0),
  path_count(0), rank_count(0), from_count(0),
  order(source.k()),
  next(mapper.alpha.sigma + 1, 0), next_from(mapper.alpha.sigma + 1, 0)
{
  std::ofstream path_file(this->path_name.c_str(), std::ios_base::binary);
  if(!path_file)
  {
    std::cerr << "MergedGraph::MergedGraph(): Cannot open output file " << this->path_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::ofstream rank_file(this->rank_name.c_str(), std::ios_base::binary);
  if(!rank_file)
  {
    std::cerr << "MergedGraph::MergedGraph(): Cannot open output file " << this->rank_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::ofstream from_file(this->from_name.c_str(), std::ios_base::binary);
  if(!from_file)
  {
    std::cerr << "MergedGraph::MergedGraph(): Cannot open output file " << this->from_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  /*
     Initialize next[comp] to be the the rank of the first kmer starting with
     the corresponding character. Later, next[comp] is transformed into the rank
     of the first path with firstLabel(0) >= next[comp].
  */
  for(size_type comp = 0; comp < mapper.alpha.sigma; comp++)
  {
    this->next[comp] = mapper.charRange(comp).first;
  }
  this->next[mapper.alpha.sigma] = ~(size_type)0;
  this->next_from[mapper.alpha.sigma] = ~(size_type)0;

  PathGraphMerger merger(source);
  std::vector<range_type> curr_from;
  size_type curr_comp = 0;  // Used to transform next.
  for(range_type range = merger.firstRange(); !(merger.atEnd(range)); range = merger.nextRange(range))
  {
    merger.mergePathNodes(range, curr_from, this->path_count);
    merger.buffer[range.second].serialize(path_file, rank_file, this->rank_count);
    if(curr_from.size() > 0)
    {
      from_file.write((char*)(curr_from.data()), curr_from.size() * sizeof(range_type));
    }
    while(merger.buffer[range.second].firstLabel(0) >= this->next[curr_comp])
    {
      this->next[curr_comp] = this->path_count;
      this->next_from[curr_comp] = this->from_count;
      curr_comp++;
    }
    this->path_count++;
    this->rank_count += merger.buffer[range.second].node.ranks();
    this->from_count += curr_from.size();
    merger.clearUntil(range);
  }
  merger.close();
  path_file.close(); rank_file.close(); from_file.close();

  this->paths = (PathNode*)(this->map(this->path_name, "path", this->path_d, this->path_bytes()));
  this->labels = (PathNode::rank_type*)(this->map(this->rank_name, "rank", this->rank_d, this->rank_bytes()));
  this->from_nodes = (range_type*)(this->map(this->from_name, "from", this->from_d, this->from_bytes()));

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "MergedGraph::MergedGraph(): " << this->size() << " paths with "
            << this->ranks() << " ranks and " << this->extra() << " additional from nodes" << std::endl;
  std::cerr << "MergedGraph::MergedGraph(): " << inGigabytes(this->bytes()) << " GB" << std::endl;
#endif
}

MergedGraph::~MergedGraph()
{
  this->clear();
}

void
MergedGraph::clear()
{
  // Unmap.
  if(this->paths != 0) { munmap(this->paths, this->path_bytes()); this->paths = 0; }
  if(this->labels != 0) { munmap(this->labels, this->rank_bytes()); this->labels = 0; }
  if(this->from_nodes != 0) { munmap(this->from_nodes, this->from_bytes()); this->from_nodes = 0; }

  // Close.
  if(this->path_d != NO_FILE) { close(this->path_d); this->path_d = NO_FILE; }
  if(this->rank_d != NO_FILE) { close(this->rank_d); this->rank_d = NO_FILE; }
  if(this->from_d != NO_FILE) { close(this->from_d); this->from_d = NO_FILE; }

  // Remove.
  remove(this->path_name.c_str()); this->path_name = "";
  remove(this->rank_name.c_str()); this->rank_name = "";
  remove(this->from_name.c_str()); this->from_name = "";

  this->path_count = 0; this->rank_count = 0; this->from_count = 0;
  this->order = 0;

  for(size_type i = 0; i < this->next.size(); i++) { this->next[i] = 0; }
}

void*
MergedGraph::map(const std::string& filename, const std::string& file_type, int& fd, size_type n)
{
  if((fd = open(filename.c_str(), O_RDONLY)) == NO_FILE)
  {
    std::cerr << "MergedGraph::map(): Cannot open " << file_type << " file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }

  void* pointer = mmap(0, n, PROT_READ, MAP_FILE | MAP_SHARED, fd, 0);
  if(pointer == 0)
  {
    std::cerr << "MergedGraph::map(): Cannot memory map " << file_type << " file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  madvise(pointer, n, MADV_SEQUENTIAL);

  return pointer;
}

//------------------------------------------------------------------------------

} // namespace gcsa
