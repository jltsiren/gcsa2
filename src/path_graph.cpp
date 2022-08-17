/*
  Copyright (c) 2018, 2019 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#include <gcsa/path_graph.h>

#include <deque>

namespace gcsa
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type PathLabel::LABEL_LENGTH;
constexpr PathLabel::rank_type PathLabel::NO_RANK;

constexpr size_type PathNode::LABEL_LENGTH;

constexpr size_type PathGraph::UNKNOWN;

constexpr size_type MergedGraph::UNKNOWN;

//------------------------------------------------------------------------------

// Other class variables.

const std::string PathGraph::PREFIX = "gcsa";

const std::string MergedGraph::PREFIX = "gcsa";

//------------------------------------------------------------------------------

PathNode::PathNode(const KMer& kmer, WriteBuffer<PathNode::rank_type>& labels)
{
  this->from = kmer.from; this->to = kmer.to;
  this->fields = 0;

  if(kmer.sorted()) { this->makeSorted(); }
  this->setPredecessors(Key::predecessors(kmer.key));
  this->setOrder(1); this->setLCP(1);

  this->setPointer(labels.size());
  labels.push_back(Key::label(kmer.key));
  labels.push_back(0);  // Dummy value; the last label is not in use.
}

PathNode::PathNode(const PathNode& source,
    const std::vector<PathNode::rank_type>& old_labels, std::vector<PathNode::rank_type>& new_labels)
{
  this->from = source.from; this->to = source.to;
  this->fields = source.fields;

  this->setPointer(new_labels.size());
  for(size_type i = 0, j = source.pointer(); i < source.ranks(); i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
}

PathNode::PathNode(const PathNode& left, const PathNode& right,
    const std::vector<PathNode::rank_type>& old_labels, std::vector<PathNode::rank_type>& new_labels)
{
  this->from = left.from; this->to = right.to;
  if(right.sorted()) { this->makeSorted(); }

  this->fields = 0;
  this->setPredecessors(left.predecessors());

  size_type left_order = left.order();
  size_type new_order = left_order + right.order();
  this->setOrder(new_order); this->setLCP(left_order + right.lcp());

  this->setPointer(new_labels.size());
  for(size_type i = 0, j = left.pointer(); i < left_order; i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
  for(size_type i = 0, j = right.pointer(); i < right.ranks(); i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
}

PathNode::PathNode()
{
  this->from = 0; this->to = 0;
  this->fields = 0;
}

PathNode::PathNode(const PathNode& source)
{
  this->copy(source);
}

PathNode::PathNode(PathNode&& source)
{
  *this = std::move(source);
}

PathNode::~PathNode()
{
}

PathNode&
PathNode::operator= (const PathNode& source)
{
  if(&source != this)
  {
    this->copy(source);
  }
  return *this;
}

void
PathNode::copy(const PathNode& source)
{
  this->from = source.from; this->to = source.to;
  this->fields = source.fields;
}

PathNode&
PathNode::operator= (PathNode&& source)
{
  if(&source != this)
  {
    this->from = std::move(source.from);
    this->to = std::move(source.to);
    this->fields = std::move(source.fields);
  }
  return *this;
}

void
PathNode::print(std::ostream& out, const std::vector<PathNode::rank_type>& labels) const
{
  this->print(out, labels.data());
}

void
PathNode::print(std::ostream& out, const rank_type* labels) const
{
  out << "(" << Node::decode(this->from) << " -> " << Node::decode(this->to);
  out << "; o" << this->order();
  out << "; l" << this->lcp();
  for(size_type i = 0; i < this->order(); i++)
  {
    out << (i == 0 ? "; [" : ", ") << this->firstLabel(i, labels);
  }
  for(size_type i = 0; i < this->order(); i++)
  {
    out << (i == 0 ? " to " : ", ") << this->lastLabel(i, labels);
  }
  out << "])";
}

//------------------------------------------------------------------------------

LCP::LCP()
{
}

LCP::LCP(const std::vector<key_type>& keys, size_type _kmer_length)
{
  this->kmer_length = _kmer_length;
  this->total_keys = keys.size();

  sdsl::int_vector<8> buffer(keys.size(), 0);
  for(size_type i = 1; i < keys.size(); i++)
  {
    buffer[i] = Key::lcp(keys[i - 1], keys[i], this->kmer_length);
  }
  directConstruct(this->kmer_lcp, buffer);
}

range_type
LCP::min_lcp(const PathNode& a, const PathNode& b, const std::vector<LCP::rank_type>& labels) const
{
  return this->min_lcp(a, b, labels.data(), labels.data());
}

range_type
LCP::max_lcp(const PathNode& a, const PathNode& b, const std::vector<LCP::rank_type>& labels) const
{
  return this->max_lcp(a, b, labels.data(), labels.data());
}

range_type
LCP::min_lcp(const PathNode& a, const PathNode& b,
  const LCP::rank_type* a_labels, const LCP::rank_type* b_labels) const
{
  size_type order = std::min(a.order(), b.order());
  range_type lcp(0, 0);
  while(lcp.first < order && a.firstLabel(lcp.first, a_labels) == b.lastLabel(lcp.first, b_labels))
  {
    lcp.first++;
  }
  if(lcp.first < order)
  {
    size_type left = a.firstLabel(lcp.first, a_labels) + 1;
    size_type right = std::min((size_type)(b.lastLabel(lcp.first, b_labels)), this->total_keys - 1);
    lcp.second = sdsl::quantile_freq(this->kmer_lcp, left, right, 0).first;
  }
  return lcp;
}

range_type
LCP::max_lcp(const PathNode& a, const PathNode& b,
  const LCP::rank_type* a_labels, const LCP::rank_type* b_labels) const
{
  size_type order = std::min(a.order(), b.order());
  range_type lcp(0, 0);
  while(lcp.first < order && a.lastLabel(lcp.first, a_labels) == b.firstLabel(lcp.first, b_labels))
  {
    lcp.first++;
  }
  if(lcp.first < order)
  {
    size_type left = a.lastLabel(lcp.first, a_labels) + 1;
    size_type right = b.firstLabel(lcp.first, b_labels);
    lcp.second = sdsl::quantile_freq(this->kmer_lcp, left, right, 0).first;
  }
  return lcp;
}

void
LCP::swap(LCP& another)
{
  std::swap(this->kmer_length, another.kmer_length);
  std::swap(this->total_keys, another.total_keys);
  this->kmer_lcp.swap(another.kmer_lcp);
}

//------------------------------------------------------------------------------

/*
  This structure combines a PathNode and its label. It also stores the identifier of
  its source file.
*/

struct PriorityNode
{
  typedef PathLabel::rank_type rank_type;

  constexpr static size_type LABEL_LENGTH = PathLabel::LABEL_LENGTH;
  constexpr static rank_type NO_RANK = PathLabel::NO_RANK;

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

  inline rank_type firstLabel(size_type i) const { return this->label[i]; }

  inline size_type bytes() const { return this->node.bytes(); }
};

constexpr size_type PriorityNode::LABEL_LENGTH;
constexpr PriorityNode::rank_type PriorityNode::NO_RANK;

//------------------------------------------------------------------------------

/*
  This structure builds a PathGraph. It expects a stream of PriorityNodes. When all paths
  for a certain file have been written, call sort() for that file.
*/

struct PathGraphBuilder
{
  PathGraph graph;
  std::vector<WriteBuffer<PathNode>> path_files;
  std::vector<WriteBuffer<PathNode::rank_type>> rank_files;
  size_type limit;  // Bytes of disk space.

  constexpr static size_type WRITE_BUFFER_SIZE = MEGABYTE;  // PathNodes per thread.

  PathGraphBuilder(size_type file_count, size_type path_order, size_type step, size_type size_limit);
  void close();

  /*
    The file number is assumed to be valid.
    The first call is not thread safe, while the bulk write() is.
  */
  void write(PriorityNode& path);
  void write(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file);

  void sort(size_type file);
};

constexpr size_type PathGraphBuilder::WRITE_BUFFER_SIZE;

PathGraphBuilder::PathGraphBuilder(size_type file_count, size_type path_order, size_type step, size_type size_limit) :
  graph(file_count, path_order, step),
  path_files(file_count), rank_files(file_count),
  limit(size_limit)
{
  for(size_type file = 0; file < file_count; file++)
  {
    this->path_files[file].open(this->graph.path_names[file]);
    this->rank_files[file].open(this->graph.rank_names[file]);
  }
}

void
PathGraphBuilder::close()
{
  for(size_type file = 0; file < this->path_files.size(); file++)
  {
    this->path_files[file].close();
    this->rank_files[file].close();
  }
}

inline void
writePath(PathNode& path, const PathNode::rank_type* labels,
  WriteBuffer<PathNode>& path_file, WriteBuffer<PathNode::rank_type>& rank_file)
{
  size_type old_ptr = path.pointer();
  size_type limit = old_ptr + path.ranks();

  path.setPointer(rank_file.size());
  path_file.push_back(path);

  for(size_type i = old_ptr; i < limit; i++) { rank_file.push_back(labels[i]); }
}

void
PathGraphBuilder::write(PriorityNode& path)
{
  if(this->graph.bytes() + path.bytes() > this->limit)
  {
    std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  writePath(path.node, path.label, this->path_files[path.file], this->rank_files[path.file]);

  this->graph.path_counts[path.file]++; this->graph.path_count++;
  this->graph.rank_counts[path.file] += path.node.ranks();
  this->graph.rank_count += path.node.ranks();
}

void
PathGraphBuilder::write(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file)
{
  size_type bytes_required = paths.size() * sizeof(PathNode) + labels.size() * sizeof(PathNode::rank_type);
  #pragma omp critical
  {
    if(bytes_required + this->graph.bytes() > this->limit)
    {
      std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for(size_type i = 0; i < paths.size(); i++)
    {
      writePath(paths[i], labels.data(), this->path_files[file], this->rank_files[file]);
    }
    this->graph.path_counts[file] += paths.size(); this->graph.path_count += paths.size();
    this->graph.rank_counts[file] += labels.size(); this->graph.rank_count += labels.size();
  }
  paths.clear(); labels.clear();
}

void
PathGraphBuilder::sort(size_type file)
{
  this->path_files[file].close();
  this->rank_files[file].close();

  std::vector<PathNode> paths;
  std::vector<PathNode::rank_type> labels;
  this->graph.read(paths, labels, file);

  PathFirstComparator first_c(labels);
  parallelQuickSort(paths.begin(), paths.end(), first_c);

  this->path_files[file].open(this->graph.path_names[file]);
  this->rank_files[file].open(this->graph.rank_names[file]);
  for(size_type i = 0; i < paths.size(); i++)
  {
    writePath(paths[i], labels.data(), this->path_files[file], this->rank_files[file]);
  }

  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "PathGraphBuilder::sort(): File " << file << ": Sorted "
              << this->graph.path_counts[file] << " paths" << std::endl;
  }
}

//------------------------------------------------------------------------------

struct PathGraphMerger;

struct PathRange
{
  size_type  from, to;
  range_type left_lcp, range_lcp, right_lcp;

  inline range_type range() const { return range_type(this->from, this->to); }
  inline size_type length() const { return this->to + 1 - this->from; }

  PathRange(size_type start, size_type stop, range_type _left_lcp, PathGraphMerger& merger);
};

/*
  This structure reads a buffered stream of PriorityNodes in sorted order and outputs a
  stream of ranges of PriorityNodes with the same label.
*/

struct PathGraphMerger
{
  const PathGraph&                              graph;
  const LCP&                                    lcp;

  // Buffers.
  std::deque<PathRange>                         ranges;
  BufferWindow<PriorityNode>                    buffer;

  // Priority queue.
  std::vector<ReadBuffer<PathNode>>             path_files;
  std::vector<ReadBuffer<PathNode::rank_type>>  rank_files;
  std::vector<size_type>                        offsets;
  PriorityQueue<PriorityNode>                   inputs;

  PathGraphMerger(const PathGraph& path_graph, const LCP& kmer_lcp);
  void close();

  inline size_type size() const { return this->graph.size(); }

  /*
    Iterates through ranges of paths with the same label.
  */
  range_type first();
  range_type next();
  inline bool atEnd(range_type range) const { return (range.first >= this->size()); }

  inline range_type range_lcp(size_type i, size_type j) const // i < j
  {
    return this->lcp.min_lcp(this->buffer[i].node, this->buffer[j].node,
      this->buffer[i].label, this->buffer[j].label);
  }

  inline range_type border_lcp(size_type i, size_type j) const // i < j
  {
    return this->lcp.max_lcp(this->buffer[i].node, this->buffer[j].node,
      this->buffer[i].label, this->buffer[j].label);
  }

  /*
    Extends the front range forward into a maximal range of paths such that

    (a) comp(next_range) returns true for each additional next_range; and
    (b) the paths share a common prefix no other path has.
  */
  template<class FromComparator>
  range_type extendRange(FromComparator& comp);

  /*
    Merges the path nodes into paths[range.second] of the front range.
  */
  void mergePathNodes();

  // Find the rightmost path with the same label.
  size_type rangeEnd(size_type start);

  // Add the next PriorityNode to buffer.
  void bufferNext();

  // Read the next PriorityNode from the file.
  void read(PriorityNode& path);
};

PathGraphMerger::PathGraphMerger(const PathGraph& path_graph, const LCP& kmer_lcp) :
  graph(path_graph), lcp(kmer_lcp),
  path_files(path_graph.files()), rank_files(path_graph.files()),
  offsets(path_graph.files()), inputs(path_graph.files())
{
  for(size_type file = 0; file < path_graph.files(); file++)
  {
    this->path_files[file].open(path_graph.path_names[file]);
    this->rank_files[file].open(path_graph.rank_names[file]);
    this->offsets[file] = 0;
    this->inputs[file].file = file; this->read(this->inputs[file]);
  }
  this->inputs.heapify();
}

void
PathGraphMerger::close()
{
  sdsl::util::clear(this->ranges);
  sdsl::util::clear(this->buffer);

  for(size_type file = 0; file < this->graph.files(); file++)
  {
    this->path_files[file].close();
    this->rank_files[file].close();
  }
  this->path_files.clear();
  this->rank_files.clear();
  this->offsets.clear();
  this->inputs.clear();
}

range_type
PathGraphMerger::first()
{
  if(this->buffer.offset > 0)
  {
    std::cerr << "PathGraphMerger::first(): Cannot seek back in the buffer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->ranges.clear();

  this->ranges.push_back(PathRange(0, this->rangeEnd(0), range_type(0, 0), *this));
  return this->ranges.front().range();
}

range_type
PathGraphMerger::next()
{
  PathRange temp = this->ranges.front(); this->ranges.pop_front();
  if(this->ranges.empty())
  {
    this->ranges.push_back(
      PathRange(temp.to + 1, this->rangeEnd(temp.to + 1), temp.right_lcp, *this));
  }
  this->buffer.seek(this->ranges.front().from);
  return this->ranges.front().range();
}

template<class FromComparator>
range_type
PathGraphMerger::extendRange(FromComparator& comp)
{
  PathRange range = this->ranges.front();
  size_type curr = 1;
  range_type parent_lcp = range.range_lcp;

  while(true)
  {
    // Find the next range.
    if(curr >= this->ranges.size())
    {
      const PathRange& temp = this->ranges.back();
      this->ranges.push_back(PathRange(temp.to + 1, this->rangeEnd(temp.to + 1), temp.right_lcp, *this));
    }
    const PathRange& next_range = this->ranges[curr];
    if(next_range.from >= this->size()) { break; }

    // Is this a suffix tree node with the same start nodes as range and lcp > range.left_lcp?
    if(!(comp(next_range.range()))) { break; }
    parent_lcp = std::min(parent_lcp, next_range.left_lcp);
    if(parent_lcp <= range.left_lcp) { break; }
    if(next_range.right_lcp >= parent_lcp) { curr++; continue; }

    // Replace range with [range.from, next_range.to].
    range.to = next_range.to;
    range.range_lcp = parent_lcp; range.right_lcp = next_range.right_lcp;
    for(size_type i = 0; i <= curr; i++) { this->ranges.pop_front(); }
    this->ranges.push_front(range); curr = 1;
  }

  return range.range();
}

void
PathGraphMerger::mergePathNodes()
{
  PathRange& range = this->ranges.front();
  if(range.length() == 1)
  {
    this->buffer[range.to].node.makeSorted();
    return;
  }

  size_type order = range.range_lcp.first;
  if(range.range_lcp.second > 0)
  {
    range_type ranks(this->buffer[range.from].node.firstLabel(order, this->buffer[range.from].label),
      this->buffer[range.to].node.lastLabel(order, this->buffer[range.to].label));
    this->buffer[range.to].label[order] = ranks.first;
    this->buffer[range.to].label[order + 1] = ranks.second;
    order++;
  }

  this->buffer[range.to].node.makeSorted();
  this->buffer[range.to].node.setOrder(order);
  this->buffer[range.to].node.setLCP(range.range_lcp.first);
  for(size_type i = range.from; i < range.to; i++)
  {
    this->buffer[range.to].node.addPredecessors(this->buffer[i].node);
  }
}

size_type
PathGraphMerger::rangeEnd(size_type start)
{
  if(!(this->buffer.buffered(start))) { this->bufferNext(); }

  size_type stop = start;
  while(stop + 1 < this->size())
  {
    if(!(this->buffer.buffered(stop + 1))) { this->bufferNext(); }
    if(this->buffer[start] < this->buffer[stop + 1]) { break; }
    stop++;
  }
  return stop;
}

void
PathGraphMerger::bufferNext()
{
  this->buffer.push_back(this->inputs[0]);  // Add to the buffer.
  this->read(this->inputs[0]);  // Read the next.
  this->inputs.down(0); // Restore heap order.
}

void
PathGraphMerger::read(PriorityNode& path)
{
  if(this->offsets[path.file] >= this->graph.path_counts[path.file])
  {
    path.node.setOrder(1);
    path.node.setLCP(1);
    path.label[0] = PriorityNode::NO_RANK;
  }
  else
  {
    this->path_files[path.file].seek(this->offsets[path.file]);
    path.node = this->path_files[path.file][this->offsets[path.file]];
    this->rank_files[path.file].seek(path.node.pointer());
    for(size_type i = 0; i < path.node.ranks(); i++)
    {
      path.label[i] = this->rank_files[path.file][path.node.pointer() + i];
    }
    path.node.setPointer(0);  // Label is now stored in the PriorityNode.
    this->offsets[path.file]++;
  }
}

PathRange::PathRange(size_type start, size_type stop, range_type _left_lcp, PathGraphMerger& merger) :
  from(start), to(stop),
  left_lcp(_left_lcp),
  range_lcp(0, 0),
  right_lcp(0, 0)
{
  if(stop < merger.size()) { this->range_lcp = merger.range_lcp(start, stop); }
  if(stop + 1 < merger.size()) { this->right_lcp = merger.border_lcp(stop, stop + 1); }
}

//------------------------------------------------------------------------------

PathGraph::PathGraph(const InputGraph& source, sdsl::int_vector<0>& distinct_labels)
{
  this->path_count = 0; this->rank_count = 0; this->range_count = 0;
  this->order = source.k(); this->doubling_steps = 0;
  this->unique = UNKNOWN; this->redundant = UNKNOWN;
  this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;
  this->delete_files = true;

  for(size_type file = 0; file < source.files(); file++)
  {
    std::string path_name = TempFile::getName(PREFIX);
    this->path_names.push_back(path_name);
    std::string rank_name = TempFile::getName(PREFIX);
    this->rank_names.push_back(rank_name);

    // Read KMers, sort them, and convert the keys labels to the ranks of those labels.
    std::vector<KMer> kmers;
    source.read(kmers, file);
    parallelQuickSort(kmers.begin(), kmers.end());
    size_type current_rank = 0;
    for(size_type i = 0; i < kmers.size(); i++)
    {
      while(Key::label(kmers[i].key) > distinct_labels[current_rank]) { current_rank++; }
      kmers[i].key = Key::replace(kmers[i].key, current_rank);
    }

    // Convert the KMers to PathNodes.
    WriteBuffer<PathNode> path_buffer(path_name);
    WriteBuffer<PathNode::rank_type> rank_buffer(rank_name);
    for(size_type i = 0; i < kmers.size(); i++)
    {
      path_buffer.push_back(PathNode(kmers[i], rank_buffer));
    }
    this->path_counts.push_back(path_buffer.size()); this->path_count += path_buffer.size();
    this->rank_counts.push_back(rank_buffer.size()); this->rank_count += rank_buffer.size();
    path_buffer.close(); rank_buffer.close();
  }

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::PathGraph(): " << this->size() << " paths with "
              << this->ranks() << " ranks" << std::endl;
    std::cerr << "PathGraph::PathGraph(): " << inGigabytes(this->bytes()) << " GB in "
              << this->files() << " file(s)" << std::endl;
  }
}

PathGraph::PathGraph(size_type file_count, size_type path_order, size_type steps) :
  path_names(file_count), rank_names(file_count), path_counts(file_count, 0), rank_counts(file_count, 0),
  path_count(0), rank_count(0), range_count(0), order(path_order), doubling_steps(steps),
  unique(0), redundant(0), unsorted(0), nondeterministic(0),
  delete_files(true)
{
  for(size_type file = 0; file < this->files(); file++)
  {
    this->path_names[file] = TempFile::getName(PREFIX);
    this->rank_names[file] = TempFile::getName(PREFIX);
  }
}

PathGraph::PathGraph(const std::string& path_name, const std::string& rank_name)
{
  this->path_count = 0; this->rank_count = 0; this->range_count = 0;
  this->order = 0; this->doubling_steps = 0;
  this->unique = 0; this->redundant = 0;
  this->unsorted = 0; this->nondeterministic = 0;
  this->delete_files = false;

  this->path_names.push_back(path_name);
  std::ifstream path_file(path_name, std::ios_base::binary);
  if(!path_file)
  {
    std::cerr << "PathGraph::PathGraph(): Cannot open path file " << path_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->path_count = fileSize(path_file) / sizeof(PathNode);
  this->path_counts.push_back(this->path_count);
  path_file.close();

  this->rank_names.push_back(rank_name);
  std::ifstream rank_file(rank_name, std::ios_base::binary);
  if(!rank_file)
  {
    std::cerr << "PathGraph::PathGraph(): Cannot open rank file " << rank_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->rank_count = fileSize(rank_file) / sizeof(PathNode::rank_type);
  this->rank_counts.push_back(this->rank_count);
  rank_file.close();
}

PathGraph::~PathGraph()
{
  this->clear();
}

void
PathGraph::clear()
{
  if(this->delete_files)
  {
    for(size_type file = 0; file < this->files(); file++)
    {
      TempFile::remove(this->path_names[file]);
      TempFile::remove(this->rank_names[file]);
    }
  }
  this->path_names.clear();
  this->rank_names.clear();
  this->path_counts.clear();
  this->rank_counts.clear();

  this->path_count = 0; this->rank_count = 0;
  this->order = 0;
  this->unique = UNKNOWN; this->redundant = UNKNOWN;
  this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;
}

void
PathGraph::swap(PathGraph& another)
{
  this->path_names.swap(another.path_names);
  this->rank_names.swap(another.rank_names);
  this->path_counts.swap(another.path_counts);
  this->rank_counts.swap(another.rank_counts);

  std::swap(this->path_count, another.path_count);
  std::swap(this->rank_count, another.rank_count);
  std::swap(this->range_count, another.range_count);
  std::swap(this->order, another.order);
  std::swap(this->doubling_steps, another.doubling_steps);

  std::swap(this->unique, another.unique);
  std::swap(this->redundant, another.redundant);
  std::swap(this->unsorted, another.unsorted);
  std::swap(this->nondeterministic, another.nondeterministic);
}

void
PathGraph::open(std::ifstream& path_file, std::ifstream& rank_file, size_type file) const
{
  if(file >= this->files())
  {
    std::cerr << "PathGraph::open(): Invalid file number: " << file << std::endl;
    std::exit(EXIT_FAILURE);
  }

  path_file.open(this->path_names[file].c_str(), std::ios_base::binary);
  if(!path_file)
  {
    std::cerr << "PathGraph::open(): Cannot open path file " << this->path_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }

  rank_file.open(this->rank_names[file].c_str(), std::ios_base::binary);
  if(!rank_file)
  {
    std::cerr << "PathGraph::open(): Cannot open rank file " << this->rank_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//------------------------------------------------------------------------------

struct SameFromFile
{
  const PathGraphMerger& merger;
  node_type              from;
  size_type              file;
  bool                   same_from, same_file;

  SameFromFile(const PathGraphMerger& source, range_type range) :
    merger(source), from(source.buffer[range.first].node.from), file(source.buffer[range.first].file),
    same_from(true), same_file(true)
  {
    for(size_type i = range.first + 1; i <= range.second; i++)
    {
      if(this->merger.buffer[i].node.from != this->from) { this->same_from = false; }
      if(this->merger.buffer[i].file != this->file) { this->same_file = false; }
    }
  }

  inline bool operator() (range_type range) const
  {
    for(size_type i = range.first; i <= range.second; i++)
    {
      if(this->merger.buffer[i].node.from != this->from || this->merger.buffer[i].file != this->file)
      {
        return false;
      }
    }
    return true;
  }
};

void
PathGraph::prune(const LCP& lcp, size_type size_limit)
{
  size_type old_path_count = this->size();

  PathGraphMerger merger(*this, lcp);
  PathGraphBuilder builder(this->files(), this->k(), this->step(), size_limit);
  for(range_type range = merger.first(); !(merger.atEnd(range)); range = merger.next())
  {
    SameFromFile same_from(merger, range);
    if(same_from.same_from)
    {
      if(same_from.same_file)
      {
        range = merger.extendRange(same_from);
        merger.mergePathNodes();
        builder.write(merger.buffer[range.second]);
        builder.graph.unique++;
      }
      else
      {
        // FIXME Later: Write just one path per file.
        for(size_type i = range.first; i <= range.second; i++)
        {
          merger.buffer[i].node.makeSorted();
          builder.write(merger.buffer[i]);
        }
        builder.graph.redundant += Range::length(range);
      }
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
    builder.graph.range_count++;
  }
  merger.close(); builder.close();
  this->clear(); this->swap(builder.graph);

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::prune(): " << old_path_count << " -> " << this->size() << " paths ("
              << this->ranges() << " ranges)" << std::endl;
    std::cerr << "PathGraph::prune(): "
              << this->unique << " unique, "
              << this->redundant << " redundant, "
              << this->unsorted << " unsorted, "
              << this->nondeterministic << " nondeterministic paths" << std::endl;
    std::cerr << "PathGraph::prune(): " << inGigabytes(this->bytes()) << " GB in "
              << this->files() << " file(s)" << std::endl;
  }
}

//------------------------------------------------------------------------------

void
PathGraph::extend(size_type size_limit)
{
  size_type old_path_count = this->size();

  PathGraphBuilder builder(this->files(), 2 * this->k(), this->step() + 1, size_limit);
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

    // Create thread-specific buffers.
    std::vector<std::vector<PathNode>> temp_nodes(threads);
    std::vector<std::vector<PathNode::rank_type>> temp_labels(threads);
    for(size_type thread = 0; thread < threads; thread++)
    {
      temp_nodes[thread].reserve(PathGraphBuilder::WRITE_BUFFER_SIZE);
      temp_labels[thread].reserve(((1 << builder.graph.step()) + 1) * PathGraphBuilder::WRITE_BUFFER_SIZE);
    }

    // Create the next generation.
    size_type chunk_size = getChunkSize(paths.size(), MEGABYTE);
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for(size_type i = 0; i < paths.size(); i++)
    {
      size_type thread = omp_get_thread_num();
      if(paths[i].sorted())
      {
        temp_nodes[thread].push_back(PathNode(paths[i], labels, temp_labels[thread]));
        if(temp_nodes[thread].size() >= PathGraphBuilder::WRITE_BUFFER_SIZE)
        {
          builder.write(temp_nodes[thread], temp_labels[thread], file);
        }
      }
      else
      {
        size_type first = from_index.find(paths[i].to);
        for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++)
        {
          temp_nodes[thread].push_back(PathNode(paths[i], paths[j], labels, temp_labels[thread]));
          if(temp_nodes[thread].size() >= PathGraphBuilder::WRITE_BUFFER_SIZE)
          {
            builder.write(temp_nodes[thread], temp_labels[thread], file);
          }
        }
      }
    }
    for(size_type thread = 0; thread < threads; thread++)
    {
      builder.write(temp_nodes[thread], temp_labels[thread], file);
    }
    sdsl::util::clear(paths); sdsl::util::clear(labels);

    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "PathGraph::extend(): File " << file << ": Created " << builder.graph.path_counts[file]
                << " order-" << builder.graph.k() << " paths" << std::endl;
    }

    // Sort the next generation.
    builder.sort(file);
  }
  builder.close();
  this->clear(); this->swap(builder.graph);

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::extend(): " << old_path_count << " -> " << this->size() << " paths ("
              << this->ranks() << " ranks)" << std::endl;
    std::cerr << "PathGraph::extend(): " << inGigabytes(this->bytes()) << " GB in "
              << this->files() << " file(s)" << std::endl;
  }
}

void
PathGraph::debugExtend()
{
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

    // Count the paths in the next generation.
    size_type sorted_paths = 0, unsorted_paths = 0;
    for(size_type i = 0; i < paths.size(); i++)
    {
      if(paths[i].sorted())
      {
        sorted_paths++;
      }
      else
      {
        size_type start = from_index.find(paths[i].to);
        size_type limit = start;
        while(limit < paths.size() && paths[limit].from == paths[i].to) { limit++; }
        unsorted_paths += limit - start;
        if(limit - start >= 100 && Verbosity::level >= Verbosity::FULL)
        {
          std::cerr << "PathGraph::debugExtend(): File " << file << ", path " << i
                    << " has " << (limit - start) << " extensions" << std::endl;
          std::cerr << "PathGraph::debugExtend(): The path is: ";
          paths[i].print(std::cerr, labels); std::cerr << std::endl;
        }
      }
    }

    if(Verbosity::level >= Verbosity::EXTENDED)
    {
      std::cerr << "PathGraph::debugExtend(): File " << file << ": " << this->size() << " -> "
                << (sorted_paths + unsorted_paths) << " paths ("
                << sorted_paths << " sorted, " << unsorted_paths << " unsorted)" << std::endl;
    }
  }
}

//------------------------------------------------------------------------------

void
PathGraph::read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file) const
{
  paths.resize(this->path_counts[file]);
  labels.resize(this->rank_counts[file]);

  std::ifstream path_file, rank_file;
  this->open(path_file, rank_file, file);
  if(!DiskIO::read(path_file, paths.data(), this->path_counts[file]))
  {
    std::cerr << "PathGraph::read(): Unexpected EOF in " << this->path_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(!DiskIO::read(rank_file, labels.data(), this->rank_counts[file]))
  {
    std::cerr << "PathGraph::read(): Unexpected EOF in " << this->rank_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  path_file.close(); rank_file.close();

  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "PathGraph::read(): File " << file << ": Read " << paths.size()
              << " order-" << this->k()<< " paths" << std::endl;
  }
}

//------------------------------------------------------------------------------

struct SameFromSet
{
  const PathGraphMerger& merger;
  std::vector<node_type> nodes, buffer;

  SameFromSet(const PathGraphMerger& source) :
    merger(source)
  {
  }

  inline void fromNodes(range_type range, std::vector<node_type>& to) const
  {
    to.clear();
    node_type prev = ~(node_type)0;
    for(size_type i = range.first; i <= range.second; i++)
    {
      node_type curr = this->merger.buffer[i].node.from;
      if(curr != prev) { to.push_back(curr); prev = curr; }
    }
    if(to.size() > 1) { removeDuplicates(to, false); }
  }

  inline bool operator() (range_type range)
  {
    this->fromNodes(range, this->buffer);

    // Manual comparison guarantees using a single thread.
    if(this->buffer.size() != this->nodes.size()) { return false; }
    for(size_type i = 0; i < this->buffer.size(); i++)
    {
      if(this->buffer[i] != this->nodes[i]) { return false; }
    }
    return true;
  }

  inline void select(range_type range)
  {
    this->fromNodes(range, this->nodes);
  }
};

MergedGraph::MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper, const LCP& kmer_lcp, size_type size_limit) :
  path_name(TempFile::getName(PREFIX)), rank_name(TempFile::getName(PREFIX)),
  from_name(TempFile::getName(PREFIX)), lcp_name(TempFile::getName(PREFIX)),
  path_count(0), rank_count(0), from_count(0),
  order(source.k()),
  next(mapper.alpha.sigma + 1, 0), next_from(mapper.alpha.sigma + 1, 0)
{
  WriteBuffer<PathNode>            path_file(this->path_name);
  WriteBuffer<PathNode::rank_type> rank_file(this->rank_name);
  WriteBuffer<range_type>          from_file(this->from_name);
  WriteBuffer<uint8_t>             lcp_file(this->lcp_name);

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

  PathGraphMerger merger(source, kmer_lcp);
  SameFromSet same_from_set(merger);
  size_type curr_comp = 0;  // Used to transform next.

  size_type bytes = 0;
  for(range_type range = merger.first(); !(merger.atEnd(range)); range = merger.next())
  {
    range_type path_lcp = merger.ranges.front().left_lcp; // Write this to the LCP array.
    same_from_set.select(range);
    range = merger.extendRange(same_from_set);
    merger.mergePathNodes();
    PriorityNode& curr = merger.buffer[range.second];
    curr.node.from = same_from_set.nodes[0];

    // Write the actual data
    bytes += curr.node.bytes() + (same_from_set.nodes.size() - 1) * sizeof(range_type) + 1;
    if(bytes > size_limit)
    {
      std::cerr << "MergedGraph::MergedGraph(): Size limit exceeded, construction aborted" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    writePath(curr.node, curr.label, path_file, rank_file);
    for(size_type i = 1; i < same_from_set.nodes.size(); i++)
    {
      from_file.push_back(range_type(this->path_count, same_from_set.nodes[i]));
    }
    lcp_file.push_back(path_lcp.first * mapper.order() + path_lcp.second);

    // Update the counts and the pointers to paths starting with each comp value.
    while(curr.firstLabel(0) >= this->next[curr_comp])
    {
      this->next[curr_comp] = this->path_count;
      this->next_from[curr_comp] = this->from_count;
      curr_comp++;
    }
    this->path_count++;
    this->rank_count += curr.node.ranks();
    this->from_count += same_from_set.nodes.size() - 1;
  }
  merger.close();
  path_file.close(); rank_file.close(); from_file.close(); lcp_file.close();

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "MergedGraph::MergedGraph(): " << this->size() << " paths with "
              << this->ranks() << " ranks and " << this->extra() << " additional start nodes" << std::endl;
    std::cerr << "MergedGraph::MergedGraph(): " << inGigabytes(this->bytes()) << " GB" << std::endl;
  }
}

MergedGraph::~MergedGraph()
{
  this->clear();
}

void
MergedGraph::clear()
{
  TempFile::remove(this->path_name);
  TempFile::remove(this->rank_name);
  TempFile::remove(this->from_name);
  TempFile::remove(this->lcp_name);

  this->path_count = 0; this->rank_count = 0; this->from_count = 0;
  this->order = 0;

  for(size_type i = 0; i < this->next.size(); i++) { this->next[i] = 0; }
  for(size_type i = 0; i < this->next_from.size(); i++) { this->next_from[i] = 0; }
}

//------------------------------------------------------------------------------

} // namespace gcsa
