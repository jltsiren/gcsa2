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

#ifndef _GCSA_PATH_GRAPH_H
#define _GCSA_PATH_GRAPH_H

#include "files.h"

namespace gcsa
{

//------------------------------------------------------------------------------

/*
  A path graph is a set of files, each of them containing the paths derived from one
  of the input files. The PathNodes in each file are sorted by their labels, and the
  read() member functions will also return the PathNodes in sorted order.
*/

struct PathGraph
{
  std::vector<std::string> filenames;
  std::vector<size_type>   sizes, rank_counts;

  size_type path_count, rank_count;
  size_type order;

  size_type unique, unsorted, nondeterministic;

  const static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // .gcsa

  PathGraph(const InputGraph& source, sdsl::sd_vector<>& key_exists);
  PathGraph(size_type file_count, size_type path_order);
  ~PathGraph();

  void clear();
  void swap(PathGraph& another);
  void open(std::ifstream& input, size_type file) const;

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type k() const { return this->order; }
  inline size_type files() const { return this->filenames.size(); }

  inline size_type bytes() const
  {
    return this->size() * sizeof(PathNode) + this->ranks() * sizeof(PathNode::rank_type);
  }

  // Size limits are in bytes.
  void prune(const LCP& lcp, size_type size_limit);
  void extend(size_type size_limit);

  /*
    Setting append = true has unpredictable side effects if done outside the member
    functions of PathGraph.
  */
  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels) const;
  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels,
            size_type file, bool append = false) const;

  PathGraph(const PathGraph&) = delete;
  PathGraph& operator= (const PathGraph&) = delete;
};

//------------------------------------------------------------------------------

/*
  A merged graph is a path graph with the path nodes in one file and the additional
  from nodes in another file. The PathNodes are sorted by their labels.

  FIXME Implement
*/

struct MergedGraph
{
  std::string path_name, from_name;

  size_type path_count, rank_count, from_count;
  size_type order;

  const static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // .gcsa

  explicit MergedGraph(const PathGraph& source);
  ~MergedGraph();

  void clear();

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type extras() const { return this->from_count; }
  inline size_type k() const { return this->order; }

  inline size_type bytes() const
  {
    return this->size() * sizeof(PathNode)
         + this->ranks() * sizeof(PathNode::rank_type)
         + this->extras() * sizeof(range_type);
  }

  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels,
    std::vector<range_type>& from_nodes) const;

  MergedGraph(const MergedGraph&) = delete;
  MergedGraph& operator= (const MergedGraph&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_UTILS_H
