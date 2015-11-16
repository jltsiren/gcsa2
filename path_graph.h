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

#include "dbg.h"
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
  A merged graph is a path graph with the path nodes in one file, the labels in another
  file, and the additional from nodes in a third file. The PathNodes are sorted by their
  labels and their label pointers are set correctly.

  After writing everything to disk, the constructor opens the files and memory maps them.
*/

struct MergedGraph
{
  std::string path_name, rank_name, from_name;

  // File descriptors.
  int path_d, rank_d, from_d;

  // Memory mapped files.
  PathNode*            paths;
  PathNode::rank_type* labels;
  range_type*          from_nodes;

  size_type path_count, rank_count, from_count;
  size_type order;

  std::vector<size_type> next;  // paths[next[comp]] is the first path starting with comp.
  std::vector<size_type> next_from; // Where to find the corresponding additional from nodes.

  const static int NO_FILE = -1;
  const static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // .gcsa

  /*
    FIXME: Constructor should open the files and memory map them.
    clear() should unmap the files and close them.
  */

  MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper);
  ~MergedGraph();

  void clear();

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type extra() const { return this->from_count; }
  inline size_type k() const { return this->order; }

  inline size_type path_bytes() const { return this->size() * sizeof(PathNode); }
  inline size_type rank_bytes() const { return this->ranks() * sizeof(PathNode::rank_type); }
  inline size_type from_bytes() const { return this->extra() * sizeof(range_type); }

  inline size_type bytes() const
  {
    return this->path_bytes() + this->rank_bytes() + this->from_bytes();
  }

  void read(std::vector<PathNode>& _paths, std::vector<PathNode::rank_type>& _labels,
    std::vector<range_type>& _from_nodes) const;

  void* map(const std::string& filename, const std::string& file_type, int& fd, size_type n);

  MergedGraph(const MergedGraph&) = delete;
  MergedGraph& operator= (const MergedGraph&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_UTILS_H
