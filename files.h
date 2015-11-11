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

#ifndef _GCSA_FILES_H
#define _GCSA_FILES_H

#include "support.h"

namespace gcsa
{

//------------------------------------------------------------------------------

struct GraphFileHeader
{
  size_type flags;
  size_type kmer_count;
  size_type kmer_length;

  GraphFileHeader();
  GraphFileHeader(size_type kmers, size_type length);
  explicit GraphFileHeader(std::istream& in);
  ~GraphFileHeader();

  size_type serialize(std::ostream& out);
};

/*
  These functions read the input until eof. They do not close the input stream. The
  return value is kmer length.
*/
size_type readBinary(std::istream& in, std::vector<KMer>& kmers, bool append = false);
size_type readText(std::istream& in, std::vector<KMer>& kmers, bool append = false);

// FIXME Later: writeText()
void writeBinary(std::ostream& out, std::vector<KMer>& kmers, size_type kmer_length);
void writeKMers(const std::string& base_name, std::vector<KMer>& kmers, size_type kmer_length);

//------------------------------------------------------------------------------

/*
  An input graph is just a set of input files.
*/

struct InputGraph
{
  std::vector<std::string> filenames;
  std::vector<size_type>   sizes;

  bool binary;
  size_type kmer_count, kmer_length;

  const static size_type UNKNOWN = ~(size_type)0;
  const static std::string BINARY_EXTENSION;  // .graph
  const static std::string TEXT_EXTENSION;    // .gcsa2

  InputGraph(size_type file_count, char** base_names, bool binary_format);

  void open(std::ifstream& input, size_type file) const;
  void setK(size_type new_k, size_type file);
  void checkK(size_type new_k, size_type file) const;

  inline size_type size() const { return this->kmer_count; }
  inline size_type k() const { return this->kmer_length; }
  inline size_type files() const { return this->filenames.size(); }

  /*
    Setting append = true has unpredictable side effects if done outside the member
    functions of InputGraph.
  */
  void read(std::vector<KMer>& kmers) const;
  void read(std::vector<KMer>& kmers, size_type file, bool append = false) const;
  void read(std::vector<key_type>& keys) const;
};

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

  const static std::string PREFIX;  // .gcsa

  PathGraph(const InputGraph& source, sdsl::sd_vector<>& key_exists);
  ~PathGraph();

  void open(std::ifstream& input, size_type file) const;
  void clear();

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type k() const { return this->order; }
  inline size_type files() const { return this->filenames.size(); }

  inline size_type bytes() const
  {
    return this->size() * sizeof(PathNode) + this->ranks() * sizeof(PathNode::rank_type);
  }

  /*
    Setting append = true has unpredictable side effects if done outside the member
    functions of PathGraph.
  */
  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels) const;
  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels,
            size_type file, bool append = false) const;
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_UTILS_H
