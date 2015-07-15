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

#ifndef _GCSA_UTILS_H
#define _GCSA_UTILS_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

#include <sdsl/wavelet_trees.hpp>

#include <omp.h>

namespace gcsa
{

//------------------------------------------------------------------------------

typedef std::uint64_t size_type;
typedef std::uint8_t  char_type;
typedef std::uint8_t  comp_type;
typedef std::uint8_t  byte_type;

const size_type WORD_BITS = 64;

const size_type KILOBYTE     = 1024;
const size_type MILLION      = 1000000;
const size_type MEGABYTE     = KILOBYTE * KILOBYTE;
const size_type GIGABYTE     = KILOBYTE * MEGABYTE;

const double KILOBYTE_DOUBLE = 1024.0;
const double MILLION_DOUBLE  = 1000000.0;
const double MEGABYTE_DOUBLE = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
const double GIGABYTE_DOUBLE = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

//------------------------------------------------------------------------------

/*
  range_type stores a closed range [first, second]. Empty ranges are indicated by
  first > second. The emptiness check uses +1 to handle a common special case
  [0, -1].
*/

typedef std::pair<size_type, size_type> range_type;

struct Range
{
  inline static size_type length(range_type range)
  {
    return range.second + 1 - range.first;
  }

  inline static bool empty(range_type range)
  {
    return (range.first + 1 > range.second + 1);
  }
};

template<class A, class B>
std::ostream& operator<<(std::ostream& stream, const std::pair<A, B>& data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

//------------------------------------------------------------------------------

template<class IntegerType>
inline size_type
bit_length(IntegerType val)
{
  return sdsl::bits::hi(val) + 1;
}

//------------------------------------------------------------------------------

const size_type FNV_OFFSET_BASIS = 0xcbf29ce484222325UL;
const size_type FNV_PRIME        = 0x100000001b3UL;

inline size_type fnv1a_hash(byte_type b, size_type seed)
{
  return (seed ^ b) * FNV_PRIME;
}

inline size_type fnv1a_hash(size_type val, size_type seed)
{
  byte_type* chars = (byte_type*)&val;
  for(size_type i = 0; i < 8; i++) { seed = fnv1a_hash(chars[i], seed); }
  return seed;
}

//------------------------------------------------------------------------------

inline double
inMegabytes(size_type bytes)
{
  return bytes / MEGABYTE_DOUBLE;
}

inline double
inBPC(size_type bytes, size_type size)
{
  return (8.0 * bytes) / size;
}

inline double
inMicroseconds(double seconds)
{
  return seconds * MILLION_DOUBLE;
}

const size_type DEFAULT_INDENT = 18;

void printHeader(const std::string& header, size_type indent = DEFAULT_INDENT);
void printSize(const std::string& header, size_type bytes, size_type data_size, size_type indent = DEFAULT_INDENT);
void printTime(const std::string& header, size_type queries, double seconds, size_type indent = DEFAULT_INDENT);

//------------------------------------------------------------------------------

double readTimer();
size_type memoryUsage(); // Peak memory usage in bytes.

// Returns the total length of the rows, excluding line ends.
size_type readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

std::string tempFile(const std::string& name_part);

template<class Iterator, class Comparator>
void
parallelQuickSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::balanced_quicksort_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelQuickSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::balanced_quicksort_tag());
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
parallelMergeSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelMergeSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
sequentialSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
sequentialSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last);
#endif
}

template<class Element>
void
removeDuplicates(std::vector<Element>& vec, bool parallel)
{
  if(parallel) { parallelQuickSort(vec.begin(), vec.end()); }
  else         { sequentialSort(vec.begin(), vec.end()); }
  vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
}

//------------------------------------------------------------------------------

/*
  GCSA uses a contiguous byte alphabet [0, sigma - 1] internally. Array C is based on the
  number of occurrences of each character in the BWT.
*/

template<class AlphabetType>
inline bool
hasChar(const AlphabetType& alpha, comp_type comp)
{
  return (alpha.C[comp + 1] > alpha.C[comp]);
}

template<class AlphabetType>
inline range_type
charRange(const AlphabetType& alpha, comp_type comp)
{
  return range_type(alpha.C[comp], alpha.C[comp + 1] - 1);
}

template<class AlphabetType>
comp_type
findChar(const AlphabetType& alpha, size_type bwt_pos)
{
  comp_type comp = 0;
  while(alpha.C[comp + 1] <= bwt_pos) { comp++; }
  return comp;
}

template<class BWTType, class AlphabetType>
inline size_type
LF(const BWTType& bwt, const AlphabetType& alpha, size_type i, comp_type comp)
{
  return alpha.C[comp] + bwt.rank(i, comp);
}

template<class BWTType, class AlphabetType>
inline range_type
LF(const BWTType& bwt, const AlphabetType& alpha, range_type range, comp_type comp)
{
  return range_type(LF(bwt, alpha, range.first, comp), LF(bwt, alpha, range.second + 1, comp) - 1);
}

//------------------------------------------------------------------------------

/*
  Some SDSL extensions.
*/

template<class Element>
size_type
write_vector(const std::vector<Element>& vec, std::ostream& out, sdsl::structure_tree_node* v, std::string name)
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(vec));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(vec.size(), out, child, "size");
  out.write((char*)(vec.data()), vec.size() * sizeof(Element));
  written_bytes += vec.size() * sizeof(Element);
  sdsl::structure_tree::add_size(v, written_bytes);
  return written_bytes;
}

template<class Element>
void
read_vector(std::vector<Element>& vec, std::istream& in)
{
  sdsl::util::clear(vec);
  size_type size = 0;
  sdsl::read_member(size, in);
  std::vector<Element> temp(size);
  in.read((char*)(temp.data()), temp.size() * sizeof(Element));
  vec.swap(temp);
}

/*
  Extracts the given range from source, overwriting target.
*/
template<class VectorType>
void
extractBits(const VectorType& source, range_type range, sdsl::bit_vector& target)
{
  if(Range::empty(range) || range.second >= source.size()) { return; }

  target = sdsl::bit_vector(Range::length(range), 0);
  for(size_type i = 0; i < target.size(); i += WORD_BITS)
  {
    size_type len = std::min(WORD_BITS, target.size() - i);
    target.set_int(i, source.get_int(range.first + i, len), len);
  }
}

/*
  Generic in-memory construction from int_vector_buffer<8> and size. Not very space-efficient, as it
  duplicates the data.
*/
template<class Type>
void
directConstruct(Type& structure, const sdsl::int_vector<8>& data)
{
  std::string ramfile = sdsl::ram_file_name(sdsl::util::to_string(&structure));
  sdsl::store_to_file(data, ramfile);
  {
    sdsl::int_vector_buffer<8> buffer(ramfile); // Must remove the buffer before removing the ramfile.
    Type temp(buffer, data.size());
    structure.swap(temp);
  }
  sdsl::ram_fs::remove(ramfile);
}

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_UTILS_H
