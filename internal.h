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

#ifndef _GCSA_INTERNAL_H
#define _GCSA_INTERNAL_H

#include <map>

#include "utils.h"

namespace gcsa
{

/*
  internal.h: Internal support structures.
*/

//------------------------------------------------------------------------------

template<class ValueType, class Getter>
struct ValueIndex
{
  sdsl::sd_vector<>               values;     // Marks the values that are present.
  sdsl::sd_vector<>::rank_1_type  value_rank;

  sdsl::bit_vector                first_occ;  // Marks the first occurrence of each rank.
  sdsl::bit_vector::select_1_type first_select;

  ValueIndex(const std::vector<ValueType>& input)
  {
    std::vector<size_type> buffer;
    this->first_occ = sdsl::bit_vector(input.size(), 0);

    size_type prev = ~(size_type)0;
    for(size_type i = 0; i < input.size(); i++)
    {
      size_type curr = Getter::get(input[i]);
      if(curr != prev)
      {
        buffer.push_back(curr);
        this->first_occ[i] = 1;
        prev = curr;
      }
    }

    // Fills in values, but only works if there are any values to fill
    if(buffer.size() > 0)
    {
      sdsl::sd_vector<> temp(buffer.begin(), buffer.end());
      this->values.swap(temp);
      sdsl::util::clear(buffer);
    }

    sdsl::util::init_support(this->value_rank, &(this->values));
    sdsl::util::init_support(this->first_select, &(this->first_occ));
  }

  // Finds the first occurrence of the value.
  size_type find(size_type value) const
  {
    if(value >= this->values.size() || this->values[value] == 0) { return this->first_occ.size(); }
    return this->first_select(this->value_rank(value) + 1);
  }

  ValueIndex(const ValueIndex&) = delete;
  ValueIndex& operator= (const ValueIndex&) = delete;
};

//------------------------------------------------------------------------------

/*
  A simple byte array that stores large values in an std::map. Values start as 0s.
  Supports access and increment().
*/
struct SLArray
{
  std::vector<byte_type> data;
  std::map<size_type, size_type> large_values;

  const static byte_type LARGE_VALUE = 255;

  explicit SLArray(size_type n);

  inline bool size() const { return data.size(); }

  inline size_type operator[] (size_type i) const
  {
    return (this->data[i] == LARGE_VALUE ? this->large_values.at(i) : this->data[i]);
  }

  inline void increment(size_type i)
  {
    if(this->data[i] == LARGE_VALUE) { this->large_values[i]++; }
    else
    {
      this->data[i]++;
      if(this->data[i] == LARGE_VALUE) { this->large_values[i] = LARGE_VALUE; }
    }
  }

  void clear();
};

//------------------------------------------------------------------------------

template<class Element>
struct PriorityQueue
{
  std::vector<Element> data;

  PriorityQueue();
  explicit PriorityQueue(size_type n) : data(n) { }

  void resize(size_type n) { this->data.resize(n); }

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

  inline Element& operator[] (size_type i) { return this->data[i]; }
  inline const Element& operator[] (size_type i) const { return this->data[i]; }

  void heapify();
};

template<class Element>
void
PriorityQueue<Element>::heapify()
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

/*
  A simple wrapper for reading a file of Elements. Use with ReadBuffer.
*/

template<class Element>
struct ElementReader
{
  ElementReader() { this->elements = 0; }
  ~ElementReader() { this->close(); }

  void init(const std::string& filename);
  void close() { this->file.close(); this->elements = 0; }

  inline size_type size() const { return this->elements; }
  inline void seek(size_type i) { this->file.seekg(i * sizeof(Element), std::ios_base::beg); }
  inline void read(std::vector<Element>& buffer)
  {
    DiskIO::read(this->file, buffer.data(), buffer.size());
  }

  const static bool SEEKABLE = true;

  std::ifstream file;
  size_type elements;

  ElementReader(const ElementReader&) = delete;
  ElementReader& operator= (const ElementReader&) = delete;
};

template<class Element>
void
ElementReader<Element>::init(const std::string& filename)
{
  if(this->file.is_open())
  {
    std::cerr << "ElementReader::init(): The file is already open" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  this->file.open(filename.c_str(), std::ios_base::binary);
  if(!(this->file))
  {
    std::cerr << "ElementReader::init(): Cannot open input file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->elements = fileSize(this->file);
}

//------------------------------------------------------------------------------

/*
  A buffer for reading a file of Elements sequentially. The buffer contains Elements
  offset to offset + buffer.size() - 1. Reader must implement the same interface as
  ElementReader.
*/

template<class Element, class Reader = ElementReader<Element>>
struct ReadBuffer
{
  Reader              reader;
  size_type           offset;
  std::deque<Element> buffer;

  // After seek(), buffer size should in [MINIMUM_SIZE, BUFFER_SIZE].
  const static size_type BUFFER_SIZE = MEGABYTE;
  const static size_type MINIMUM_SIZE = BUFFER_SIZE / 2;

  ReadBuffer() { this->offset = 0; }
  ~ReadBuffer() { this->clear(); }

  template<class... ReaderParameters> void init(ReaderParameters... parameters)
  {
    this->reader.init(parameters...);
  }
  void clear() { sdsl::util::clear(this->buffer); this->reader.close(); }

  inline size_type size() const { return this->reader.size(); }
  inline bool buffered(size_type i) const
  {
    return (i >= this->offset && i < this->offset + this->buffer.size());
  }

  void seek(size_type i);     // Set offset to i.
  void prefetch(size_type i); // Prefetch [offset, i].

  /*
    Beware: Any buffer manipulations (seek(), fill(), read(), popUntil()) may invalidate
    the reference. The same may also happen when calling operator[] with a non-buffered
    position.
  */
  inline Element& operator[] (size_type i)
  {
    if(!(this->buffered(i))) { if(Reader::SEEKABLE) { this->seek(i); } else { this->prefetch(i); } }
    return this->buffer[i - this->offset];
  }

  void resize(size_type target_size = BUFFER_SIZE);
  inline void pop() { this->buffer.pop_front(); this->offset++; }

  ReadBuffer(const ReadBuffer&) = delete;
  ReadBuffer& operator= (const ReadBuffer&) = delete;
};

template<class Element, class Reader>
void
ReadBuffer<Element, Reader>::seek(size_type i)
{
  if(i >= this->size()) { return; }
  if(this->buffered(i))
  {
    while(this->offset < i) { this->pop(); }
    if(this->buffer.size() < MINIMUM_SIZE) { this->resize(); }
    return;
  }

  this->offset += this->buffer.size(); this->buffer.clear();
  if(Reader::SEEKABLE) { this->reader.seek(i); }
  else
  {
    std::vector<Element> temp(BUFFER_SIZE);
    while(this->offset < i)
    {
      temp.resize(std::min(i - this->offset, BUFFER_SIZE));
      this->reader.read(temp); this->offset += temp.size();
    }
  }
  this->offset = i; this->resize();
}

template<class Element, class Reader>
void
ReadBuffer<Element, Reader>::prefetch(size_type i)
{
  if(i < this->offset + this->buffer.size()) { return; }
  this->resize(i + 1 - this->offset);
}

template<class Element, class Reader>
void
ReadBuffer<Element, Reader>::resize(size_type target_size)
{
  target_size = std::min(target_size, this->size() - this->offset);
  if(this->buffer.size() >= target_size) { return; }

  size_type count = target_size - this->buffer.size();
  std::vector<Element> temp(count);
  this->reader.read(temp);
  this->buffer.insert(this->buffer.end(), temp.begin(), temp.end());
}

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_INTERNAL_H
