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

#ifndef _GCSA_INTERNAL_H
#define _GCSA_INTERNAL_H

#include <map>

// C++ threads for DiskIO, ReadBuffer.
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "utils.h"

namespace gcsa
{

/*
  internal.h: Internal support structures.
*/

//------------------------------------------------------------------------------

/*
  Utility methods for disk I/O and read/write volume measurement. These methods don't use
  mutexes / critical sections for performance reasons.
*/

struct DiskIO
{
  static std::atomic<size_type> read_volume, write_volume;

  template<class Element>
  inline static void read(std::istream& in, Element* data, size_type n = 1)
  {
    read_volume += n * sizeof(Element);
    in.read((char*)data, n * sizeof(Element));
  }

  template<class Element>
  inline static void write(std::ostream& out, const Element* data, size_type n = 1)
  {
    write_volume += n * sizeof(Element);
    out.write((const char*)data, n * sizeof(Element));
  }
};

//------------------------------------------------------------------------------

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
  A simple counter array that uses an int_vector<0> for most counters and stores large values
  in an std::map.
*/
struct CounterArray
{
  sdsl::int_vector<0>            data;
  std::map<size_type, size_type> large_values;
  size_type                      width, large_value;
  size_type                      total;

  CounterArray();
  CounterArray(size_type n, size_type data_width);

  inline size_type size() const { return this->data.size(); }
  inline size_type sum() const { return this->total; }

  inline size_type operator[] (size_type i) const
  {
    return (this->data[i] == this->large_value ? this->large_values.at(i) : this->data[i]);
  }

  inline void increment(size_type i)
  {
    if(this->data[i] == this->large_value) { this->large_values[i]++; }
    else
    {
      this->data[i]++;
      if(this->data[i] == this->large_value) { this->large_values[i] = this->large_value; }
    }
    this->total++;
  }

  inline void increment(size_type i, size_type val)
  {
    if(this->data[i] == this->large_value) { this->large_values[i] += val; }
    else if(this->data[i] + val >= this->large_value)
    {
      this->large_values[i] = this->data[i] + val;
      this->data[i] = this->large_value;
    }
    else { this->data[i] += val; }
    this->total += val;
  }

  void clear();
  void swap(CounterArray& another);
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
  void append(std::deque<Element>& buffer, size_type n, std::vector<Element>& read_buffer);
  void read(std::vector<Element>& read_buffer);
  void seek(size_type i);
  void finish() { this->offset = this->size(); }
  bool finished() { return (this->offset >= this->size()); }

  std::ifstream file;
  size_type elements, offset;

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
  this->elements = fileSize(this->file) / sizeof(Element);
  this->offset = 0;
}

template<class Element>
void
ElementReader<Element>::append(std::deque<Element>& buffer, size_type n, std::vector<Element>& read_buffer)
{
  n = std::min(this->size() - this->offset, n);
  size_type block_size = std::min(n, (size_type)(read_buffer.capacity()));

  read_buffer.resize(block_size);
  this->read(read_buffer);
  buffer.insert(buffer.end(), read_buffer.begin(), read_buffer.end());
  read_buffer.clear();
}

template<class Element>
void
ElementReader<Element>::read(std::vector<Element>& read_buffer)
{
  if(read_buffer.size() > this->size() - this->offset) { read_buffer.resize(this->size() - this->offset); }
  DiskIO::read(this->file, read_buffer.data(), read_buffer.size());
  this->offset += read_buffer.size();
}

template<class Element>
void
ElementReader<Element>::seek(size_type i)
{
  i = std::min(i, this->size());
  this->file.seekg(i * sizeof(Element), std::ios_base::beg);
  this->offset = i;
}

//------------------------------------------------------------------------------

/*
  A buffer for reading a file of Elements sequentially. The buffer contains Elements
  offset to offset + buffer.size() - 1. The offset is moved by calling seek() or when
  accessing an element before the current offset.

  A separate thread is spawned for reading in the background.

  The Reader must implement the following interface:

  Reader()                        Default constructor.
  init(parameters...)             Initialize the reader.
  close()                         Close the reader.
  size()                          Return the size of the source data.
  append(buffer, n, read_buffer)  Append n elements to buffer. read_buffer can be used as
                                  a buffer if necessary. Its capacity should be > 0.
  read(read_buffer)               Read read_buffer.size() elements.
  seek(i)                         Go to offset i.
  finish()                        Finish reading.
  finished()                      Has reading finished (at end or by calling finish())?
*/

template<class Element, class Reader = ElementReader<Element>>
struct ReadBuffer
{
  // Main thread.
  size_type               offset;
  std::deque<Element>     buffer;
  Reader                  reader;

  // Reader thread.
  std::vector<Element>    read_buffer;
  std::mutex              mtx;
  std::condition_variable empty;
  std::thread             reader_thread;

  // Minimum size when refilling the buffer.
  const static size_type BUFFER_SIZE = MEGABYTE;

  // Refill the buffer if its size falls below this threshold.
  const static size_type MINIMUM_SIZE = BUFFER_SIZE / 2;

  // Read this many elements at once.
  const static size_type READ_BUFFER_SIZE = BUFFER_SIZE;

  ReadBuffer();
  ~ReadBuffer();

  template<class... ReaderParameters> void init(ReaderParameters... parameters);
  void close();

  inline size_type size() const { return this->reader.size(); }
  inline bool buffered(size_type i) const
  {
    return (i >= this->offset && i < this->offset + this->buffer.size());
  }

  /*
    Beware: Calling seek() may invalidate the reference. The same may also happen when
    calling operator[] with a non-buffered position.
  */
  inline Element& operator[] (size_type i)
  {
    if(!(this->buffered(i))) { this->read(i); }
    return this->buffer[i - this->offset];
  }
  void seek(size_type i); // Set offset to i.

  // Internal functions.
  bool fill();            // Fill the read buffer.
  void read(size_type i); // Read i into buffer, possibly seeking backwards.
  void addBlock();        // Add a block into buffer, assuming that this read holds mtx.

  ReadBuffer(const ReadBuffer&) = delete;
  ReadBuffer& operator= (const ReadBuffer&) = delete;
};

template<class Element, class Reader>
ReadBuffer<Element, Reader>::ReadBuffer()
{
  this->offset = 0;
  this->read_buffer.reserve(READ_BUFFER_SIZE);
}

template<class Element, class Reader>
ReadBuffer<Element, Reader>::~ReadBuffer()
{
  this->close();
}

template<class Element, class Reader>
void
readerThread(ReadBuffer<Element, Reader>* buffer)
{
  while(!(buffer->fill()));
}

template<class Element, class Reader>
template<class... ReaderParameters>
void
ReadBuffer<Element, Reader>::init(ReaderParameters... parameters)
{
  this->reader.init(parameters...);
  this->reader_thread = std::thread(readerThread<Element, Reader>, this);
}

template<class Element, class Reader>
void
ReadBuffer<Element, Reader>::close()
{
  // We need to stop the reader thread.
  this->mtx.lock();
  this->reader.finish();
  this->read_buffer.clear();
  this->empty.notify_one();
  this->mtx.unlock();
  if(this->reader_thread.joinable()) { this->reader_thread.join(); }

  sdsl::util::clear(this->buffer);
  sdsl::util::clear(this->read_buffer);
  this->reader.close();
}

template<class Element, class Reader>
void
ReadBuffer<Element, Reader>::seek(size_type i)
{
  if(i >= this->size()) { return; }

  if(this->buffered(i))
  {
    while(this->offset < i) { this->buffer.pop_front(); this->offset++; }
  }
  else
  {
    std::unique_lock<std::mutex> lock(this->mtx);
    this->buffer.clear(); this->read_buffer.clear();
    this->reader.seek(i); this->offset = i;
  }
  if(this->buffer.size() < MINIMUM_SIZE)
  {
    std::unique_lock<std::mutex> lock(this->mtx);
    this->addBlock();
    this->empty.notify_one();
  }
}

template<class Element, class Reader>
bool
ReadBuffer<Element, Reader>::fill()
{
  std::unique_lock<std::mutex> lock(this->mtx);
  this->empty.wait(lock, [this]() { return read_buffer.empty(); } );

  this->read_buffer.resize(READ_BUFFER_SIZE);
  this->reader.read(this->read_buffer);

  return this->reader.finished();
}

template<class Element, class Reader>
void
ReadBuffer<Element, Reader>::read(size_type i)
{
  if(i >= this->size()) { return; }
  if(i < this->offset) { this->seek(i); return; }

  std::unique_lock<std::mutex> lock(this->mtx);
  while(!(this->buffered(i))) { this->addBlock(); }
  this->empty.notify_one();
}

template<class Element, class Reader>
void
ReadBuffer<Element, Reader>::addBlock()
{
  if(this->read_buffer.empty())
  {
    this->reader.append(this->buffer, READ_BUFFER_SIZE, this->read_buffer);
  }
  else
  {
    this->buffer.insert(this->buffer.end(), this->read_buffer.begin(), this->read_buffer.end());
    this->read_buffer.clear();
  }
}

//------------------------------------------------------------------------------

/*
  A simple wrapper for buffered writing of elementary types.
*/

template<class Element>
struct WriteBuffer
{
  explicit WriteBuffer(const std::string& filename, size_type _buffer_size = MEGABYTE);
  ~WriteBuffer();
  void close();

  inline void push_back(Element value)
  {
    this->buffer.push_back(value);
    if(buffer.size() >= this->buffer_size)
    {
      DiskIO::write(this->file, this->buffer.data(), this->buffer.size());
      this->buffer.clear();
    }
  }

  std::ofstream        file;
  std::vector<Element> buffer;
  size_type            buffer_size;

  WriteBuffer(const WriteBuffer&) = delete;
  WriteBuffer& operator= (const WriteBuffer&) = delete;
};

template<class Element>
WriteBuffer<Element>::WriteBuffer(const std::string& filename, size_type _buffer_size) :
  file(filename.c_str(), std::ios_base::binary), buffer(), buffer_size(_buffer_size)
{
  if(!(this->file))
  {
    std::cerr << "WriteBuffer::WriteBuffer(): Cannot open output file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->buffer.reserve(this->buffer_size);
}

template<class Element>
WriteBuffer<Element>::~WriteBuffer()
{
  this->close();
}

template<class Element>
void
WriteBuffer<Element>::close()
{
  if(this->buffer.size() > 0)
  {
    DiskIO::write(this->file, this->buffer.data(), this->buffer.size());
  }
  this->file.close();
  sdsl::util::clear(this->buffer);
}

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_INTERNAL_H
