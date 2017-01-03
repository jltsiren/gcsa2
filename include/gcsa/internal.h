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
#include <condition_variable>
#include <mutex>
#include <thread>

#include <gcsa/utils.h>

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

  inline void clear() { this->data.clear(); }
  inline void resize(size_type n) { this->data.resize(n); }

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
  A buffer that keeps Elements [offset, offset + size - 1] in memory. It supports
  seeking to a new position and adding Elements to the end.
*/

template<class Element>
struct BufferWindow
{
  std::deque<Element> data;
  size_type           offset;

  BufferWindow();
  ~BufferWindow();

  inline size_type size() const { return this->data.size(); }
  inline bool buffered(size_type i) const
  {
    return (i >= this->offset && i < this->offset + this->size());
  }

  /*
    Beware: Calling seek() may invalidate the reference. The same may also happen when
    calling operator[] with a non-buffered position.
  */
  inline Element& operator[] (size_type i) { return this->data[i - this->offset]; }
  inline const Element& operator[] (size_type i) const { return this->data[i - this->offset]; }

  void clear();
  void seek(size_type i);
  void swap(BufferWindow<Element>& another);

  inline void push_back(const Element& element) { this->data.push_back(element); }

  template<class Iterator>
  void insert(Iterator from, Iterator to) { this->data.insert(this->data.end(), from, to); }
};

template<class Element>
BufferWindow<Element>::BufferWindow() :
  offset(0)
{
}

template<class Element>
BufferWindow<Element>::~BufferWindow()
{
}

template<class Element>
void
BufferWindow<Element>::clear()
{
  sdsl::util::clear(this->data);
  this->offset = 0;
}

template<class Element>
void
BufferWindow<Element>::seek(size_type i)
{
  if(this->buffered(i))
  {
    for(size_type j = this->offset; j < i; j++) { this->data.pop_front(); }
  }
  else { this->data.clear(); }
  this->offset = i;
}

template<class Element>
void
BufferWindow<Element>::swap(BufferWindow<Element>& another)
{
  if(this != &another)
  {
    this->data.swap(another.data);
    std::swap(this->offset, another.offset);
  }
}

//------------------------------------------------------------------------------

/*
  A buffer for reading a file of Elements sequentially. Function seek() moves the start
  of the buffer to the new position. Accesses before the current start call seek(), while
  accesses after it expand the buffer until the requested position is contained in it.

  A separate thread is spawned for reading in the background. The reader thread stops
  when it reaches the end of the file.
*/

template<class Element>
struct ReadBuffer
{
  // Main thread.
  BufferWindow<Element>   buffer;

  // File
  std::ifstream           file;
  size_type               elements, file_offset;

  // Reader thread.
  std::vector<Element>    read_buffer;
  std::mutex              mtx;
  std::condition_variable empty;  // Is read_buffer empty?
  std::thread             reader_thread;

  // Read this many elements at once.
  const static size_type READ_BUFFER_SIZE = MEGABYTE;

  // Refill the buffer if its size falls below this threshold.
  const static size_type MINIMUM_SIZE = READ_BUFFER_SIZE / 2;

  ReadBuffer();
  ~ReadBuffer();

  void open(const std::string& filename);
  void close();

  inline size_type size() const { return this->elements; }

  /*
    Beware: Calling seek() may invalidate the reference. The same may also happen when
    calling operator[] with a non-buffered position.
  */
  inline Element& operator[] (size_type i)
  {
    if(!(this->buffer.buffered(i))) { this->read(i); }
    return this->buffer[i];
  }

  void seek(size_type i); // Set offset to i.

  // Internal functions.
  bool fill();            // Fill the read buffer.
  void read(size_type i); // Read i into buffer, possibly seeking backwards.
  void forceRead();       // Add elements into buffer, assuming that the current thread holds the mutex.

  ReadBuffer(const ReadBuffer&) = delete;
  ReadBuffer& operator= (const ReadBuffer&) = delete;
};

template<class Element>
ReadBuffer<Element>::ReadBuffer()
{
  this->elements = 0; this->file_offset = 0;
  size_type buffer_size = READ_BUFFER_SIZE;  // avoid direct use
  this->read_buffer.reserve(buffer_size);
}

template<class Element>
ReadBuffer<Element>::~ReadBuffer()
{
  this->close();
}

template<class Element>
void
readerThread(ReadBuffer<Element>* buffer)
{
  while(!(buffer->fill()));
}

template<class Element>
void
ReadBuffer<Element>::open(const std::string& filename)
{
  if(this->file.is_open())
  {
    std::cerr << "ReadBuffer::open(): The file is already open" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  this->file.open(filename.c_str(), std::ios_base::binary);
  if(!(this->file))
  {
    std::cerr << "ReadBuffer::open(): Cannot open input file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->elements = fileSize(this->file) / sizeof(Element);
  this->file_offset = 0;

  this->reader_thread = std::thread(readerThread<Element>, this);
}

template<class Element>
void
ReadBuffer<Element>::close()
{
  // We need to stop the reader thread.
  this->mtx.lock();
  this->file_offset = this->size();
  this->read_buffer.clear();
  this->empty.notify_one();
  this->mtx.unlock();
  if(this->reader_thread.joinable()) { this->reader_thread.join(); }

  this->file.close();
  this->elements = 0;

  sdsl::util::clear(this->buffer);
}

template<class Element>
void
ReadBuffer<Element>::seek(size_type i)
{
  if(i >= this->size()) { return; }

  // Move the buffer to the new position.
  // Clear the buffer and seek in the file if file offset is wrong.
  this->buffer.seek(i);
  if(!(this->buffer.buffered(i)))
  {
    std::unique_lock<std::mutex> lock(this->mtx);
    if(this->file_offset != i + this->read_buffer.size())
    {
      this->read_buffer.clear();
      this->file.seekg(i * sizeof(Element), std::ios_base::beg);
      this->file_offset = i;
    }
  }

  // Force read but only if there is still something to read.
  if(this->buffer.size() < MINIMUM_SIZE && i + this->buffer.size() < this->size())
  {
    std::unique_lock<std::mutex> lock(this->mtx);
    this->forceRead();
    this->empty.notify_one();
  }
}

template<class Element>
bool
ReadBuffer<Element>::fill()
{
  std::unique_lock<std::mutex> lock(this->mtx);
  this->empty.wait(lock, [this]() { return read_buffer.empty(); } );

  size_type read_buffer_size = READ_BUFFER_SIZE;  // avoid direct use
  this->read_buffer.resize(std::min(read_buffer_size, this->size() - this->file_offset));
  DiskIO::read(this->file, this->read_buffer.data(), this->read_buffer.size());
  this->file_offset += this->read_buffer.size();

  return (this->file_offset >= this->size());
}

template<class Element>
void
ReadBuffer<Element>::read(size_type i)
{
  if(i >= this->size()) { return; }
  if(i < this->buffer.offset) { this->seek(i); return; }

  std::unique_lock<std::mutex> lock(this->mtx);
  while(!(this->buffer.buffered(i))) { this->forceRead(); }
  this->empty.notify_one();
}

template<class Element>
void
ReadBuffer<Element>::forceRead()
{
  if(this->read_buffer.empty())
  {
    size_type read_buffer_size = READ_BUFFER_SIZE;  // avoid direct use
    this->read_buffer.resize(std::min(read_buffer_size, this->size() - this->file_offset));
    DiskIO::read(this->file, this->read_buffer.data(), this->read_buffer.size());
    this->file_offset += this->read_buffer.size();
  }

  this->buffer.insert(this->read_buffer.begin(), this->read_buffer.end());
  this->read_buffer.clear();
}

//------------------------------------------------------------------------------

/*
  A simple wrapper for buffered writing of elementary types.
*/

template<class Element>
struct WriteBuffer
{
  WriteBuffer();
  explicit WriteBuffer(const std::string& filename, size_type _buffer_size = MEGABYTE);
  ~WriteBuffer();

  void open(const std::string& filename, size_type _buffer_size = MEGABYTE);
  void close();

  inline size_type size() const { return this->elements; }

  inline void push_back(Element value)
  {
    this->buffer.push_back(value); this->elements++;
    if(buffer.size() >= this->buffer_size)
    {
      DiskIO::write(this->file, this->buffer.data(), this->buffer.size());
      this->buffer.clear();
    }
  }

  std::ofstream        file;
  std::vector<Element> buffer;
  size_type            buffer_size, elements;

  WriteBuffer(const WriteBuffer&) = delete;
  WriteBuffer& operator= (const WriteBuffer&) = delete;
};

template<class Element>
WriteBuffer<Element>::WriteBuffer() :
  buffer_size(0), elements(0)
{
}

template<class Element>
WriteBuffer<Element>::WriteBuffer(const std::string& filename, size_type _buffer_size)
{
  this->open(filename, _buffer_size);
}

template<class Element>
WriteBuffer<Element>::~WriteBuffer()
{
  this->close();
}

template<class Element>
void
WriteBuffer<Element>::open(const std::string& filename, size_type _buffer_size)
{
  this->file.open(filename.c_str(), std::ios_base::binary);
  if(!(this->file))
  {
    std::cerr << "WriteBuffer::open(): Cannot open output file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }

  this->buffer_size = _buffer_size; this->elements = 0;
  this->buffer.reserve(this->buffer_size);
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
  this->buffer_size = 0;
  this->elements = 0;
}

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_INTERNAL_H
