/*
  Copyright (c) 2019 Jouni Sirén
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

#include <gcsa/internal.h>

namespace gcsa
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type DiskIO::block_size;

template<class Element> constexpr size_type ReadBuffer<Element>::READ_BUFFER_SIZE;
template<class Element> constexpr size_type ReadBuffer<Element>::MINIMUM_SIZE;

//------------------------------------------------------------------------------

// Other class variables.

std::atomic<size_type> DiskIO::read_volume(0);
std::atomic<size_type> DiskIO::write_volume(0);

//------------------------------------------------------------------------------

CounterArray::CounterArray() :
  width(8), large_value(sdsl::bits::lo_set[width]),
  total(0)
{
}

CounterArray::CounterArray(size_type n, size_type data_width) :
  data(n, 0, data_width),
  width(data_width), large_value(sdsl::bits::lo_set[width]),
  total(0)
{
}

void
CounterArray::clear()
{
  sdsl::util::clear(this->data);
  sdsl::util::clear(this->large_values);
  this->total = 0;
}

void
CounterArray::swap(CounterArray& another)
{
  if(this != &another)
  {
    this->data.swap(another.data);
    this->large_values.swap(another.large_values);
    std::swap(this->width, another.width);
    std::swap(this->large_value, another.large_value);
    std::swap(this->total, another.total);
  }
}

//------------------------------------------------------------------------------

} // namespace gcsa
