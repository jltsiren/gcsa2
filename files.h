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

extern const std::string BINARY_EXTENSION;  // .graph
extern const std::string TEXT_EXTENSION;    // .gcsa2

size_type readKMers(size_type files, char** filenames, std::vector<KMer>& kmers,
  bool binary, bool print = false);

struct GraphFileHeader
{
  size_type flags;
  size_type kmer_count;
  size_type kmer_length;

  GraphFileHeader();
  explicit GraphFileHeader(std::istream& in);
  ~GraphFileHeader();

  size_type serialize(std::ostream& out);
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_UTILS_H
