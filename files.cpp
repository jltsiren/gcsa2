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

#include "files.h"

namespace gcsa
{

//------------------------------------------------------------------------------

bool
tokenize(const std::string& line, std::vector<std::string>& tokens)
{
  {
    std::string token;
    std::istringstream ss(line);
    while(std::getline(ss, token, '\t'))
    {
      tokens.push_back(token);
    }
    if(tokens.size() != 5)
    {
      std::cerr << "tokenize(): The kmer line must contain 5 tokens" << std::endl;
      std::cerr << "tokenize(): The line was: " << line << std::endl;
      return false;
    }
  }

  // Split the list of successor positions into separate tokens.
  std::string destinations = tokens[4], token;
  std::istringstream ss(destinations);
  tokens.resize(4);
  while(std::getline(ss, token, ','))
  {
    tokens.push_back(token);
  }

  return true;
}

size_type
readText(std::istream& in, std::vector<KMer>& kmers, bool append)
{
  if(!append) { sdsl::util::clear(kmers); }

  Alphabet alpha;
  size_type kmer_length = ~(size_type)0;
  while(true)
  {
    std::string line;
    std::getline(in, line);
    if(in.eof()) { break; }

    std::vector<std::string> tokens;
    if(!tokenize(line, tokens)) { continue; }
    if(kmer_length == ~(size_type)0)
    {
      kmer_length = tokens[0].length();
      if(kmer_length == 0 || kmer_length > Key::MAX_LENGTH)
      {
        std::cerr << "readText(): Invalid kmer length: " << kmer_length << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else if(tokens[0].length() != kmer_length)
    {
      std::cerr << "readText(): Invalid kmer length: " << tokens[0].length()
                << " (expected " << kmer_length << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    for(size_type successor = 4; successor < tokens.size(); successor++)
    {
      kmers.push_back(KMer(tokens, alpha, successor));
    }
  }

  return kmer_length;
}

//------------------------------------------------------------------------------

size_type
readBinary(std::istream& in, std::vector<KMer>& kmers, bool append)
{
  if(!append) { sdsl::util::clear(kmers); }

  size_type kmer_length = ~(size_type)0;
  size_type section = 0;
  while(true)
  {
    GraphFileHeader header(in);
    if(in.eof()) { break; }
    if(header.flags != 0)
    {
      std::cerr << "readBinary(): Invalid flags in section " << section << ": " << header.flags << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(kmer_length == ~(size_type)0)
    {
      kmer_length = header.kmer_length;
      if(kmer_length == 0 || kmer_length > Key::MAX_LENGTH)
      {
        std::cerr << "readBinary(): Invalid kmer length in section " << section << ": "
                  << kmer_length << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else if(header.kmer_length != kmer_length)
    {
      std::cerr << "readBinary(): Invalid kmer length in section " << section << ": "
                << header.kmer_length << " (expected " << kmer_length << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    size_type old_size = kmers.size();
    kmers.resize(old_size + header.kmer_count);
    DiskIO::read(in, kmers.data() + old_size, header.kmer_count);
    section++;
  }

  return kmer_length;
}

//------------------------------------------------------------------------------

void
writeBinary(std::ostream& out, std::vector<KMer>& kmers, size_type kmer_length)
{
  GraphFileHeader header(kmers.size(), kmer_length);
  header.serialize(out);
  DiskIO::write(out, kmers.data(), header.kmer_count);
}

void
writeKMers(const std::string& base_name, std::vector<KMer>& kmers, size_type kmer_length)
{
  std::string filename = base_name + InputGraph::BINARY_EXTENSION;
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "writeKMers(): Cannot open output file " << filename << std::endl;
    return;
  }

  writeBinary(output, kmers, kmer_length);
  output.close();

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "writeKMers(): Wrote " << kmers.size() << " " << kmer_length << "-mers" << std::endl;
#endif
}

//------------------------------------------------------------------------------

GraphFileHeader::GraphFileHeader() :
  flags(0), kmer_count(0), kmer_length(0)
{
}

GraphFileHeader::GraphFileHeader(size_type kmers, size_type length) :
  flags(0), kmer_count(kmers), kmer_length(length)
{
}

GraphFileHeader::GraphFileHeader(std::istream& in)
{
  DiskIO::read(in, this);
}

GraphFileHeader::~GraphFileHeader()
{
}

size_type
GraphFileHeader::serialize(std::ostream& out)
{
  DiskIO::write(out, this);
  return sizeof(*this);
}

//------------------------------------------------------------------------------

range_type
parseKMer(const std::string& kmer_line)
{
  std::string extensions;
  range_type result(InputGraph::UNKNOWN, 0);
  {
    std::string token;
    std::istringstream ss(kmer_line);
    if(!std::getline(ss, token, '\t')) { return result; }
    result.first = token.length();  // kmer length
    for(size_type i = 0; i < 3; i++)
    {
      if(!std::getline(ss, token, '\t')) { return result; }
    }
    if(!std::getline(ss, extensions, '\t')) { return result; }
  }

  {
    std::string token;
    std::istringstream ss(extensions);
    while(std::getline(ss, token, ',')) { result.second++; }
  }

  return result;
}

/*
  If the kmer ends with an endmarker, it cannot be extended, and we mark it
  sorted. If its label is not unique, it will be treated as a nondeterministic
  path. We also mark nodes ending with the source marker sorted, because there
  should be only one such path, which is duplicated in each input file.
*/
void
markSourceSinkNodes(std::vector<KMer>& kmers)
{
  for(size_type i = 0; i < kmers.size(); i++)
  {
    if(Key::last(kmers[i].key) == Alphabet::SOURCE_COMP ||
       Key::last(kmers[i].key) == Alphabet::SINK_COMP)
    {
      kmers[i].makeSorted();
    }
  }
}

//------------------------------------------------------------------------------

const std::string InputGraph::BINARY_EXTENSION = ".graph";
const std::string InputGraph::TEXT_EXTENSION = ".gcsa2";

InputGraph::InputGraph(size_type file_count, char** base_names, bool binary_format)
{
  this->binary = binary_format;
  this->kmer_count = 0; this->kmer_length = UNKNOWN;
  for(size_type file = 0; file < file_count; file++)
  {
    std::string filename = std::string(base_names[file]) + (binary ? BINARY_EXTENSION : TEXT_EXTENSION);
    this->filenames.push_back(filename);
  }
  this->sizes = std::vector<size_type>(file_count, 0);

  // Read the files and determine kmer_count, kmer_length.
  for(size_type file = 0; file < this->files(); file++)
  {
    std::ifstream input; this->open(input, file);
    if(this->binary)
    {
      while(true)
      {
        GraphFileHeader header(input);
        if(input.eof()) { break; }
        this->kmer_count += header.kmer_count;
        this->setK(header.kmer_length, file);
        this->sizes[file] += header.kmer_count;
        input.seekg(header.kmer_count * sizeof(KMer), std::ios_base::cur);
      }
    }
    else
    {
      while(true)
      {
        std::string line;
        std::getline(input, line);
        if(input.eof()) { break; }
        range_type temp = parseKMer(line);
        this->setK(temp.first, file);
        this->kmer_count += temp.second; this->sizes[file] += temp.second;
      }

    }
    input.close();
  }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "InputGraph::InputGraph(): " << this->size() << " kmers in "
            << this->files() << " file(s)" << std::endl;
#endif
}

void
InputGraph::open(std::ifstream& input, size_type file) const
{
  if(file >= this->files())
  {
    std::cerr << "InputGraph::open(): Invalid file number: " << file << std::endl;
    std::exit(EXIT_FAILURE);
  }

  input.open(this->filenames[file].c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "InputGraph::open(): Cannot open graph file " << this->filenames[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void
InputGraph::setK(size_type new_k, size_type file)
{
  if(this->k() == UNKNOWN) { this->kmer_length = new_k; }
  this->checkK(new_k, file);
}

void
InputGraph::checkK(size_type new_k, size_type file) const
{
  if(new_k != this->k())
  {
    std::cerr << "InputGraph::checkK(): Invalid kmer length in graph file " << this->filenames[file] << std::endl;
    std::cerr << "InputGraph::checkK(): Expected " << this->k() << ", got " << new_k << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//------------------------------------------------------------------------------

void
InputGraph::read(std::vector<KMer>& kmers) const
{
  sdsl::util::clear(kmers);
  kmers.reserve(this->size());

  for(size_type file = 0; file < this->files(); file++)
  {
    this->read(kmers, file, true);
  }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "InputGraph::read(): Read " << kmers.size() << " " << this->k() << "-mers" << std::endl;
#endif

  markSourceSinkNodes(kmers);
}

void
InputGraph::read(std::vector<KMer>& kmers, size_type file, bool append) const
{
  if(!append) { sdsl::util::clear(kmers); }

  std::ifstream input; this->open(input, file);
  if(!append) { kmers.reserve(this->sizes[file]); }
  size_type new_k = (this->binary ? readBinary(input, kmers, append) : readText(input, kmers, append));
  this->checkK(new_k, file);
  input.close();

#ifdef VERBOSE_STATUS_INFO
  if(!append)
  {
    std::cerr << "InputGraph::read(): Read " << kmers.size() << " " << this->k() << "-mers"
              << " from " << this->filenames[file] << std::endl;
  }
#endif

  if(!append) { markSourceSinkNodes(kmers); }
}

void
InputGraph::read(std::vector<key_type>& keys) const
{
  sdsl::util::clear(keys);
  keys.reserve(this->size());

  // Read the keys.
  for(size_type file = 0; file < this->files(); file++)
  {
    std::vector<KMer> kmers;
    this->read(kmers, file, false);
    for(size_type i = 0; i < this->sizes[file]; i++)
    {
      keys.push_back(kmers[i].key);
    }
  }

  // Sort the keys and merge the ones sharing the same label.
  parallelQuickSort(keys.begin(), keys.end());
  size_type i = 0;
  for(size_type j = 1; j < keys.size(); j++)
  {
    if(Key::label(keys[i]) == Key::label(keys[j])) { keys[i] = Key::merge(keys[i], keys[j]); }
    else { i++; keys[i] = keys[j]; }
  }
  keys.resize(i + 1);

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "InputGraph::read(): " << keys.size() << " unique keys" << std::endl;
#endif
}

//------------------------------------------------------------------------------

} // namespace gcsa
