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

const std::string BINARY_EXTENSION = ".graph";
const std::string TEXT_EXTENSION = ".gcsa2";

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
    in.read((char*)(kmers.data() + old_size), header.kmer_count * sizeof(KMer));
    section++;
  }

  return kmer_length;
}

//------------------------------------------------------------------------------

size_type
readKMers(size_type files, char** base_names, std::vector<KMer>& kmers, bool binary)
{
  size_type kmer_length = ~(size_type)0;
  sdsl::util::clear(kmers);

  for(size_type i = 0; i < files; i++)
  {
    std::string filename = std::string(base_names[i]) + (binary ? BINARY_EXTENSION : TEXT_EXTENSION);
    std::ifstream input(filename.c_str(), std::ios_base::binary);
    if(!input)
    {
      std::cerr << "readKMers(): Cannot open graph file " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }

    size_type new_length = (binary ? readBinary(input, kmers, true) : readText(input, kmers, true));
    if(kmer_length == ~(size_type)0) { kmer_length = new_length; }
    else if(new_length != kmer_length)
    {
      std::cerr << "readKMers(): Invalid kmer length in graph file " << filename << std::endl;
      std::cerr << "readKMers(): Expected " << kmer_length << ", got " << new_length << std::endl;
      std::exit(EXIT_FAILURE);
    }
    input.close();
  }

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "readKMers(): Read " << kmers.size() << " kmers of length " << kmer_length << std::endl;
#endif

  /*
    If the kmer includes one or more endmarkers, the successor position is past
    the GCSA sink node. Those kmers are marked as sorted, as they cannot be
    extended.
  */
  size_type sink_node = ~(size_type)0;
  for(size_type i = 0; i < kmers.size(); i++)
  {
    if(Key::label(kmers[i].key) == 0) { sink_node = Node::id(kmers[i].from); break; }
  }
  for(size_type i = 0; i < kmers.size(); i++)
  {
    if(Node::id(kmers[i].to) == sink_node && Node::offset(kmers[i].to) > 0) { kmers[i].makeSorted(); }
  }

  return kmer_length;
}

//------------------------------------------------------------------------------

void
writeBinary(std::ostream& out, std::vector<KMer>& kmers, size_type kmer_length)
{
  GraphFileHeader header(kmers.size(), kmer_length);
  header.serialize(out);
  out.write((char*)(kmers.data()), header.kmer_count * sizeof(KMer));
}

void
writeKMers(const std::string& base_name, std::vector<KMer>& kmers, size_type kmer_length)
{
  std::string filename = base_name + BINARY_EXTENSION;
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "writeKMers(): Cannot open output file " << filename << std::endl;
    return;
  }

  writeBinary(output, kmers, kmer_length);
  output.close();

#ifdef VERBOSE_STATUS_INFO
  std::cerr << "writeKMers(): Wrote " << kmers.size() << " kmers of length " << kmer_length << std::endl;
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
  in.read((char*)this, sizeof(*this));
}

GraphFileHeader::~GraphFileHeader()
{
}

size_type
GraphFileHeader::serialize(std::ostream& out)
{
  out.write((char*)this, sizeof(*this));
  return sizeof(*this);
}

//------------------------------------------------------------------------------

} // namespace gcsa
