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
readText(size_type files, char** filenames, std::vector<KMer>& kmers)
{
  Alphabet alpha;
  size_type kmer_length = ~(size_type)0;
  sdsl::util::clear(kmers);

  for(size_type i = 0; i < files; i++)
  {
    std::string filename = std::string(filenames[i]) + TEXT_EXTENSION;
    std::ifstream input(filename.c_str(), std::ios_base::binary);
    if(!input)
    {
      std::cerr << "readText(): Cannot open graph file " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }

    while(input)
    {
      std::string line;
      std::getline(input, line);
      if(line.length() == 0) { continue; }

      std::vector<std::string> tokens;
      if(!tokenize(line, tokens)) { continue; }
      if(kmer_length == ~(size_type)0)
      {
        kmer_length = tokens[0].length();
        if(kmer_length == 0 || kmer_length > Key::MAX_LENGTH)
        {
          std::cerr << "readText(): Invalid kmer length in graph file " << filename
                    << ": " << kmer_length << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      else if(tokens[0].length() != kmer_length)
      {
        std::cerr << "readText(): Invalid kmer length in graph file " << filename
                  << ": expected " << kmer_length << ", got " << tokens[0].length() << std::endl;
        std::exit(EXIT_FAILURE);
      }

      for(size_type successor = 4; successor < tokens.size(); successor++)
      {
        kmers.push_back(KMer(tokens, alpha, successor));
      }
    }
    input.close();
  }

  return kmer_length;
}

//------------------------------------------------------------------------------

size_type
readBinary(size_type files, char** filenames, std::vector<KMer>& kmers)
{
  sdsl::util::clear(kmers);
  size_type kmer_length = ~(size_type)0;

  for(size_type i = 0; i < files; i++)
  {
    std::string filename = std::string(filenames[i]) + BINARY_EXTENSION;
    std::ifstream input(filename.c_str(), std::ios_base::binary);
    if(!input)
    {
      std::cerr << "readBinary(): Cannot open graph file " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }

    size_type section = 0;
    while(true)
    {
      GraphFileHeader header(input);
      if(input.eof()) { break; }
      if(header.flags != 0)
      {
        std::cerr << "readBinary(): Invalid flags in graph file " << filename
                  << ", section " << section << std::endl;
        std::cerr << "readBinary(): Expected 0, got " << header.flags << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if(kmer_length == ~(size_type)0)
      {
        kmer_length = header.kmer_length;
        if(kmer_length == 0 || kmer_length > Key::MAX_LENGTH)
        {
          std::cerr << "readBinary(): Invalid kmer length in graph file " << filename
                    << ", section " << section << ": " << kmer_length << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      else if(header.kmer_length != kmer_length)
      {
        std::cerr << "readBinary(): Invalid kmer length in graph file " << filename
                  << ", section " << section << std::endl;
        std::cerr << "readBinary(): Expected " << kmer_length << ", got " << header.kmer_length << std::endl;
        std::exit(EXIT_FAILURE);
      }

      size_type old_size = kmers.size();
      kmers.resize(old_size + header.kmer_count);
      input.read((char*)(kmers.data() + old_size), header.kmer_count * sizeof(KMer));
      section++;
    }
    input.close();
  }

  return kmer_length;
}

//------------------------------------------------------------------------------

size_type
readKMers(size_type files, char** filenames, std::vector<KMer>& kmers, bool binary, bool print)
{
  size_type kmer_length = (binary ? readBinary(files, filenames, kmers) : readText(files, filenames, kmers));
  if(print)
  {
    std::cout << "Read " << kmers.size() << " kmers of length " << kmer_length << std::endl;
  }

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

void
writeKMers(std::vector<KMer>& kmers, size_type kmer_length, const std::string& base_name, bool print)
{
  std::string filename = base_name + BINARY_EXTENSION;
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "writeKMers(): Cannot open output file " << filename << std::endl;
    return;
  }

  GraphFileHeader header(kmers.size(), kmer_length);
  header.serialize(output);
  output.write((char*)(kmers.data()), header.kmer_count * sizeof(KMer));
  output.close();

  if(print)
  {
    std::cout << "Wrote " << header.kmer_count << " kmers of length " << header.kmer_length << std::endl;
  }
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
