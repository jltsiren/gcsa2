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

#include <unistd.h>

#include <gcsa/gcsa.h>

using namespace gcsa;

//------------------------------------------------------------------------------

void identifyGCSA(const std::string& input_name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: gcsa_format input" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string input_name = argv[1];
  identifyGCSA(input_name);

  return 0;
}

//------------------------------------------------------------------------------

/*
  This is GCSA file format header version 0. It was used in pre-releases 0.1 to 0.4.
*/

struct GCSAHeader_0
{
  uint64_t path_nodes;
  uint64_t order;

  const static uint32_t VERSION = 0;

  GCSAHeader_0();

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check() const;
};

GCSAHeader_0::GCSAHeader_0() :
  path_nodes(0), order(0)
{
}

size_type
GCSAHeader_0::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(this->path_nodes, out, child, "path_nodes");
  written_bytes += sdsl::write_member(this->order, out, child, "order");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSAHeader_0::load(std::istream& in)
{
  sdsl::read_member(this->path_nodes, in);
  sdsl::read_member(this->order, in);
}

bool
GCSAHeader_0::check() const
{
  return true;
}

std::ostream& operator<<(std::ostream& stream, const GCSAHeader_0& header)
{
  return stream << "GCSA header version " << GCSAHeader_0::VERSION << ": "
                << header.path_nodes << " path nodes, order " << header.order;
}

//------------------------------------------------------------------------------

bool
tryVersion0(std::ifstream& input)
{
  std::streampos pos = input.tellg();
  GCSAHeader_0 header;
  header.load(input);
  input.seekg(pos);

  if(header.check())
  {
    std::cout << header << std::endl;
    return true;
  }
  return false;
}

bool
tryCurrentVersion(std::ifstream& input)
{
  std::streampos pos = input.tellg();
  GCSAHeader header;
  header.load(input);
  input.seekg(pos);

  for(uint32_t version = GCSAHeader::VERSION; version >= GCSAHeader::MIN_VERSION; version--)
  {
    if(header.check(version))
    {
      std::cout << header << std::endl;
      return true;
    }
  }
  if(header.checkNew())
  {
    std::cout << header << std::endl;
    std::cout << "Note: The file is newer than this version of GCSA" << std::endl;
    return true;
  }

  return false;
}

void
identifyGCSA(const std::string& input_name)
{
  std::cout << "GCSA file format inspector" << std::endl;
  std::cout << std::endl;

  std::cout << "Input: " << input_name << std::endl;
  std::cout << std::endl;

  std::ifstream input(input_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "identifyGCSA(): Cannot open input file " << input_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(!tryCurrentVersion(input))
  {
    std::cout << "File format cannot be identified, trying version 0" << std::endl;
    tryVersion0(input);
  }
  input.close();

  std::cout << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
