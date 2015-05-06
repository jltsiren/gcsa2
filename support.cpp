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

#include <sstream>
#include <vector>

#include "support.h"

namespace gcsa
{

//------------------------------------------------------------------------------

/*
  The default alphabet interprets \0 and $ as endmarkers, ACGT and acgt as ACGT,
  # as a the label of the source node, and the and the remaining characters as N.
*/

const sdsl::int_vector<8> Alphabet::DEFAULT_CHAR2COMP =
{
  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 6,  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

const sdsl::int_vector<8> Alphabet::DEFAULT_COMP2CHAR = { '$', 'A', 'C', 'G', 'T', 'N', '#' };

//------------------------------------------------------------------------------

Alphabet::Alphabet() :
  char2comp(DEFAULT_CHAR2COMP), comp2char(DEFAULT_COMP2CHAR),
  C(sdsl::int_vector<64>(DEFAULT_COMP2CHAR.size() + 1, 0)),
  sigma(DEFAULT_COMP2CHAR.size())
{
}

Alphabet::Alphabet(const Alphabet& a)
{
  this->copy(a);
}

Alphabet::Alphabet(Alphabet&& a)
{
  *this = std::move(a);
}

Alphabet::Alphabet(const sdsl::int_vector<64>& counts,
  const sdsl::int_vector<8>& _char2comp, const sdsl::int_vector<8>& _comp2char) :
  char2comp(_char2comp), comp2char(_comp2char),
  C(sdsl::int_vector<64>(_comp2char.size() + 1, 0)),
  sigma(_comp2char.size())
{
  for(size_type i = 0; i < counts.size(); i++) { this->C[i + 1] = this->C[i] + counts[i]; }
}

Alphabet::~Alphabet()
{
}

void
Alphabet::copy(const Alphabet& a)
{
  this->char2comp = a.char2comp;
  this->comp2char = a.comp2char;
  this->C = a.C;
  this->sigma = a.sigma;
}

void
Alphabet::swap(Alphabet& a)
{
  if(this != &a)
  {
    this->char2comp.swap(a.char2comp);
    this->comp2char.swap(a.comp2char);
    this->C.swap(a.C);
    std::swap(this->sigma, a.sigma);
  }
}

Alphabet&
Alphabet::operator=(const Alphabet& a)
{
  if(this != &a) { this->copy(a); }
  return *this;
}

Alphabet&
Alphabet::operator=(Alphabet&& a)
{
  if(this != &a)
  {
    this->char2comp = std::move(a.char2comp);
    this->comp2char = std::move(a.comp2char);
    this->C = std::move(a.C);
    this->sigma = a.sigma;
  }
  return *this;
}

Alphabet::size_type
Alphabet::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->char2comp.serialize(out, child, "char2comp");
  written_bytes += this->comp2char.serialize(out, child, "comp2char");
  written_bytes += this->C.serialize(out, child, "C");
  written_bytes += sdsl::write_member(this->sigma, out, child, "sigma");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
Alphabet::load(std::istream& in)
{
  this->char2comp.load(in);
  this->comp2char.load(in);
  this->C.load(in);
  sdsl::read_member(this->sigma, in);
}

//------------------------------------------------------------------------------

std::string
Key::decode(const Alphabet& alpha, key_type key, size_type kmer_length)
{
  key = kmer(key);
  kmer_length = std::min(kmer_length, (size_type)16);

  std::string res(kmer_length, '\0');
  for(size_type i = 1; i <= kmer_length; i++)
  {
    res[kmer_length - i] = alpha.comp2char[key & 0x7];
    key >>= 3;
  }

  return res;
}

//------------------------------------------------------------------------------

node_type
Node::encode(const std::string& token)
{
  size_t separator = 0;
  size_type node = std::stoul(token, &separator);
  if(separator >= token.length())
  {
    std::cerr << "Node::encode(): Invalid position token " << token << std::endl;
    return 0;
  }

  std::string temp = token.substr(separator + 1);
  size_type offset = std::stoul(temp);
  if(offset > OFFSET_MASK)
  {
    std::cerr << "Node::encode(): Offset " << offset << " too large!" << std::endl;
    return 0;
  }

  return encode(node, offset);
}

std::string
Node::decode(node_type node)
{
  std::ostringstream ss;
  ss << id(node) << ':' << offset(node);
  return ss.str();
}

//------------------------------------------------------------------------------

KMer::KMer()
{
}

KMer::KMer(const std::vector<std::string>& tokens, const Alphabet& alpha, size_type successor)
{
  byte_type predecessors = chars(tokens[2], alpha);
  byte_type successors = chars(tokens[3], alpha);
  this->key = Key::encode(alpha, tokens[0], predecessors, successors);
  this->from = Node::encode(tokens[1]);
  this->to = Node::encode(tokens[successor]);
}

bool
KMer::tokenize(const std::string& line, std::vector<std::string>& tokens)
{
  {
    std::string token;
    std::istringstream ss(line);
    while(std::getline(ss, token, '\t'))
    {
      tokens.push_back(token);
    }
    if(tokens.size() < 4 || tokens.size() > 5)
    {
      std::cerr << "KMer::tokenize(): The kmer line must contain 4 or 5 tokens." << std::endl;
      std::cerr << "KMer::tokenize(): The line was: " << line << std::endl;
      return false;
    }
  }

  // Split the list of successor positions into separate tokens.
  if(tokens.size() == 5)
  {
    std::string destinations = tokens[4], token;
    std::istringstream ss(destinations);
    tokens.resize(4);
    while(std::getline(ss, token, ','))
    {
      tokens.push_back(token);
    }
  }
  else  // Use the source node as the destination.
  {
    tokens.push_back(tokens[1]);
  }

  return true;
}

byte_type
KMer::chars(const std::string& token, const Alphabet& alpha)
{
  byte_type val = 0;
  for(size_type i = 0; i < token.length(); i += 2) { val |= 1 << alpha.char2comp[token[i]]; }
  return val;
}

std::ostream&
operator<< (std::ostream& out, const KMer& kmer)
{
  out << "(key " << Key::kmer(kmer.key)
      << ", in " << (size_type)(Key::predecessors(kmer.key))
      << ", out " << (size_type)(Key::successors(kmer.key))
      << ", from " << Node::decode(kmer.from)
      << ", to " << Node::decode(kmer.to) << ")";
  return out;
}

//------------------------------------------------------------------------------

} // namespace gcsa
