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
Key::decode(key_type key, size_type kmer_length, const Alphabet& alpha)
{
  key = label(key);
  kmer_length = std::min(kmer_length, MAX_LENGTH);

  std::string res(kmer_length, '\0');
  for(size_type i = 1; i <= kmer_length; i++)
  {
    res[kmer_length - i] = alpha.comp2char[key & CHAR_MASK];
    key >>= CHAR_WIDTH;
  }

  return res;
}

//------------------------------------------------------------------------------

node_type
Node::encode(const std::string& token)
{
  size_t separator = 0;
  size_type node = std::stoul(token, &separator);
  if(separator + 1 >= token.length())
  {
    std::cerr << "Node::encode(): Invalid position token " << token << std::endl;
    return 0;
  }

  bool reverse_complement = false;
  if(token[separator + 1] == '-')
  {
    reverse_complement = true;
    separator++;
  }

  std::string temp = token.substr(separator + 1);
  size_type offset = std::stoul(temp);
  if(offset > OFFSET_MASK)
  {
    std::cerr << "Node::encode(): Offset " << offset << " too large" << std::endl;
    return 0;
  }

  return encode(node, offset, reverse_complement);
}

std::string
Node::decode(node_type node)
{
  std::ostringstream ss;
  ss << id(node) << ':';
  if(rc(node)) { ss << '-'; }
  ss << offset(node);
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
  out << "(key " << Key::label(kmer.key)
      << ", in " << (size_type)(Key::predecessors(kmer.key))
      << ", out " << (size_type)(Key::successors(kmer.key))
      << ", from " << Node::decode(kmer.from)
      << ", to " << Node::decode(kmer.to) << ")";
  return out;
}

void
uniqueKeys(std::vector<KMer>& kmers, std::vector<key_type>& keys, sdsl::int_vector<0>& last_char, bool print)
{
  if(kmers.empty()) { return; }
  parallelQuickSort(kmers.begin(), kmers.end());

  // Pass 1: Count the number of unique keys.
  size_type total_keys = 1;
  for(size_type i = 1; i < kmers.size(); i++)
  {
    if(Key::label(kmers[i].key) != Key::label(kmers[i - 1].key)) { total_keys++; }
  }
  if(print)
  {
    std::cout << "Unique keys: " << total_keys << std::endl;
  }

  // Pass 2: Create the merged key array and the last character array for edge generation.
  // Replace the kmer values with ranks in the key array.
  keys = std::vector<key_type>(total_keys, 0);
  last_char = sdsl::int_vector<0>(total_keys, 0, Key::CHAR_WIDTH);
  keys[0] = kmers[0].key; last_char[0] = Key::last(kmers[0].key);
  kmers[0].key = Key::replace(kmers[0].key, 0);
  for(size_type kmer = 1, key = 0; kmer < kmers.size(); kmer++)
  {
    if(Key::label(kmers[kmer].key) == Key::label(keys[key]))
    {
      keys[key] = Key::merge(keys[key], kmers[kmer].key);
    }
    else
    {
      key++; keys[key] = kmers[kmer].key; last_char[key] = Key::last(kmers[kmer].key);
    }
    kmers[kmer].key = Key::replace(kmers[kmer].key, key);
  }
}

//------------------------------------------------------------------------------

std::vector<PathNode::rank_type>
PathNode::dummyRankVector()
{
  std::vector<rank_type> temp;
  temp.reserve(LABEL_LENGTH + 1);
  return temp;
}

PathNode::PathNode(const KMer& kmer, std::vector<PathNode::rank_type>& labels)
{
  this->from = kmer.from; this->to = kmer.to;
  this->fields = 0;

  if(kmer.sorted()) { this->makeSorted(); }
  this->setPredecessors(Key::predecessors(kmer.key));
  this->setOrder(1); this->setLCP(1);

  this->setPointer(labels.size());
  labels.push_back(Key::label(kmer.key));
  labels.push_back(0);  // Dummy value; the last label is not in use.
}

PathNode::PathNode(const PathNode& source,
    const std::vector<PathNode::rank_type>& old_labels, std::vector<PathNode::rank_type>& new_labels)
{
  this->from = source.from; this->to = source.to;
  this->fields = source.fields;

  this->setPointer(new_labels.size());
  for(size_type i = 0, j = source.pointer(); i < source.ranks(); i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
}

PathNode::PathNode(const PathNode& left, const PathNode& right,
    const std::vector<PathNode::rank_type>& old_labels, std::vector<PathNode::rank_type>& new_labels)
{
  this->from = left.from; this->to = right.to;
  if(right.sorted()) { this->makeSorted(); }

  this->fields = 0;
  this->setPredecessors(left.predecessors());

  size_type left_order = left.order();
  size_type new_order = left_order + right.order();
  this->setOrder(new_order); this->setLCP(left_order + right.lcp());

  this->setPointer(new_labels.size());
  for(size_type i = 0, j = left.pointer(); i < left_order; i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
  for(size_type i = 0, j = right.pointer(); i < right.ranks(); i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
}

PathNode::PathNode(std::ifstream& in, std::vector<PathNode::rank_type>& labels)
{
  in.read((char*)this, sizeof(*this));
  this->setPointer(labels.size());

  rank_type buffer[LABEL_LENGTH + 1];
  in.read((char*)buffer, this->ranks() * sizeof(rank_type));
  labels.insert(labels.end(), buffer, buffer + this->ranks());
}

size_type
PathNode::serialize(std::ostream& out, const std::vector<rank_type>& labels) const
{
  size_type bytes = sizeof(*this);
  out.write((char*)this, sizeof(*this));

  size_type rank_bytes = this->ranks() * sizeof(rank_type);
  out.write((char*)(labels.data() + this->pointer()), rank_bytes);
  bytes += rank_bytes;

  return bytes;
}

PathNode::PathNode()
{
  this->from = 0; this->to = 0;
  this->fields = 0;
}

PathNode::PathNode(std::vector<rank_type>& labels)
{
  this->from = 0; this->to = 0;
  this->fields = 0;
  this->setPointer(labels.size());
}

PathNode::PathNode(const PathNode& source)
{
  this->copy(source);
}

PathNode::PathNode(PathNode&& source)
{
  *this = std::move(source);
}

PathNode::~PathNode()
{
}

PathNode&
PathNode::operator= (const PathNode& source)
{
  if(&source != this)
  {
    this->copy(source);
  }
  return *this;
}

void
PathNode::copy(const PathNode& source)
{
  this->from = source.from; this->to = source.to;
  this->fields = source.fields;
}

PathNode&
PathNode::operator= (PathNode&& source)
{
  if(&source != this)
  {
    this->from = std::move(source.from);
    this->to = std::move(source.to);
    this->fields = std::move(source.fields);
  }
  return *this;
}

void
PathNode::print(std::ostream& out, const std::vector<PathNode::rank_type>& labels) const
{
  out << "(" << Node::decode(this->from) << " -> " << Node::decode(this->to);
  out << "; o" << this->order();
  for(size_type i = 0; i < this->order(); i++)
  {
    out << (i == 0 ? "; [" : ", ") << this->firstLabel(i, labels);
  }
  for(size_type i = 0; i < this->order(); i++)
  {
    out << (i == 0 ? " to " : ", ") << this->lastLabel(i, labels);
  }
  out << "])";
}

//------------------------------------------------------------------------------

bool
PathNode::intersect(const PathLabel& first, const PathLabel& last,
  const std::vector<PathNode::rank_type>& labels) const
{
  PathLabel my_first = this->firstLabel(labels);

  if(my_first.compare(first))
  {
    PathLabel my_last = this->lastLabel(labels);
    return !(my_last.compare(first));
  }
  else
  {
    return !(last.compare(my_first));
  }
}

size_type
PathNode::min_lcp(const PathNode& another, const std::vector<PathNode::rank_type>& labels) const
{
  size_type ord = std::min(this->order(), another.order());
  for(size_type i = 0; i < ord; i++)
  {
    if(this->firstLabel(i, labels) != another.lastLabel(i, labels)) { return i; }
  }
  return ord;
}

size_type
PathNode::max_lcp(const PathNode& another, const std::vector<PathNode::rank_type>& labels) const
{
  size_type ord = std::min(this->order(), another.order());
  for(size_type i = 0; i < ord; i++)
  {
    if(this->lastLabel(i, labels) != another.firstLabel(i, labels)) { return i; }
  }
  return ord;
}

//------------------------------------------------------------------------------

LCP::LCP()
{
}

LCP::LCP(const std::vector<key_type>& keys, size_type _kmer_length)
{
  this->kmer_length = _kmer_length;
  this->total_keys = keys.size();
  {
    sdsl::int_vector<0> temp(keys.size(), 0, bit_length(this->kmer_length - 1));
    for(size_type i = 1; i < keys.size(); i++)
    {
      temp[i] = Key::lcp(keys[i - 1], keys[i], this->kmer_length);
    }
    this->kmer_lcp.swap(temp);
  }
  sdsl::util::assign(this->lcp_rmq, rmq_type(&(this->kmer_lcp)));
}

range_type
LCP::min_lcp(const PathNode& a, const PathNode& b, const std::vector<LCP::rank_type>& labels) const
{
  size_type order = std::min(a.order(), b.order());
  range_type lcp(a.min_lcp(b, labels), 0);
  if(lcp.first < order)
  {
    size_type right = std::min((size_type)(b.lastLabel(lcp.first, labels)), this->total_keys - 1);
    lcp.second =
      this->kmer_lcp[this->lcp_rmq(a.firstLabel(lcp.first, labels) + 1, right)];
  }
  return lcp;
}

range_type
LCP::max_lcp(const PathNode& a, const PathNode& b, const std::vector<LCP::rank_type>& labels) const
{
  size_type order = std::min(a.order(), b.order());
  range_type lcp(a.max_lcp(b, labels), 0);
  if(lcp.first < order)
  {
    lcp.second =
      this->kmer_lcp[this->lcp_rmq(a.lastLabel(lcp.first, labels) + 1, b.firstLabel(lcp.first, labels))];
  }
  return lcp;
}

void
LCP::swap(LCP& another)
{
  std::swap(this->kmer_length, another.kmer_length);
  std::swap(this->total_keys, another.total_keys);
  this->kmer_lcp.swap(another.kmer_lcp);
  this->lcp_rmq.swap(another.lcp_rmq);
}

//------------------------------------------------------------------------------

} // namespace gcsa
