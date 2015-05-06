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

#ifndef _GCSA_SUPPORT_H
#define _GCSA_SUPPORT_H

#include "utils.h"

namespace gcsa
{

//------------------------------------------------------------------------------

template<class ByteVector>
void characterCounts(const ByteVector& sequence, const sdsl::int_vector<8>& char2comp, sdsl::int_vector<64>& counts);

/*
  This replaces the SDSL byte_alphabet. The main improvements are:
    - The alphabet can be built from an existing sequence.
    - The comp order does not need to be the same as character order, as long as \0 is the first character.
*/

class Alphabet
{
public:
  typedef gcsa::size_type size_type;
  const static size_type MAX_SIGMA = 256;

  const static sdsl::int_vector<8> DEFAULT_CHAR2COMP;
  const static sdsl::int_vector<8> DEFAULT_COMP2CHAR;

  Alphabet();
  Alphabet(const Alphabet& a);
  Alphabet(Alphabet&& a);
  ~Alphabet();

  /*
    ByteVector only has to support operator[] and size(). If there is a clearly faster way for
    sequential access, function characterCounts() should be specialized.
  */
  template<class ByteVector>
  explicit Alphabet(const ByteVector& sequence,
    const sdsl::int_vector<8>& _char2comp = DEFAULT_CHAR2COMP,
    const sdsl::int_vector<8>& _comp2char = DEFAULT_COMP2CHAR) :
    char2comp(_char2comp), comp2char(_comp2char),
    C(sdsl::int_vector<64>(_comp2char.size() + 1, 0)),
    sigma(_comp2char.size())
  {
    if(sequence.size() == 0) { return; }

    characterCounts(sequence, this->char2comp, this->C);
    for(size_type i = 0, sum = 0; i < this->C.size(); i++)
    {
      size_type temp = this->C[i]; this->C[i] = sum; sum += temp;
    }
  }

  /*
    The counts array holds character counts for all comp values.
  */
  explicit Alphabet(const sdsl::int_vector<64>& counts,
    const sdsl::int_vector<8>& _char2comp = DEFAULT_CHAR2COMP,
    const sdsl::int_vector<8>& _comp2char = DEFAULT_COMP2CHAR);

  void swap(Alphabet& a);
  Alphabet& operator=(const Alphabet& a);
  Alphabet& operator=(Alphabet&& a);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  sdsl::int_vector<8>  char2comp, comp2char;
  sdsl::int_vector<64> C;
  size_type            sigma;

private:
  void copy(const Alphabet& a);
};  // class Alphabet

template<class ByteVector>
void
characterCounts(const ByteVector& sequence, const sdsl::int_vector<8>& char2comp, sdsl::int_vector<64>& counts)
{
  for(size_type c = 0; c < counts.size(); c++) { counts[c] = 0; }
  for(size_type i = 0; i < sequence.size(); i++) { counts[char2comp[sequence[i]]]++; }
}

//------------------------------------------------------------------------------

/*
  This interface is intended for indexing kmers of length 16 or less on an alphabet of size
  8 or less. The kmer is encoded as an 64-bit integer (most significant bit first):
    - 16x3 bits for the kmer, with high-order characters 0s when necessary
    - 8 bits for marking which predecessors are present
    - 8 bits for marking which successors are present
*/

typedef std::uint64_t key_type;

struct Key
{
  const static size_type CHAR_WIDTH = 3;
  const static size_type CHAR_MASK = 0x7;
  const static size_type MAX_LENGTH = 16;

  inline static key_type encode(const Alphabet& alpha, const std::string& label,
    byte_type pred, byte_type succ)
  {
    key_type value = 0;
    for(size_type i = 0; i < label.length(); i++)
    {
      value = (value << CHAR_WIDTH) | alpha.char2comp[label[i]];
    }
    value = (value << 8) | pred;
    value = (value << 8) | succ;
    return value;
  }

  static std::string decode(const Alphabet& alpha, key_type key, size_type kmer_length);

  inline static size_type kmer(key_type key) { return (key >> 16); }
  inline static byte_type predecessors(key_type key) { return (key >> 8) & 0xFF; }
  inline static byte_type successors(key_type key) { return key & 0xFF; }
  inline static char_type last(key_type key) { return (key >> 16) & CHAR_MASK; }

  inline static key_type merge(key_type key1, key_type key2) { return (key1 | (key2 & 0xFFFF)); }
  inline static key_type replace(key_type key, size_type kmer_val)
  {
    return (kmer_val << 16) | (key & 0xFFFF);
  }
};

//------------------------------------------------------------------------------

typedef std::uint64_t node_type;

struct Node
{
  const static size_type OFFSET_BITS = 6;
  const static size_type OFFSET_MASK = 63;

  inline static node_type encode(size_type node_id, size_type node_offset)
  {
    return (node_id << OFFSET_BITS) | node_offset;
  }

  static node_type encode(const std::string& token);
  static std::string decode(node_type node);

  inline static size_type id(node_type node) { return node >> OFFSET_BITS; }
  inline static size_type offset(node_type node) { return node & OFFSET_MASK; }
};

//------------------------------------------------------------------------------

struct KMer
{
  key_type  key;
  node_type from, to;

  KMer();
  KMer(const std::vector<std::string>& tokens, const Alphabet& alpha, size_type successor);

  inline bool
  operator< (const KMer& another) const
  {
    return (Key::kmer(this->key) < Key::kmer(another.key));
  }

  inline bool sorted() const { return (this->from == this->to); }
  inline void makeSorted() { this->to = this->from; }

  static bool tokenize(const std::string& line, std::vector<std::string>& tokens);
  static byte_type chars(const std::string& token, const Alphabet& alpha);
};

std::ostream& operator<< (std::ostream& out, const KMer& kmer);

inline bool
operator< (key_type key, const KMer& kmer)
{
  return (Key::kmer(key) < Key::kmer(kmer.key));
}

/*
  This function does several things:

  1. Sorts the kmer array by the kmer labels encoded in the key
  2. Builds an array of unique kmer labels, with the predecessor and successor
     fields merged from the original kmers.
  3. Stores the last character of each unique kmer label in an array.
  4. Replaces the kmer labels in the keys by their ranks.
*/
void uniqueKeys(std::vector<KMer>& kmers, std::vector<key_type>& keys, sdsl::int_vector<0>& last_chars);

//------------------------------------------------------------------------------

/*
  The node type used during doubling. As in the original GCSA, from and to are nodes
  in the original graph, denoting a path as a semiopen range [from, to). If
  from == to, the path will not be extended, because it already has a unique label.
  rank_type is the integer type used to store ranks of the original kmers.
  During edge generation, to will be used to store the number of outgoing edges.

  FIXME Remove the template, move code to .cpp?
*/

template<class rank_type = std::uint32_t, size_type label_length = 8>
struct DoublingNode
{
  node_type from, to;
  rank_type label[label_length];

  inline bool operator< (const DoublingNode& another) const
  {
    if(this->from != another.from) { return (this->from < another.from); }
    return (this->to < another.to);
  }

  inline bool sorted() const { return (this->from == this->to); }
  inline void makeSorted() { this->to = this->from; }

  inline size_type outdegree() const { return this->to; }

  DoublingNode(const KMer& kmer)
  {
    this->from = kmer.from; this->to = kmer.to;
    this->label[0] = Key::kmer(kmer.key);
    for(size_type i = 1; i < label_length; i++) { label[i] = 0; }
  }

  DoublingNode(const DoublingNode& left, const DoublingNode& right, size_type old_order)
  {
    this->from = left.from; this->to = right.to;
    for(size_type i = 0; i < old_order; i++) { this->label[i] = left.label[i]; }
    for(size_type i = old_order; i < 2 * old_order; i++) { this->label[i] = right.label[i - old_order]; }
    for(size_type i = 2 * old_order; i < label_length; i++) { this->label[i] = 0; }
  }

  explicit DoublingNode(std::ifstream& in)
  {
    sdsl::read_member(&(this->from), in); sdsl::read_member(&(this->to), in);
    in.read((char*)(this->label), label_length * sizeof(rank_type));
  }

  size_type serialize(std::ostream& out) const
  {
    size_type bytes = 0;
    bytes += sdsl::write_member(this->from, out); bytes += sdsl::write_member(this->to, out);
    out.write((char*)(this->label), label_length * sizeof(rank_type));
    bytes += label_length * sizeof(rank_type);
    return bytes;
  }

  DoublingNode() {}
  DoublingNode(const DoublingNode& source) { this->copy(source); }
  DoublingNode(DoublingNode&& source) { *this = std::move(source); }
  ~DoublingNode() {}

  void copy(const DoublingNode& source)
  {
    if(&source != this)
    {
      this->from = source.from; this->to = source.to;
      for(size_type i = 0; i < label_length; i++) { this->label[i] = source.label[i]; }
    }
  }

  void swap(DoublingNode& another)
  {
    if(&another != this)
    {
      std::swap(this->from, another.from); std::swap(this->to, another.to);
      for(size_type i = 0; i < label_length; i++) { std::swap(this->label[i], another.label[i]); }
    }
  }

  DoublingNode& operator= (const DoublingNode& source)
  {
    this->copy(source);
    return *this;
  }

  DoublingNode& operator= (DoublingNode&& source)
  {
    if(&source != this)
    {
      this->from = std::move(source.from);
      this->to = std::move(source.to);
      for(size_type i = 0; i < label_length; i++) { this->label[i] = std::move(source.label[i]); }
    }
    return *this;
  }
};

template<class rank_type = std::uint32_t, size_type label_length = 8>
struct NodeLabelComparator
{
  typedef DoublingNode<rank_type, label_length> dn_type;

  size_type max_length;

  NodeLabelComparator(size_type len = label_length) : max_length(len) {}

  inline bool operator() (const dn_type& a, const dn_type& b) const
  {
    for(size_t i = 0; i < this->max_length; i++)
    {
      if(a.label[i] != b.label[i]) { return (a.label[i] < b.label[i]); }
    }
    return false;
  }
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_SUPPORT_H
