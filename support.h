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
    - 16x3 bits for the label, with high-order characters 0s when necessary
    - 8 bits for marking which predecessors are present
    - 8 bits for marking which successors are present
*/

typedef std::uint64_t key_type;

struct Key
{
  const static size_type CHAR_WIDTH = 3;
  const static size_type CHAR_MASK = 0x7;
  const static size_type MAX_LENGTH = 16;

  inline static key_type encode(const Alphabet& alpha, const std::string& kmer,
    byte_type pred, byte_type succ)
  {
    key_type value = 0;
    for(size_type i = 0; i < kmer.length(); i++)
    {
      value = (value << CHAR_WIDTH) | alpha.char2comp[kmer[i]];
    }
    value = (value << 8) | pred;
    value = (value << 8) | succ;
    return value;
  }

  static std::string decode(key_type key, size_type kmer_length, const Alphabet& alpha);

  inline static size_type label(key_type key) { return (key >> 16); }
  inline static byte_type predecessors(key_type key) { return (key >> 8) & 0xFF; }
  inline static byte_type successors(key_type key) { return key & 0xFF; }
  inline static comp_type last(key_type key) { return (key >> 16) & CHAR_MASK; }

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
  const static size_type OFFSET_BITS = 10;
  const static size_type OFFSET_MASK = 0x3FF;

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
    return (Key::label(this->key) < Key::label(another.key));
  }

  inline bool sorted() const { return (this->to == ~(node_type)0); }
  inline void makeSorted() { this->to = ~(node_type)0; }

  static byte_type chars(const std::string& token, const Alphabet& alpha);
};

std::ostream& operator<< (std::ostream& out, const KMer& kmer);

inline bool
operator< (key_type key, const KMer& kmer)
{
  return (Key::label(key) < Key::label(kmer.key));
}

/*
  This function does several things:

  1. Sorts the kmer array by the kmer labels encoded in the key
  2. Builds an array of unique kmer labels, with the predecessor and successor
     fields merged from the original kmers.
  3. Stores the last character of each unique kmer label in an array.
  4. Replaces the kmer labels in the keys by their ranks.
*/
void uniqueKeys(std::vector<KMer>& kmers, std::vector<key_type>& keys, sdsl::int_vector<0>& last_char,
  bool print = false);

//------------------------------------------------------------------------------

/*
  The node type used during doubling. As in the original GCSA, from and to are nodes
  in the original graph, denoting a path as a semiopen range [from, to). If
  from == -1, the path will not be extended, because it already has a unique label.
  rank_type is the integer type used to store ranks of the original kmers.
  During edge generation, to will be used to store indegree and the outdegree.

  FIXME Later: Store the labels of all PathNodes in a single array. 4+4 bits in the
  fields tell the actual lengths of the first/last labels, while another 4 bits tell
  the lcp of the two labels.
*/

struct PathNode
{
  typedef std::uint32_t rank_type;

  // This should be at least 1 << GCSA::DOUBLING_STEPS.
  const static size_type LABEL_LENGTH = 8;

  const static rank_type LOWER_PADDING = 0;
  const static rank_type UPPER_PADDING = ~(rank_type)0;

  node_type from, to;
  rank_type first_label[LABEL_LENGTH];
  rank_type last_label[LABEL_LENGTH];


  /*
    From low-order to high-order bits:

    8 bits   which predecessor comp values exist
    4 bits   length of the labels
    52 bits  unused
  */
  size_type fields;

  inline bool sorted() const { return (this->to == ~(node_type)0); }
  inline void makeSorted() { this->to = ~(node_type)0; }

  inline byte_type predecessors() const { return (this->fields & 0xFF); }
  inline void setPredecessors(byte_type preds)
  {
    this->fields &= ~(size_type)0xFF;
    this->fields |= (size_type)preds;
  }
  inline bool hasPredecessor(comp_type comp) const
  {
    return (this->fields & (1 << comp));
  }
  inline void addPredecessors(const PathNode& another)
  {
    this->fields |= another.predecessors();
  }

  // Convention: The labels contain padding beyond the order of the path node.
  inline size_type order() const { return ((this->fields >> 8) & 0xF); }
  inline void setOrder(size_type new_order)
  {
    this->fields &= ~(size_type)0xF00;
    this->fields |= new_order << 8;
  }

  /*
    We reuse the to field for indegree (upper 32 bits) and outdegree (lower 32 bits).
  */
  inline void initDegree() { this->to = 0; }
  inline void incrementOutdegree() { this->to++; }
  inline size_type outdegree() const { return (this->to & 0xFFFFFFFF); }
  inline void incrementIndegree() { this->to += ((size_type)1) << 32; }
  inline size_type indegree() const { return (this->to >> 32); }

//------------------------------------------------------------------------------

  /*
    Convention: If a.first_label is a proper prefix of b.first_label, a < b. If
    a.last_label is a proper prefix of b.last_label, a.last_label > b.last_label.
  */

  // Make this the union of this and another.
  void merge(const PathNode& another);

  // Do the two path nodes intersect?
  bool intersect(const PathNode& another) const;

  inline bool operator< (const PathNode& another) const
  {
    size_type ord = std::min(this->order(), another.order());
    for(size_type i = 0; i < ord; i++)
    {
      if(this->first_label[i] != another.first_label[i])
      {
        return (this->first_label[i] < another.first_label[i]);
      }
    }
    return (this->order() < another.order());
  }

  // Like operator<, but for the last labels.
  inline bool compareLast(const PathNode& another) const
  {
    size_type ord = std::min(this->order(), another.order());
    for(size_type i = 0; i < ord; i++)
    {
      if(this->last_label[i] != another.last_label[i])
      {
        return (this->last_label[i] < another.last_label[i]);
      }
    }
    return (another.order() < this->order());
  }

//------------------------------------------------------------------------------

  explicit PathNode(const KMer& kmer);
  PathNode(const PathNode& left, const PathNode& right);
  explicit PathNode(std::ifstream& in);

  size_type serialize(std::ostream& out) const;

  PathNode();
  PathNode(const PathNode& source);
  PathNode(PathNode&& source);
  ~PathNode();

  void copy(const PathNode& source);
  void swap(PathNode& another);

  PathNode& operator= (const PathNode& source);
  PathNode& operator= (PathNode&& source);
};

struct PathFromComparator
{
  inline bool operator() (const PathNode& a, const PathNode& b) const
  {
    return (a.from < b.from);
  }
};

std::ostream& operator<< (std::ostream& stream, const PathNode& pn);

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_SUPPORT_H
