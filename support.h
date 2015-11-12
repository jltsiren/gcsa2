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

#include <sdsl/rmq_support.hpp>

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
  const static key_type  CHAR_MASK = 0x7;
  const static size_type MAX_LENGTH = 16;
  const static key_type  PRED_SUCC_MASK = 0xFFFF;

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

  inline static key_type merge(key_type key1, key_type key2) { return (key1 | (key2 & PRED_SUCC_MASK)); }
  inline static key_type replace(key_type key, size_type kmer_val)
  {
    return (kmer_val << 16) | (key & PRED_SUCC_MASK);
  }

  inline static size_type lcp(key_type a, key_type b, size_type kmer_length)
  {
    size_type res = 0;
    key_type mask = CHAR_MASK << (CHAR_WIDTH * kmer_length);
    a = label(a); b = label(b);

    while(mask > 0)
    {
      mask >>= CHAR_WIDTH;
      if((a & mask) != (b & mask)) { break; }
      res++;
    }

    return res;
  }
};

//------------------------------------------------------------------------------

/*
  node_type is the identifier of a node in the input graph. It corresponds to a
  position in the original graph.

  The identifier contains 53 bits for node id, 1 bit for orientation (forward or
  reverse complement, and 10 bits for node offset. The string representation of
  a node_type is id:offset for forward positions and id:-offset for reverse
  complement positions. If the forward offsets are 0 to k, the corresponding
  reverse complement offsets are -k to -0 (in the same order).
*/

typedef std::uint64_t node_type;

struct Node
{
  const static size_type ID_OFFSET        = 11;
  const static size_type ORIENTATION_MASK = 0x400;
  const static size_type OFFSET_MASK      = 0x3FF;

  inline static node_type encode(size_type node_id, size_type node_offset)
  {
    return (node_id << ID_OFFSET) | node_offset;
  }

  inline static node_type encode(size_type node_id, size_type node_offset, bool reverse_complement)
  {
    return encode(node_id, node_offset) | (reverse_complement ? ORIENTATION_MASK : 0);
  }

  static node_type encode(const std::string& token);
  static std::string decode(node_type node);

  inline static size_type id(node_type node) { return node >> ID_OFFSET; }
  inline static bool rc(node_type node) { return node & ORIENTATION_MASK; }
  inline static size_type offset(node_type node) { return node & OFFSET_MASK; }
};

//------------------------------------------------------------------------------

struct KMer
{
  key_type  key;
  node_type from, to;

  KMer();
  KMer(const std::vector<std::string>& tokens, const Alphabet& alpha, size_type successor);

  KMer(key_type _key, node_type _from, node_type _to) :
    key(_key), from(_from), to(_to)
  {
  }

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
  Convention: If a is a proper prefix of b, then
  - a < b, if a is a first label; and
  - a > b, if a is a last label.
*/

struct PathLabel
{
  typedef std::uint32_t rank_type;

  // This should be at least 1 << GCSA::DOUBLING_STEPS.
  const static size_type LABEL_LENGTH = 8;
  const static size_type LENGTH_MASK  = 0xF;
  const static size_type FIRST_MASK   = 0x10; // First label or last label.

  // Labels starting with NO_RANK will be after real labels in lexicographic order.
  const static rank_type NO_RANK = ~(rank_type)0;

  size_type fields;
  rank_type label[LABEL_LENGTH];

  PathLabel() : fields(0) {}

  inline size_type length() const { return this->fields & LENGTH_MASK; }
  inline void setLength(size_type len) { this->fields &= ~LENGTH_MASK; this->fields |= len & LENGTH_MASK; }

  inline bool first() const { return this->fields & FIRST_MASK; }
  inline void setFirst() { this->fields |= FIRST_MASK; }
  inline void setLast()  { this->fields &= ~FIRST_MASK; }

  // Is *this < another?
  inline bool compare(const PathLabel& another) const
  {
    size_type order = std::min(this->length(), another.length());
    for(size_type i = 0; i < order; i++)
    {
      if(this->label[i] != another.label[i]) { return (this->label[i] < another.label[i]); }
    }
    if(this->first()) { return (this->length() < another.length()); }
    else { return (another.length() < this->length()); }
  }
};

//------------------------------------------------------------------------------

/*
  The node type used during doubling. As in the original GCSA, from and to are nodes
  in the original graph, denoting a path as a semiopen range [from, to). If
  from == -1, the path will not be extended, because it already has a unique label.
  rank_type is the integer type used to store ranks of the original kmers.
  During edge generation, 'to' node will be used to store indegree and the outdegree.

  The rank sequences are stored in a separate array at position 'pointer()'. The stored
  sequence consists of the first label ('order()' ranks) followed by one rank for the
  diverging last rank of the last label. If the first and the last ranks are identical,
  the last rank is a dummy value.
*/

struct PathNode
{
  typedef PathLabel::rank_type rank_type;

  const static size_type LABEL_LENGTH = PathLabel::LABEL_LENGTH;

  node_type from, to;

  inline bool sorted() const { return (this->to == ~(node_type)0); }
  inline void makeSorted() { this->to = ~(node_type)0; }

  /*
    We reuse the 'to' field for indegree (upper 2 bits) and outdegree (lower 30 bits).
    The only relevant states for indegree are 0, 1, and 2+, which we encode as 0, 1,
    and 3+.
  */
  inline void initDegree() { this->to = 0; }

  inline size_type outdegree() const { return (this->to & 0x3FFFFFFF); }
  inline void incrementOutdegree() { this->to++; }

  inline size_type indegree() const { return (this->to >> 30); }
  inline void incrementIndegree()
  {
    this->to |= ((this->indegree() << 1) | 1) << 30;
  }

//------------------------------------------------------------------------------

  /*
    From low-order to high-order bits:

    8 bits   which predecessor comp values exist
    4 bits   length of the kmer rank sequences representing the path label range
    4 bits   lcp of the above sequences
    8 bits   unused
    40 bits  pointer to the label data
  */
  size_type fields;

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

  // Order is the length of the kmer rank sequences representing the path label range.
  inline size_type order() const { return ((this->fields >> 8) & 0xF); }
  inline void setOrder(size_type new_order)
  {
    this->fields &= ~(size_type)0xF00;
    this->fields |= new_order << 8;
  }

  // LCP is the length of the common prefix of kmer rank sequences.
  inline size_type lcp() const { return ((this->fields >> 12) & 0xF); }
  inline void setLCP(size_type new_lcp)
  {
    this->fields &= ~(size_type)0xF000;
    this->fields |= new_lcp << 12;
  }

  inline size_type pointer() const { return (this->fields >> 24); }
  inline void setPointer(size_type new_pointer)
  {
    this->fields &= 0xFFFFFF;
    this->fields |= new_pointer << 24;
  }

//------------------------------------------------------------------------------

  inline size_type ranks() const { return this->order() + 1; }
  inline size_type bytes() const { return sizeof(*this) + this->ranks() * sizeof(rank_type); }

  inline PathLabel firstLabel(const std::vector<rank_type>& labels) const
  {
    PathLabel res;
    res.setLength(this->order());
    for(size_type i = 0, j = this->pointer(); i < res.length(); i++, j++)
    {
      res.label[i] = labels[j];
    }
    res.setFirst();
    return res;
  }

  inline PathLabel lastLabel(const std::vector<rank_type>& labels) const
  {
    PathLabel res;

    res.setLength(this->lcp());
    for(size_type i = 0, j = this->pointer(); i < res.length(); i++, j++)
    {
      res.label[i] = labels[j];
    }

    if(this->lcp() < this->order())
    {
      res.label[res.length()] = labels[this->pointer() + this->order()];
      res.setLength(this->lcp() + 1);
    }

    res.setLast();
    return res;
  }

  inline rank_type firstLabel(size_type i, const std::vector<rank_type>& labels) const
  {
    return labels[this->pointer() + i];
  }

  inline rank_type firstLabel(size_type i, const rank_type* labels) const
  {
    return labels[i];
  }

  inline rank_type lastLabel(size_type i, const std::vector<rank_type>& labels) const
  {
    return this->lastLabel(i, labels.data() + this->pointer());
  }

  inline rank_type lastLabel(size_type i, const rank_type* labels) const
  {
    if(i < this->lcp()) { return labels[i]; }
    else { return labels[this->order()]; }
  }

  // Do the two path nodes intersect?
  bool intersect(const PathLabel& first, const PathLabel& last, const std::vector<rank_type>& labels) const;

//------------------------------------------------------------------------------

  static std::vector<rank_type> dummyRankVector();

  PathNode(const KMer& kmer, std::vector<rank_type>& labels);

  PathNode(const PathNode& source,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  PathNode(const PathNode& left, const PathNode& right,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  /*
    Warning: Do not mix the versions using a std::vector and a pointer.
  */
  PathNode(std::istream& in, std::vector<rank_type>& labels);
  PathNode(std::istream& in, rank_type* labels);
  void serialize(std::ostream& out, const std::vector<rank_type>& labels) const;
  void serialize(std::ostream& out, const rank_type* labels) const;

  void print(std::ostream& out, const std::vector<rank_type>& labels) const;

  PathNode();
  explicit PathNode(std::vector<rank_type>& labels);
  PathNode(PathNode&& source);
  ~PathNode();

  inline void swap(PathNode& another)
  {
    if(&another != this)
    {
      std::swap(this->from, another.from); std::swap(this->to, another.to);
      std::swap(this->fields, another.fields);
    }
  }

  PathNode& operator= (PathNode&& source);

  /*
    These are dangerous, because the nodes will share the same label. Changing one will
    change the other as well.
  */
  PathNode(const PathNode& source);
  PathNode& operator= (const PathNode& source);
  void copy(const PathNode& source);
};

// Compares the first labels.
struct PathFirstComparator
{
  const std::vector<PathNode::rank_type>& labels;

  explicit PathFirstComparator(const std::vector<PathNode::rank_type>& _labels) : labels(_labels) { }

  inline bool operator() (const PathNode& a, const PathNode& b) const
  {
    size_type order = std::min(a.order(), b.order());
    for(size_type i = 0, a_ptr = a.pointer(), b_ptr = b.pointer(); i < order; i++, a_ptr++, b_ptr++)
    {
      if(labels[a_ptr] != labels[b_ptr]) { return (labels[a_ptr] < labels[b_ptr]); }
    }
    return (a.order() < b.order());
  }
};

// Compares the 'from' nodes.
struct PathFromComparator
{
  inline bool operator() (const PathNode& a, const PathNode& b) const
  {
    return (a.from < b.from);
  }
};

//------------------------------------------------------------------------------

struct LCP
{
  typedef sdsl::rmq_succinct_sada<>       rmq_type; // Faster than rmq_support_sct.
  typedef PathNode::rank_type             rank_type;
  typedef std::pair<rank_type, rank_type> rank_range;

  size_type           kmer_length, total_keys;
  sdsl::int_vector<0> kmer_lcp;
  rmq_type            lcp_rmq;

  LCP();
  LCP(const std::vector<key_type>& keys, size_type _kmer_length);

  /*
    Computes the minimal/maximal lcp of the path labels corresponding to path nodes a and b.
    a must be before b in lexicographic order, and the ranges must not overlap.
    The returned lcp value is a pair (x,y), where x is the lcp of the PathNode labels
    and y is the lcp of the first diverging kmers.

    FIXME Later: Do not use the rmq if the kmer ranks are close.
  */
  range_type min_lcp(const PathNode& a, const PathNode& b, const std::vector<rank_type>& labels) const;
  range_type max_lcp(const PathNode& a, const PathNode& b, const std::vector<rank_type>& labels) const;

  range_type min_lcp(const PathNode& a, const PathNode& b,
    const rank_type* a_labels, const rank_type* b_labels) const;
  range_type max_lcp(const PathNode& a, const PathNode& b,
    const rank_type* a_labels, const rank_type* b_labels) const;

  // Increments the lcp by 1.
  inline range_type increment(range_type lcp) const
  {
    if(lcp.second + 1 < this->kmer_length) { lcp.second++; }
    else { lcp.first++; lcp.second = 0; }
    return lcp;
  }

  /*
    Extends the given rank range into a maximal range having the given lcp.

    FIXME Later: Build an LCP interval tree to do this faster.
  */
  inline rank_range extendRange(rank_range range, size_type lcp) const
  {
    while(range.first > 0 && this->kmer_lcp[range.first] >= lcp) { range.first--; }
    while(range.second + 1 < this->total_keys && this->kmer_lcp[range.second + 1] >= lcp) { range.second++; }
    return range;
  }

  void swap(LCP& another);
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_SUPPORT_H
