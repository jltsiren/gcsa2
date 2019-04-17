/*
  Copyright (c) 2019 Jouni Sirén
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

#ifndef GCSA_ALGORITHMS_H
#define GCSA_ALGORITHMS_H

#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>

namespace gcsa
{

/*
  algorithms.h: Algorithms using GCSA. All of them use multiple threads.
*/

//------------------------------------------------------------------------------

/*
  Index verification. The index is queried with all unique kmer labels in the input, and
  the list of occurrences is expected to be the same as the set of start nodes of the
  kmers with that label. To guarantee this, the input should be the same as for the
  constructor, or a subset where all kmers with a given label are either present or
  absent. Counting queries are also verified.

  If lcp != 0, parent() queries are also verified.

  Note: Verification sorts the kmer array.

  Returns false if index verification fails, and true otherwise.
*/
bool verifyIndex(const GCSA& index, const LCPArray* lcp, std::vector<KMer>& kmers, size_type kmer_length, const NodeMapping& mapping = NodeMapping());
bool verifyIndex(const GCSA& index, const LCPArray* lcp, const InputGraph& graph);

//------------------------------------------------------------------------------

struct KMerSearchParameters
{
  size_type seed_length;  // Parallelize using seed kmers of this length.
  bool include_Ns;        // Include kmers containing Ns in the search.
  bool force;             // Force searching for kmers longer than the order of the index.
  std::string output;     // Base name for output.

  constexpr static size_type SEED_LENGTH = 5;
  const static std::string LEFT_EXTENSION;  // .left
  const static std::string RIGHT_EXTENSION; // .right

  KMerSearchParameters() : seed_length(SEED_LENGTH), include_Ns(false), force(false), output() {}
};

/*
  Kmer counting. Returns the number of kmers of the given length. If k > order(),
  prints a warning and returns immediately unless force == true.

  By default, the algorithm counts only kmers consisting of bases (comp values 1-4;
  comp values encoded in fast_bwt). If include_Ns == true, the algorithm will also
  count kmers containing Ns (comp values encoded in sparse_bwt, except 0 and sigma - 1).
*/
size_type countKMers(const GCSA& index, size_type k,
  const KMerSearchParameters& parameters = KMerSearchParameters());

/*
  As above, but counts the kmers in two indexes and reports the results as
  (shared kmers, kmers unique to left, kmers unique to right). The indexes must
  use compatible alphabets.
*/
std::array<size_type, 3> compareKMers(const GCSA& left, const GCSA& right, size_type k,
  const KMerSearchParameters& parameters = KMerSearchParameters());

//------------------------------------------------------------------------------

void printStatistics(const GCSA& gcsa, const LCPArray& lcp_array);

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_ALGORITHMS_H
