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

#ifndef _GCSA_ALGORITHMS_H
#define _GCSA_ALGORITHMS_H

#include "gcsa.h"
#include "lcp.h"

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
bool verifyIndex(const GCSA& index, const LCPArray* lcp, std::vector<KMer>& kmers, size_type kmer_length);
bool verifyIndex(const GCSA& index, const LCPArray* lcp, const InputGraph& graph);

/*
  Kmer counting. Returns the number of kmers of the given length. If k > order(),
  prints a warning and returns immediately unless force == true.

  A kmer is a path of length k containing comp values 0 <= comp < sigma - 1. Alternatively,
  the path contains either k values 0 < comp < sigma - 1, or k - 1 such values followed
  by a 0.
*/
size_type countKMers(const GCSA& index, size_type k, bool force = false);

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // _GCSA_ALGORITHMS_H
