# GCSA2

This is a reimplementation of the Generalized Compressed Suffix Array (GCSA), a BWT-based index for directed graphs. The implementation is based on the [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) (SDSL). To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory. As the implementation uses C++11, OpenMP, and libstdc++ parallel mode, you need g++ 4.7 or newer to compile.

[The old implementation](http://jltsiren.kapsi.fi/gcsa) indexed all paths in a directed acyclic graph, which had to be determinized before index construction. This implementation indexes paths of length at most 128 in any graph. The limit on path length should limit the combinatorial explosion often occurring in graphs containing regions with a lot of branching in a small area.

The input to index construction is a set of paths of length up to *k* in the input graph. The prefix-doubling algorithm transforms the input into an equivalent of order-*8k* de Bruijn graph for the paths of the input graph. As such, the index supports path queries of length up to *8k*. As each doubling step is followed by a pruning step that merges lexicographically adjacent paths starting from the same node, the resulting graph should be smaller than a de Bruijn graph.

At the moment, GCSA2 is being developed as a plugin to Erik Garrison's [variant graph tools](https://github.com/ekg/vg). The only implemented construction option is based on extracting *k*-mers from vg. Later, GCSA2 should become a more general graph indexing library.

See the wiki for some further documentation.

## Compilation options

The maximum resident size reported by `getrusage()` is in kilobytes in Linux and in bytes in OS X. By default, the implementation assumes Linux-like behavior. To get the correct memory usage reports in OS X, uncomment the line `RUSAGE_FLAGS=-DRUSAGE_IN_BYTES` in the makefile.

There are some verification/debugging options in `build_gcsa`. To disable them, comment out the line `VERIFY_FLAGS=-DVERIFY_CONSTRUCTION` in the makefile.

Index construction can be set to output some status information to `stderr` by uncommenting the line `OUTPUT_FLAGS=-DVERBOSE_STATUS_INFO` in the makefile.

## Construction

**Warning:** The optional arguments will eventually change.

The primary GCSA constructor takes up to five parameters:

* `std::vector<KMer>& kmers` contains the graph as a list of paths of a fixed length.
* `size_type kmer_length` is the length of the paths.
* `size_type doubling_steps` (optional; default 3) is the number of doubling steps to perform. The implementation currently supports 1-3 steps.
* `size_type size_limit` (optional; in gigabytes; default 200) limits the size of temporary graphs.
* `const Alphabet& _alpha` (optional) is an object that maps characters to a contiguous range of small integers.

The constructor writes several large temporary files to disk to save memory. After the index has been built, it can be serialized and loaded by using the SDSL `serialize()` / `load()` interface.

Each `KMer` object contains three integers: a key encoding the label of the path and its predecessor/successor characters, the starting node of the path, and a successor node. If the path has multiple successors, separate `KMer` objects must be created for each of them. See issue #6 for partial documentation.

## Interface

Query `find(P)` returns the lexicographic range of path nodes matching pattern *P*. The pattern can be a pair of random-access iterators (e.g. those returned by `pattern.begin()` and `pattern.end()`), a container with random-access iterators (e.g. `std::string`, `std::vector`, or `sdsl::int_vector`), or a character pointer with string length.

Query `locate(i, results)` returns the identifiers of the input nodes at the beginning of the paths represented by the path node with lexicographic rank *i*. Query `locate(range_type(sp,ep), results)` does the same for lexicographic range *[sp,ep]*.

The low-level interface and the graph navigation operations are still subject to change.

## References

Jouni Sirén, Niko Välimäki, and Veli Mäkinen: **Indexing Graphs for Path Queries with Applications in Genome Research**.
IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
[DOI: 10.1109/TCBB.2013.2297101](http://dx.doi.org/10.1109/TCBB.2013.2297101)

Alexander Bowe, Taku Onodera, Kunihiko Sadakane, and Tetsuo Shibuya: **Succinct de Bruijn Graphs**.
Proc. WABI 2012, Springer LNCS 7534, pp. 225-235, Ljubljana, Slovenia, September 10-12, 2012.
[DOI: 10.1007/978-3-642-33122-0_18](http://dx.doi.org/10.1007/978-3-642-33122-0_18)
