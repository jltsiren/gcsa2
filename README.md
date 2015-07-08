# GCSA2

This is a reimplementation of the Generalized Compressed Suffix Array (GCSA), a BWT-based index for directed graphs. The implementation is based on the [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) (SDSL). To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory. As the implementation uses C++11, OpenMP, and libstdc++ parallel mode, you need g++ 4.7 or newer to compile.

[The old implementation](http://jltsiren.kapsi.fi/gcsa) indexed all paths in a directed acyclic graph, which had to be determinized before index construction. This implementation indexes paths of length at most 128 in any graph. The limit on path length should limit the combinatorial explosion often occurring in graphs containing regions with a lot of branching in a small area.

The input to index construction is a set of paths of length up to *k* in the original graph. The prefix-doubling algorithm transforms the input into an equivalent of order-*8k* de Bruijn graph for the paths of the input graph. As such, the index supports path queries of length up to *8k*. As each doubling step is followed by a pruning step that merges lexicographically adjacent paths starting from the same node, the resulting graph should be smaller than a de Bruijn graph.

At the moment, GCSA2 is being developed as a plugin to Erik Garrison's [variant graph tools](https://github.com/ekg/vg). The only implemented construction option is based on extracting *k*-mers from vg. Later, GCSA2 should become a more general graph indexing library.

## Compilation options

The maximum resident size reported by `getrusage()` is in kilobytes in Linux and in bytes in OS X. By default, the implementation assumes Linux-like behavior. To get the correct memory usage reports in OS X, uncomment the line `RUSAGE_FLAGS=-DRUSAGE_IN_BYTES` in the makefile.

There are some verification/debugging options in `build_gcsa`. To disable them, comment out the line `VERIFY_FLAGS=-DVERIFY_CONSTRUCTION` in the makefile.

Index construction can be set to output some status information to `stderr` by uncommenting the line `OUTPUT_FLAGS=-DVERBOSE_STATUS_INFO` in the makefile.

## Data model

The input to GCSA2 is a directed graph. Each **node** of the input graph is a pair *(id,c)*, where integer *id* is the **unique identifier** of the node and character *c* is the **label** of the node. For best results, nodes on unary paths should have successive identifiers. Each node must have at least one incoming edge and one outgoing edge.

In the current implementation, the graph must have exactly one **source** node and one **sink** node. There must be an edge from the sink node to the source node. The source node must not have any other incoming edges, and the sink node must not have any other outgoing edges. The source and the sink must have unique labels, which the `Alphabet` object must map to values that are smaller than larger, respectively, than for any real character. These additional restrictions are a matter of convenience. There are no fundamental or performance reasons for having them.

The nodes of the final transformed graph are called **path nodes**, as they correspond to sets of **paths** in the original graph. Path nodes are identified by their ranks in lexicographic order. All paths having the same label are represented by the same path node. A path node matches **pattern** *P*, if either *P* is a prefix of its label, or the corresponding path in the input graph can be extended to match the pattern. (The construction guarantees that all paths represented by the same path node have the same extensions up to the maximum query length.)

## Construction

The primary GCSA constructor takes up to four parameters:

* `std::vector<KMer>& kmers` contains the graph as a list of paths of a fixed length.
* `size_type kmer_length` is the length of the paths.
* `size_type doubling_steps` (optional; default 3) is the number of doubling steps to perform. The implementation currently supports 1-3 steps.
* `const Alphabet& _alpha` (optional) is an object that maps characters to a contiguous range of small integers.

After the index has been built, it can be serialized and loaded by using the SDSL `serialize()` / `load()` interface.

Each `KMer` object contains three integers: a key encoding the label of the path and its predecessor/successor characters, the starting node of the path, and a successor node. If the path has multiple successors, separate `KMer` objects must be created for each of them. The `KMer` construction interface has not been finalized yet.

## Interface

Query `find(P)` returns the lexicographic range of path nodes matching pattern *P*. The pattern can be a pair of random-access iterators (e.g. those returned by `pattern.begin()` and `pattern.end()`), a container with random-access iterators (e.g. `std::string`, `std::vector`, or `sdsl::int_vector`), or a character pointer with string length.

Query `locate(i, results)` returns the identifiers of the input nodes at the beginning of the paths represented by the path node with lexicographic rank *i*. Query `locate(range_type(sp,ep), results)` does the same for lexicographic range *[sp,ep]*.

The construction interface, the low-level interface and the graph navigation operations are still subject to change.

## Version history

### Current version

* More space-efficient construction.
* Support for cyclic graphs.
* Binary format for input graphs.
* Fixed a rare incompatibility between `mergePaths()` and `mergeByLabel()`.

### 0.1 (2015-05-27)

* The first release.
* Index construction from paths extracted from vg.
* `find()` and `locate()` queries.

## Todo

* Optimizations
  * Multi-threaded construction.
  * More space-efficient construction. The size of a `PathNode` can be reduced to 16 bytes with a more efficient encoding for *(id,offset)* pairs. We can write the results of a doubling step to disk, and even do the merging step on disk if necessary.
  * More space-efficient index representation.
  * Sample compression. More efficient encoding for the *(id,offset)* pairs would already help.
  * Determine when nodes have to be prefix-sorted and when we can make them prefix-range-sorted.
* Generalizations
  * Make index construction work with 0 doubling steps.
  * Support for larger alphabets.
  * Support for inputs other than vg paths.
  * Sampling when starting positions are just node ids without offsets.
* Interface
  * Low-level interface.
  * Graph navigation operations.
* A paper describing the new algorithmic ideas.

## References

Jouni Sirén, Niko Välimäki, and Veli Mäkinen: **Indexing Graphs for Path Queries with Applications in Genome Research**.
IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
[DOI: 10.1109/TCBB.2013.2297101](http://dx.doi.org/10.1109/TCBB.2013.2297101)

Alexander Bowe, Taku Onodera, Kunihiko Sadakane, and Tetsuo Shibuya: **Succinct de Bruijn Graphs**.
Proc. WABI 2012, Springer LNCS 7534, pp. 225-235, Ljubljana, Slovenia, September 10-12, 2012.
[DOI: 10.1007/978-3-642-33122-0_18](http://dx.doi.org/10.1007/978-3-642-33122-0_18)
