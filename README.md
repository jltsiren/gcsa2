# GCSA2

This is a reimplementation of the Generalized Compressed Suffix Array (GCSA), a BWT-based index for directed graphs. The implementation is based on the [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) (SDSL). To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory. As the implementation uses C++11, OpenMP, and libstdc++ parallel mode, you need g++ 4.7 or newer to compile.

[The old implementation](http://jltsiren.kapsi.fi/gcsa) is still available. This new implementation should be faster. Note that while the original GCSA was a full index, this implementation stops the prefix-doubling algorithm after three steps, when path length is at most 128. This should be enough to avoid the combinatorial explosion occurring in many graphs representing genetic variation, without resorting to heuristics.

The input to index construction is a set of paths of length up to *k* in the original graph. The prefix-doubling algorithm transforms the input into an equivalent of order-*8k* de Bruijn graph for the paths of the input graph. As such, the index supports path queries of length up to *8k*. As each doubling step is followed by a pruning step that merges lexicographically adjacent paths starting from the same node, the resulting graph should be smaller than a de Bruijn graph.

At the moment, GCSA2 is being developed as a plugin to Erik Garrison's [variant graph tools](https://github.com/ekg/vg). The only implemented construction option is based on extracting *k*-mers from vg. Later, GCSA2 should become a more general graph indexing library.

## Compilation options

The maximum resident size reported by `getrusage()` is in kilobytes in Linux and in bytes in OS X. By default, the implementation assumes Linux-like behavior. To get the correct memory usage reports in OS X, uncomment the line `RUSAGE_FLAGS=-DRUSAGE_IN_BYTES` in the makefile.

There are some verification/debugging options in `build_gcsa`. To disable them, comment out the line `VERIFY_FLAGS=-DVERIFY_CONSTRUCTION` in the makefile.

Index construction can be set to output some status information to `stderr` by uncommenting the line `OUTPUT_FLAGS=-DVERBOSE_STATUS_INFO` in the makefile.

## Data model

The input to GCSA2 is a directed graph. Each **node** of the input graph is a pair *(id,c)*, where integer *id* is the unique identifier of the node and character *c* is the label of the node. For best results, nodes on unary paths should have successive identifiers. At the moment, GCSA2 assumes that the input is a directed acyclic graph. Cyclic graphs will be supported later. The graph must have exactly one **source** node with indegree 0, and exactly one **sink** node with outdegree 0. The sink node must have a unique label that the `Alphabet` object maps into value 0.

The nodes of the final transformed graph are called **path nodes**, as they correspond to sets of **paths** in the original graph. Path nodes are identified by their ranks in lexicographic order. All paths having the same label are represented by the same path node. A path node matches **pattern** *P*, if either *P* is a prefix of its label, or the corresponding path in the input graph can be extended to match the pattern. (The construction guarantees that all paths represented by the same path node have the same extensions up to the maximum query length.)

## Interface

Query `find(P)` returns the lexicographic range of path nodes matching pattern *P*. The pattern can be a pair of random-access iterators (e.g. those returned by `pattern.begin()` and `pattern.end()`), a container with random-access iterators (e.g. `std::string`, `std::vector`, or `sdsl::int_vector`), or a character pointer with string length.

Query `locate(i, results)` returns the identifiers of the input nodes at the beginning of the paths represented by the path node with lexicographic rank *i*. Query `locate(range_type(sp,ep), results)` does the same for lexicographic range *[sp,ep]*.

The construction interface, the low-level interface and the graph navigation operations are still subject to change.

## Version history

### Current version

* More space-efficient construction.
* Support for cyclic graphs.
* Binary format for input graphs.

### 0.1 (2015-05-27)

* The first release.
* Index construction from paths extracted from vg.
* `find()` and `locate()` queries.

## Todo

* Make index construction work with less than 2 doubling steps.
* Support for larger alphabets.
* Multi-threaded construction.
* More space-efficient construction.
* More space-efficient index representation.
* Sample compression (if necessary).
* Support for inputs other than vg paths.
* Sampling when starting positions are just node ids without offsets.
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
