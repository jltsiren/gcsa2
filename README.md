# GCSA2

This is a reimplementation of the Generalized Compressed Suffix Array (GCSA), a BWT-based index for directed graphs. The implementation is based on the [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) (SDSL). To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory. As the implementation uses C++11, OpenMP, and libstdc++ parallel mode, you need g++ 4.7 or newer to compile.

[The old implementation](http://jltsiren.kapsi.fi/gcsa) indexed all paths in a directed acyclic graph, which had to be determinized before index construction. This implementation indexes paths of length up to 128 in any graph. The upper bound on path length should limit the combinatorial explosion often occurring in graphs containing regions with a lot of branching.

The input to index construction is a set of paths of length *k* in the input graph. The prefix-doubling algorithm transforms the input into an equivalent of order-*8k* (order-*2k*, order-*4k*) de Bruijn graph for paths in the input graph. As such, the index supports path queries of length up to the order of the corresponding de Bruijn graph. As each `joinPaths()` step is followed by a `mergePaths()` step that merges lexicographically adjacent paths starting from the same node, the resulting graph should be smaller than a de Bruijn graph.

At the moment, GCSA2 is being developed as a plugin to Erik Garrison's [variant graph tools](https://github.com/ekg/vg). The only implemented construction option is based on extracting *k*-mers from vg. Later, GCSA2 should become a more general graph indexing library.

See [the wiki](https://github.com/jltsiren/gcsa2/wiki) for further documentation.

## Compilation options

The maximum resident size reported by `getrusage()` is in kilobytes in Linux and in bytes in OS X. By default, the implementation assumes Linux-like behavior. To get the correct memory usage reports in OS X, uncomment the line `RUSAGE_FLAGS=-DRUSAGE_IN_BYTES` in the makefile.

There are some verification/debugging options in `build_gcsa`. To disable them, comment out the line `VERIFY_FLAGS=-DVERIFY_CONSTRUCTION` in the makefile.

Index construction can be set to output some status information to `stderr` by uncommenting the line `OUTPUT_FLAGS=-DVERBOSE_STATUS_INFO` in the makefile.

## References

Jouni Sirén, Niko Välimäki, and Veli Mäkinen: **Indexing Graphs for Path Queries with Applications in Genome Research**.
IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
[DOI: 10.1109/TCBB.2013.2297101](http://dx.doi.org/10.1109/TCBB.2013.2297101)

Alexander Bowe, Taku Onodera, Kunihiko Sadakane, and Tetsuo Shibuya: **Succinct de Bruijn Graphs**.
Proc. WABI 2012, Springer LNCS 7534, pp. 225-235, Ljubljana, Slovenia, September 10-12, 2012.
[DOI: 10.1007/978-3-642-33122-0_18](http://dx.doi.org/10.1007/978-3-642-33122-0_18)
