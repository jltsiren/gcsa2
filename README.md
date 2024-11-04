# GCSA2

This is a reimplementation of the Generalized Compressed Suffix Array (GCSA), a BWT-based index for directed graphs.

See [the wiki](https://github.com/jltsiren/gcsa2/wiki) for further documentation.

## Overview

[The old implementation](https://jltsiren.kapsi.fi/gcsa) indexed all paths in a directed acyclic graph, which had to be determinized before index construction. GCSA2 approximates the graph with a de Bruijn graph (with *k = 256* or less), allowing it to index more complex graphs as well as graphs containing cycles. The order of the de Bruijn graph limits the length of the queries the index is able to answer correctly. Longer queries may result in false positives (but no false negatives).

The input to index construction is a set of paths of length *k* in the input graph. The prefix-doubling algorithm transforms the input into an order-*8k* (order-*2k*, order-*4k*, order-*16k*) pruned de Bruijn graph for paths in the input graph. A pruned de Bruijn graph differs from a de Bruijn graph in that its nodes may have shorter labels than the order of the graph, if the shorter labels uniquely determine the start nodes of the corresponding paths in the input graph. As such, pruned de Bruijn graphs are usually smaller than proper de Bruijn graphs.

GCSA2 is being developed as a part of [vg](https://github.com/vgteam/vg). The only implemented construction option is based on extracting *k*-mers from vg.

## Compiling GCSA2

This library is designed to take the compiler options from the [vgteam fork](https://github.com/vgteam/sdsl-lite) of the Succinct Data Structures Library 2.0 (SDSL). It currently requires a recent C++ compiler supporting C++17 and OpenMP. GCC is recommended, as the multithreaded `std::sort` the from libstdc++ parallel mode significantly speeds up index construction. On Apple systems, GCSA2 can be built with Apple Clang, but libomp must be installed via Macports or Homebrew.

GCSA2 is frequently tested in the following environments:

* Intel Linux (Ubuntu) with GCC.
* Intel macOS with GCC and Apple Clang.
* ARM macOS with Apple Clang.

To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory (the default is `../sdsl-lite`). Use `make` to compile the library or `install.sh` to compile it and install the headers and the library to your home directory. Another installation directory can be specified as `install.sh prefix`.

## Citing GCSA2

Jouni Sirén: **Indexing Variation Graphs**.
Proc. ALENEX 2017, SIAM, pp. 13-27, Barcelona, Spain, January 17-18, 2017.
[DOI: 10.1137/1.9781611974768.2](https://doi.org/10.1137/1.9781611974768.2)

## Other references

Jouni Sirén, Niko Välimäki, and Veli Mäkinen: **Indexing Graphs for Path Queries with Applications in Genome Research**.
IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
[DOI: 10.1109/TCBB.2013.2297101](https://doi.org/10.1109/TCBB.2013.2297101)

Alexander Bowe, Taku Onodera, Kunihiko Sadakane, and Tetsuo Shibuya: **Succinct de Bruijn Graphs**.
Proc. WABI 2012, Springer LNCS 7534, pp. 225-235, Ljubljana, Slovenia, September 10-12, 2012.
[DOI: 10.1007/978-3-642-33122-0_18](https://doi.org/10.1007/978-3-642-33122-0_18)

Jouni Sirén, Erik Garrison, Adam M. Novak, Benedict Paten, and Richard Durbin: **Haplotype-aware graph indexes**.
Bioinformatics 36(2):400-407, 2020.
DOI: [10.1093/bioinformatics/btz575](https://doi.org/10.1093/bioinformatics/btz575)
