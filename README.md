# gcsa2

This is going to be a reimplementation of the Generalized Compressed Suffix Array (GCSA). [The old implementation](http://jltsiren.kapsi.fi/gcsa) is based on the [RLCSA library](http://jltsiren.kapsi.fi/rlcsa), which is quite slow by today's standards, and not too pleasant to work with. The new implementation will use the [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite).

GCSA construction often requires heuristics to avoid combinatorial explosions. When working with graphs arising from human genetic variation, this usually happens at doubling step 8, when the length of the paths increases from 128 to 256. The necessity for (input-specific?) heuristics makes GCSA less useful than it could be.

Bowe et al. used a related structure to index de Bruijn graphs in small space. GCSA can be easily adapted to the task. By indexing all k-mers of a graph, we can still support path queries of length up to k, while index construction will be much easier. A naive de Bruijn graph of the k-mers in a variant graph would be too large for interesting values of k. Fortunately the experience from full GCSA construction suggests that merging redundant subgraphs will keep the size of the graph under control for k up to 128.

## References

Jouni Sirén, Niko Välimäki, and Veli Mäkinen: **Indexing Graphs for Path Queries with Applications in Genome Research**.
IEEE/ACM Transactions on Computational Biology and Bioinformatics 11(2):375-388, 2014.
[DOI: 10.1109/TCBB.2013.2297101](http://dx.doi.org/10.1109/TCBB.2013.2297101)

Alexander Bowe, Taku Onodera, Kunihiko Sadakane, and Tetsuo Shibuya: **Succinct de Bruijn Graphs**.
Proc. WABI 2012, Springer LNCS 7534, pp. 225-235, Ljubljana, Slovenia, September 10-12, 2012.
[DOI: 10.1007/978-3-642-33122-0_18](http://dx.doi.org/10.1007/978-3-642-33122-0_18)