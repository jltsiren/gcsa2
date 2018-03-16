SDSL_DIR=../../sdsl-lite
GCSA_DIR=..

# Multithreading with OpenMP and libstdc++ Parallel Mode.
PARALLEL_FLAGS=-fopenmp -pthread

include $(SDSL_DIR)/Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -I$(GCSA_DIR)/include
LIBOBJS=algorithms.o dbg.o files.o gcsa.o internal.o lcp.o path_graph.o support.o utils.o
SOURCES=$(wildcard *.cpp)
OBJS=$(SOURCES:.cpp=.o)
LIBS=-L$(LIB_DIR) -L$(GCSA_DIR) -lgcsa2 -lsdsl -ldivsufsort -ldivsufsort64
PROGRAMS=count_kmers query_gcsa csa_builder csa_query

all: $(PROGRAMS)

%.o:%.cpp
	$(MY_CXX) $(CXX_FLAGS) -c $<

count_kmers:count_kmers.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

query_gcsa:query_gcsa.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

csa_builder:csa_builder.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

csa_query:csa_query.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(PROGRAMS) $(OBJS)
