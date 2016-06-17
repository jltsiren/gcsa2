SDSL_DIR=../sdsl-lite

# In OS X, getrusage() returns maximum resident set size in bytes.
# In Linux, the value is in kilobytes, so this line should be commented out.
#RUSAGE_FLAGS=-DRUSAGE_IN_BYTES

# This enables various debugging options in build_gcsa.
#VERIFY_FLAGS=-DVERIFY_CONSTRUCTION

# Multithreading with OpenMP and libstdc++ Parallel Mode.
PARALLEL_FLAGS=-fopenmp -pthread
# Turn off libstdc++ parallel mode for clang
ifneq (clang,$(findstring clang,$(shell $(CXX) --version)))
PARALLEL_FLAGS+=-D_GLIBCXX_PARALLEL
endif

OTHER_FLAGS=$(RUSAGE_FLAGS) $(VERIFY_FLAGS) $(PARALLEL_FLAGS)

include $(SDSL_DIR)/Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(OTHER_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR)
LIBOBJS=algorithms.o dbg.o files.o gcsa.o internal.o lcp.o path_graph.o support.o utils.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
OBJS=$(SOURCES:.cpp=.o)
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64
LIBRARY=libgcsa2.a
PROGRAMS=build_gcsa convert_graph count_kmers gcsa_format query_gcsa
CSA_BENCHMARK=csa_builder csa_query
OBSOLETE=convert_gcsa

all: $(LIBRARY) $(PROGRAMS)
csa: $(CSA_BENCHMARK)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

build_gcsa:build_gcsa.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

gcsa_format:gcsa_format.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

convert_graph:convert_graph.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

count_kmers:count_kmers.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

query_gcsa:query_gcsa.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

csa_builder:csa_builder.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

csa_query:csa_query.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

clean:
	rm -f $(PROGRAMS) $(OBJS) $(LIBRARY) $(CSA_BENCHMARK) $(OBSOLETE)
