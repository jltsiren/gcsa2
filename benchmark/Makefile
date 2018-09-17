SDSL_DIR=../../sdsl-lite
GCSA_DIR=..
include $(SDSL_DIR)/Make.helper

# Multithreading with OpenMP and libstdc++ Parallel Mode.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -L$(GCSA_DIR) -lgcsa2 -lsdsl -ldivsufsort -ldivsufsort64

ifeq ($(shell uname -s),Darwin)
    # Our compiler might be clang that lacks -fopenmp support.
    # Sniff that
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
		# The compiler complained about fopenmp instead of its nonsense input file.
        # We need to use the hard way of getting OpenMP not bundled with the compiler.
        
		# The compiler only needs to do the preprocessing
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

		ifeq ($(shell if [ -d /opt/local/lib/libomp ];then echo 1;else echo 0;fi), 1)
            # Use /opt/local/lib/libomp if present, because Macports installs libomp there.
			# Brew is supposed to put it somewhere the compiler can find it by default.
            LIBS += -L/opt/local/lib/libomp
            # And we need to find the includes. Homebrew puts them in the normal place
            # but Macports hides them in "libomp"
            PARALLEL_FLAGS += -I/opt/local/include/libomp
        endif

        # We also need to link it
        LIBS += -lomp

    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -I$(GCSA_DIR)/include
LIBOBJS=algorithms.o dbg.o files.o gcsa.o internal.o lcp.o path_graph.o support.o utils.o
SOURCES=$(wildcard *.cpp)
OBJS=$(SOURCES:.cpp=.o)

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
