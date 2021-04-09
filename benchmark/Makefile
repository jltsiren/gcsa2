SDSL_DIR=../../sdsl-lite
GCSA_DIR=..
include $(SDSL_DIR)/Make.helper

# Multithreading with OpenMP and libstdc++ Parallel Mode.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) ../libgcsa2.a -lsdsl -ldivsufsort -ldivsufsort64

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # If HOMEBREW_PREFIX is specified, libomp probably cannot be found automatically.
        ifdef HOMEBREW_PREFIX
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        # Macports installs libomp to /opt/local/lib/libomp
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        endif

        # We also need to link it.
        LIBS += -lomp
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(GCSA_DIR)/include -I$(INC_DIR)
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
