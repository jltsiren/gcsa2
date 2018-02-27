SDSL_DIR=../sdsl-lite

# This enables various debugging options in build_gcsa.
#VERIFY_FLAGS=-DVERIFY_CONSTRUCTION

# Multithreading with OpenMP and libstdc++ Parallel Mode.
PARALLEL_FLAGS=-fopenmp -pthread
# Turn off libstdc++ parallel mode for clang
ifneq (clang,$(findstring clang,$(shell $(CXX) --version)))
PARALLEL_FLAGS+=-D_GLIBCXX_PARALLEL
endif

OTHER_FLAGS=$(VERIFY_FLAGS) $(PARALLEL_FLAGS)

include $(SDSL_DIR)/Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(OTHER_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -Iinclude
LIBOBJS=algorithms.o dbg.o files.o gcsa.o internal.o lcp.o path_graph.o support.o utils.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard include/gcsa/*.h)
OBJS=$(SOURCES:.cpp=.o)
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64
LIBRARY=libgcsa2.a
PROGRAMS=build_gcsa convert_graph gcsa_format

all: $(LIBRARY) $(PROGRAMS)

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

clean:
	rm -f $(PROGRAMS) $(OBJS) $(LIBRARY)
