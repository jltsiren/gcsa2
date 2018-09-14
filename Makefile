SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

# This enables various debugging options in build_gcsa.
#VERIFY_FLAGS=-DVERIFY_CONSTRUCTION

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64

ifeq ($(shell uname -s),Darwin)
    # Our compiler might be clang that lacks -fopenmp support.
    # Sniff that
    ifeq ($(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l), 1)
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

OTHER_FLAGS=$(VERIFY_FLAGS) $(PARALLEL_FLAGS)

CXX_FLAGS=$(MY_CXX_FLAGS) $(OTHER_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -Iinclude
LIBOBJS=algorithms.o dbg.o files.o gcsa.o internal.o lcp.o path_graph.o support.o utils.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard include/gcsa/*.h)
OBJS=$(SOURCES:.cpp=.o)

LIBRARY=libgcsa2.a
PROGRAMS=build_gcsa convert_graph gcsa_format try_extend

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

try_extend:try_extend.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

clean:
	rm -f $(PROGRAMS) $(OBJS) $(LIBRARY)
