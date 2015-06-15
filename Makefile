SDSL_DIR=../sdsl-lite

# In OS X, getrusage() returns maximum resident set size in bytes.
# In Linux, the value is in kilobytes, so this line should be commented out.
#RUSAGE_FLAGS=-DRUSAGE_IN_BYTES

# Multithreading with OpenMP and libstdc++ Parallel Mode. Requires g++ 4.7
# or newer.
PARALLEL_FLAGS=-fopenmp

# Verbose output during index construction etc.
OUTPUT_FLAGS=-DVERBOSE_STATUS_INFO

OTHER_FLAGS=$(RUSAGE_FLAGS) $(PARALLEL_FLAGS) $(OUTPUT_FLAGS)

include $(SDSL_DIR)/Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(OTHER_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR)
LIBOBJS=gcsa.o support.o utils.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
OBJS=$(SOURCES:.cpp=.o)
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64
LIBRARY=gcsa.a
PROGRAMS=build_gcsa

all: $(LIBRARY) $(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

gcsa.a:$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

build_gcsa:build_gcsa.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

clean:
	rm -f $(PROGRAMS) $(OBJS) $(LIBRARY)
