# include directories
INCLUDES += -I../gtmfe/include
INCLUDES += -I../rnascoring/include
INCLUDES += -I../iB4e/include
INCLUDES += -I/usr/include/python2.7

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O2

# library paths
LIBS += -lgmp
LIBS += -lCGAL
LIBS += -lm
LIBS += -lboost_python
LIBS += -lpython2.7
LIBS += -lboost_system
LIBS += -lboost_filesystem
LIBS += -lboost_program_options

LIBS += -L../gtmfe -lgtmfe
LIBS += -L../iB4e -liB4e
LIBS += -L../rnascoring -lrnascoring

#LDFLAGS += -static

# Optional OpenMP-related flags
HAS_OPENMP = $(shell $(CXX) -lgomp ../gtmfe/openmp-test.cc -o /dev/null 2>&1 > /dev/null ; echo $$?)
ifeq ($(HAS_OPENMP), 0)
CXXFLAGS += -fopenmp
CXXFLAGS += -DHAVE_OPENMP
LIBS += -lgomp
endif

all: gtmfe-param parametrizer rnascorer

gtmfe-param: gtmfe_wrapper.o
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

parametrizer: parametrizer_wrapper.o
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS) $(LIBS) # Some sort of ld glitch requires doubling $(LIBS) on my machine-AGD

rnascorer: rnascorer_wrapper.o
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $^ -o $@ $(LIBS)

clean:
	-rm -vf $(wildcard *.o) gtmfe-param parametrizer

.PHONY: clean
