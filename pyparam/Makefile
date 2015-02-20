# include directories
INCLUDES += -Iinclude
INCLUDES += -I../pmfe/include
INCLUDES += -I../iB4e
INCLUDES += -I/usr/include/python2.7

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O2

# library paths
LIBS += -L../pmfe -lpmfe
LIBS += -L../rnascoring -lrnascoring
LIBS += -lgmp
LIBS += -lCGAL
LIBS += -lboost_python
LIBS += -lpython2.7

# Optional OpenMP-related flags
HAS_OPENMP = $(shell $(CXX) -lgomp ../pmfe/openmp-test.cc -o /dev/null 2>&1 > /dev/null ; echo $$?)
ifeq ($(HAS_OPENMP), 0)
CXXFLAGS += -fopenmp
CXXFLAGS += -DHAVE_OPENMP
LIBS += -lgomp
endif

#LIBS += -shared

LIBRARY = pyparam_base.so

all: $(LIBRARY)

%.so: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -shared -o $@ $(LIBS)


clean:
	-rm -vf $(LIBRARY)

install:


.PHONY: clean
