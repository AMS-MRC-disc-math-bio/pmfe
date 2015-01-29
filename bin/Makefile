# Binaries
BINS = pmfe-findmfe pmfe-subopt pmfe-parametrizer pmfe-rnascorer
OLDBINS = parametrizer rnascorer pmfe-param # Old binaries that should be deleted

# include directories
INCLUDES += -I../pmfe/include
INCLUDES += -I../rnascoring/include
INCLUDES += -I../iB4e
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

LIBS += -L../pmfe -lpmfe
LIBS += -L../rnascoring -lrnascoring

# library build dependencies
LIBDEP += ../pmfe/libpmfe.a
LIBDEP += ../rnascoring/librnascoring.a

# Optional OpenMP-related flags
HAS_OPENMP = $(shell $(CXX) -lgomp ../pmfe/openmp-test.cc -o /dev/null 2>&1 > /dev/null ; echo $$?)
ifeq ($(HAS_OPENMP), 0)
CXXFLAGS += -fopenmp
CXXFLAGS += -DHAVE_OPENMP
LIBS += -lgomp
endif

all: $(BINS)

pmfe-findmfe: pmfe_wrapper.o $(LIBDEP)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

pmfe-subopt: psubopt_wrapper.o $(LIBDEP)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

pmfe-parametrizer: parametrizer_wrapper.o NLTemplate.o $(LIBDEP)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

pmfe-rnascorer: rnascorer_wrapper.o $(LIBDEP)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

%.o: %.cc $(LIBDEP)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $^ -o $@ $(LIBS)

clean:
	-rm -vf $(wildcard *.o) $(BINS) $(OLDBINS)

install:
	-install -m 0755 $(BINS) $(PREFIXBIN)

uninstall:
	-rm -vf $(addprefix $(PREFIXBIN), $(BINS) $(OLDBINS))

.PHONY: clean
