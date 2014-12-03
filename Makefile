# source files
SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.d)
HDR = $(wildcard include/*.h)

OBJ-GTMFE = $(filter-out src/iB4e_wrapper.o, $(OBJ))
OBJ-PARAM = $(filter-out src/gtmfe_wrapper.o, $(OBJ))
EXEC = gtmfe-param parametrizer

# include directories
INCLUDES += -Iinclude -Irna-scoring/

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
#CXXFLAGS += -O2
CXXFLAGS += -O0

# library paths
LIBS += -lCGAL
LIBS += -lgmp
LIBS += -lm
LIBS += -lboost_system
LIBS += -lboost_filesystem
LIBS += -lboost_program_options
LIBS += -L./rna-scoring -Wl,-R./rna-scoring -lrnascoring

# Optional OpenMP-related flags
HAS_OPENMP = $(shell $(CXX) -lgomp test.cc 2>&1 > /dev/null ; echo $$?)
$(info Exit code: $(HAS_OPENMP))
ifeq ($(HAS_OPENMP), 0)
#	$(info OpenMP detected!)
	CXXFLAGS += -fopenmp
	CXXFLAGS += -DHAVE_OPENMP
	LIBS += -lgomp
endif

all: rnascoring $(OBJ) $(EXEC)

gtmfe-param: $(OBJ-GTMFE)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(OBJ-GTMFE) -o $@ $(LIBS)

parametrizer: $(OBJ-PARAM)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(OBJ-PARAM) -o $@ $(LIBS)

-include $(DEP)

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $*.cc -o $*.o $(LIBS)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $*.cc $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

rnascoring:
	$(MAKE) -C rna-scoring

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP)

.PHONY: clean
