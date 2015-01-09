# source files
SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.d)
HDR = $(wildcard src/*.h)

LIBRARY = libgtmfe.a

# include directories
INCLUDES += -Iinclude
INCLUDES += -I/usr/include/python2.7
INCLUDES += -I../rnascoring/include

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O2

# library paths
LIBS += -L../rnascoring -lrnascoring
LIBS += -lgmp
LIBS += -lm
LIBS += -lboost_python
LIBS += -lpython2.7

# Optional OpenMP-related flags
HAS_OPENMP = $(shell $(CXX) -lgomp openmp-test.cc -o /dev/null 2>&1 > /dev/null ; echo $$?)
ifeq ($(HAS_OPENMP), 0)
CXXFLAGS += -fopenmp
CXXFLAGS += -DHAVE_OPENMP
LIBS += -lgomp
endif

all: $(OBJ) $(LIBRARY)

-include $(DEP)

$(LIBRARY): $(OBJ)
	ar rcs $@ $^
	ranlib $@

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $^ -o $@ $(LIBS)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $^ $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP) $(LIBRARY)

.PHONY: clean
