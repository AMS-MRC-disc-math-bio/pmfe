# source files.
SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.d)
HDR = $(wildcard include/*.h)
EXEC = parametrizer

# include directories
INCLUDES += -Iinclude
INCLUDES += -ICGAL/include

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
#CXXFLAGS += -O2
CXXFLAGS += -O0
CXXFLAGS += -fopenmp
CXXFLAGS += -DHAVE_OPENMP

# library paths
LIBS += -LCGAL/lib -lCGAL
LIBS += -lgmp
LIBS += -lm
LIBS += -lgomp
LIBS += -lboost_system
LIBS += -lboost_filesystem

all: $(OBJ) $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(OBJ) -o $@ $(LIBS)

-include $(DEP)

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $*.cc -o $*.o $(LIBS)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $*.cc $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP)

.PHONY: clean
