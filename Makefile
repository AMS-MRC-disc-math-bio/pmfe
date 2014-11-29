# source files.
SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.d)
HDR = $(wildcard include/*.h)

OBJ-GTMFE = $(filter-out src/iB4e_main.o, $(OBJ))
OBJ-PARAM = $(filter-out src/gtmfe_main.o, $(OBJ))
EXEC = gtmfe-param parametrizer

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

all: cgal $(OBJ) $(EXEC)

gtmfe-param: $(OBJ-GTMFE)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(OBJ-GTMFE) -o $@ $(LIBS)

parametrizer: $(OBJ-PARAM)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(OBJ-PARAM) -o $@ $(LIBS)

-include $(DEP)

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $*.cc -o $*.o $(LIBS)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $*.cc $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

CGAL/Makefile:
	cmake CGAL

cgal: CGAL/Makefile
	$(MAKE) -C CGAL

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP)
	$(MAKE) -C CGAL clean

.PHONY: clean
