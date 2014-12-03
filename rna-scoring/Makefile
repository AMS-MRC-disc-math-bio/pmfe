# source files
SRC = $(wildcard *.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.d)
HDR = $(wildcard *.h)

LIBRARY = librnascoring.so

# include directories
INCLUDES += -Iinclude

# C++ compiler flags
CXXFLAGS += -fPIC
#CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O2

# Additional flags for building library
LIBFLAGS += -shared

# library paths
LIBS += -lgmp


all: $(OBJ) $(LIBRARY)

-include $(DEP)

%.so: $(OBJ)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) $(LIBFLAGS) -Wl,-soname,$@ $^ -o $@

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $*.cc -o $*.o $(LIBS)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $*.cc $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP)

.PHONY: clean
