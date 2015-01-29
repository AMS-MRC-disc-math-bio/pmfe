# source files
SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.d)
HDR = $(wildcard include/*.h)

LIBRARY = librnascoring.a

# include directories
INCLUDES += -Iinclude

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O2

# library paths
LIBS += -lgmp

all: $(OBJ) $(LIBRARY)

-include $(DEP)

$(LIBRARY): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^
	ranlib $@

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $*.cc -o $*.o $(LIBS)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $*.cc $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP) $(LIBRARY)

install:

uninstall:

.PHONY: clean
