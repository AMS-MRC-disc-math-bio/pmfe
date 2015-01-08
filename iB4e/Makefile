# source files
SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.d)
HDR = $(wildcard include/*.h)

LIBRARY = libiB4e.a

# include directories
INCLUDES += -Iinclude
INCLUDES += -I../gtmfe/include
INCLUDES += -I/usr/include/python2.7

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O2

# library paths
LIBS += -L../gtmfe -lgtmfe
LIBS += -lgmp
LIBS += -lCGAL
LIBS += -lboost_python
LIBS += -lpython2.7

all: $(OBJ) $(LIBRARY)

-include $(DEP)

$(LIBRARY): $(OBJ)
	ar rcs $@ $^
	ranlib $@

%.o: %.cc
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $*.cc -o $*.o $(LIBS)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $*.cc $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP) $(LIBRARY) $(LIBRARY)

.PHONY: clean
