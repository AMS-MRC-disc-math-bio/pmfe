CXX = clang++

# source files
SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEP = $(SRC:.cc=.P)
HDR = $(wildcard src/*.h)

LIBRARY = libpmfe.a

BIN = pmfe-findmfe

# include directories
INCLUDES += -Iinclude
INCLUDES += -I/usr/include/python2.7

# C++ compiler flags
CXXFLAGS += -fPIC
CXXFLAGS += -Wall
CXXFLAGS += -g
CXXFLAGS += -O0

# library paths
LIBS += -lgmp
LIBS += -lCGAL
LIBS += -lm
LIBS += -lpython2.7
LIBS += -lboost_python
LIBS += -lboost_filesystem
LIBS += -lboost_program_options
LIBS += -lboost_system

# Optional OpenMP-related flags
HAS_OPENMP = $(shell $(CXX) -lgomp openmp-test.cc -o /dev/null 2>&1 > /dev/null ; echo $$?)
ifeq ($(HAS_OPENMP), 0)
CXXFLAGS += -fopenmp
CXXFLAGS += -DHAVE_OPENMP
LIBS += -lgomp
endif

all: $(OBJ) $(BIN)

-include $(DEP)

pmfe-findmfe: $(OBJ)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $^ -o $@ $(LIBS)

$(LIBRARY): $(OBJ)
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $^
	ranlib $@

%.o: %.cc
	$(CXX) -MD $(CXXFLAGS) $(INCLUDES) $(LIBS) $(LDFLAGS) -o $@ -c $<
	@cp $*.d $*.P; \
        sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
            -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
        rm -f $*.d

#	$(CXX) -MD $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -c $^ -o $@ $(LIBS)
#	$(CXX) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -MM $^ $(LIBS) | sed -e 's@^\(.*\)\.o:@src/\1.d src/\1.o:@' > $*.d

clean:
	-rm -vf $(EXEC) $(OBJ) $(DEP) $(LIBRARY) $(BIN)

install:

uninstall:

.PHONY: clean
