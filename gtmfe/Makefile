# Sage directory information (ew)
SAGE_ROOT = $(shell sage --root)
SAGE_ENV = $(SAGE_ROOT)/local/bin/sage-env
SAGE_INC = $(SAGE_ROOT)/local/include
SAGE_LIB = $(SAGE_ROOT)/local/lib
SAGE_LIB64 = $(SAGE_ROOT)/local/lib64

LD_PATH_SET = LD_LIBRARY_PATH="$(SAGE_LIB):$(SAGE_LIB64)"
SAGE_GCC = $(LD_PATH_SET) $(SAGE_ROOT)/local/bin/gcc
CC = $(SAGE_GCC)
CPP = $(SAGE_GCC)
CXX = $(SAGE_GCC)

# source files.
SRC = src/*.cc
SRC-WITH = $(SRC) gtmfe-with-openmp.cc
SRC-WITHOUT = $(SRC) gtmfe-without-openmp.cc

OBJ = *.o

OUT = _gtmfe.so

# include directories
INCLUDES += -Iinclude

# subdirectories of sage include/ to include
RINCLUDES += "python2.7" "gmp++"
# And then this nightmare
INCLUDES += -I$(SAGE_INC) $(addprefix -I,$(foreach dir,$(RINCLUDES),$(shell find $(SAGE_INC)/$(dir)  -type d -print)))

# C++ compiler flags (-g -O2 -Wall)
CXXFLAGS += -shared -fPIC -Wall -g -O2 #-nostdinc
CXXFLAGS-OPENMP = -fopenmp -DHAVE_OPENMP

# library paths
LIBS += -Wl,-rpath,$(SAGE_LIB) -Wl,-rpath,$(SAGE_LIB64)
LIBS += -L$(SAGE_LIB) -L$(SAGE_LIB64) -lgmp -lm
LIBS-OPENMP = -lgomp

gtmfe-with-openmp: $(SRC-WITH)
	$(LINK.cc) $(INCLUDES) $(CXXFLAGS-OPENMP) $(LIBS) $(LIBS-OPENMP) $^ -o $(OUT)

gtmfe-without-openmp: $(SRC-WITHOUT)
	$(LINK.cc) $(INCLUDES) $(LIBS) $^ -o $(OUT)

clean:
	-rm -vf $(OBJ)
