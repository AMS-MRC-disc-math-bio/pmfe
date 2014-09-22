# source files.
SRC = src/*.cc
SRC-WITH = $(SRC) gtmfe-with-openmp.cc
SRC-WITHOUT = $(SRC) gtmfe-without-openmp.cc

OBJ = *.o

OUT = _gtmfe.so

# include directories
INCLUDES += -Iinclude -I/usr/include/python2.7

# C++ compiler flags (-g -O2 -Wall)
CXXFLAGS += -shared -fPIC

# library paths
LIBS += -L/usr/local/lib -lgmp -lm
LIBS-WITH = $(LIBS) -lgomp

gtmfe-with-openmp: $(SRC-WITH)
	$(LINK.cc) $(INCLUDES) $(LIBS-WITH) $^ -o $(OUT)

gtmfe-without-openmp: $(SRC-WITHOUT)
	$(LINK.cc) $(INCLUDES) $(LIBS) $^ -o $(OUT)

clean:
	-rm -vf $(OBJ)
