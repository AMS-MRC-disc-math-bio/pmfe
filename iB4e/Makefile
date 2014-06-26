CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cc BBpolytope.cc euclideanvector.cc faces.cc linalg.cc stack.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=iB4e

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS)
