CC=mpicxx
CFLAGS=-c -Wall
#CFLAGS=-c
LDFLAGS=
SOURCES=main.cpp cmc.cpp cmolecule.cpp rwpara.cpp ssbf.cpp 
#aglcalc.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mcmc

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -fr *.o
