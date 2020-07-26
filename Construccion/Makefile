CPP=gcc

OBJECTS=basics.o bitrankw32int.o louds.o hash.o MemoryManager.o dacs.o dwtree.o 
BINS=builddwtree
#CPPFLAGS=-Wall #-g3 
CPPFLAGS=-Wall -O9 -g -DNDEBUG
DEST=.

%.o: %.c
	$(CPP) $(CPPFLAGS) -c $< -o $@

all: bin

bin: $(OBJECTS) $(BINS)
	
builddwtree:
	gcc $(CPPFLAGS) -o $(DEST)/builddwtree builddwtree.c $(OBJECTS) -lm

clean:
	rm -f $(OBJECTS) $(BINS)
	cd $(DEST); rm -f *.a *.o *.c~ *.h~ $(BINS)

