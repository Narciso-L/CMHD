CPP=gcc

OBJECTS=basics.o bitrankw32int.o louds.o hash.o MemoryManager.o dacs.o
BINS=search
#CPPFLAGS=-Wall -g3 
CPPFLAGS=-Wall -O9 -g -DNDEBUG
DEST=.

%.o: %.c
	$(CPP) $(CPPFLAGS) -c $< -o $@

all: bin

bin: $(OBJECTS) $(BINS)
	
search:
	gcc $(CPPFLAGS) -o $(DEST)/search search.c $(OBJECTS) -lm 

clean:
	rm -f $(OBJECTS) $(BINS)
	cd $(DEST); rm -f *.a *.o *.c~ *.h~ $(BINS)

