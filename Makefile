CC=g++
CFLAGS=-O3 -g -fopenmp -I/opt/local/var/macports/software/boost/1.46.1_0/opt/local/include -I. 
LDFLAGS=-L. -lgzstream -lz
EXE_FILES = hank
.PHONY: all clean

all: $(EXE_FILES)

clean:
	-rm $(EXE_FILES) *.o

hank: hank.cc Read.o hashtab.o edit.o libgzstream.a
	$(CC) $(CFLAGS) hank.cc Read.o hashtab.o edit.o -o hank $(LDFLAGS)

Read.o: Read.cc Read.h hashtab.o hank.o
	$(CC) $(CFLAGS) -c Read.cc

edit.o: edit.cc edit.h
	$(CC) $(CFLAGS) -c edit.cc

hashtab.o: hashtab.cc hashtab.h
	$(CC) $(CFLAGS) -c hashtab.cc

gzstream.o: gzstream.C gzstream.h
	$(CC) $(CFLAGS) -c -o gzstream.o gzstream.C

libgzstream.a: gzstream.o
	ar cr libgzstream.a gzstream.o
