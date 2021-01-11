CC=gcc
CFLAGS=-c -Wall -iquote/usr/local/lib/include

all: TDDA

TDDA: Lattice_Diffusion.o find_neighbor.o
	$(CC) Lattice_Diffusion.o find_neighbor.o -o Lattice_Diffusion -L/usr/local/gfortran/lib -lgfortran -llapack -lblas -lm

Lattice_Diffusion.o: Lattice_Diffusion.c
	$(CC) $(CFLAGS) Lattice_Diffusion.c

find_neighbor.o: find_neighbor.c
	$(CC) $(CFLAGS) find_neighbor.c

clean:
	rm *.o Lattice_Diffusion
