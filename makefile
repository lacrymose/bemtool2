INCLUDE = -I/usr/include/gmm
LIBS    = lib/quad2D.o -lgsl -lblas
CPP     = clang++

all: prog clean

#############################################

prog: prog.o
	$(CPP) -O3 $(INCLUDE) prog.o $(LIBS) -o bin/prog

prog.o: prog.cxx
	$(CPP) -c -O3 $(INCLUDE) prog.cxx -o prog.o

#############################################

prog2: prog2.o
	$(CPP) -O3 $(INCLUDE) prog2.o $(LIBS) -o bin/prog2

prog2.o: prog2.cxx
	$(CPP) -c -O3 $(INCLUDE) prog2.cxx -o prog2.o

#############################################


clean:
	rm *.o







