INCLUDE = -I/usr/include/gmm
LIBS    = src/quad2D.o -lgsl -lblas
CPP     = clang++

all: prog clean

#############################################

prog: prog.o
	$(CPP) -O3 $(INCLUDE) prog.o $(LIBS) -o prog

prog.o: prog.cxx
	$(CPP) -c -O3 $(INCLUDE) prog.cxx -o prog.o

#############################################


clean:
	rm *.o







