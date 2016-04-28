INCLUDE = -I/usr/include/gmm
LIBS    = lib/quad2D.o -lgsl -lblas
CPP     = clang++
FLAGS   = -O3


all: prog2 clean

#############################################

prog: prog.o
	$(CPP)  $(FLAGS) $(INCLUDE) prog.o $(LIBS) -o bin/prog

prog.o: prog.cxx
	$(CPP) -c $(FLAGS) $(INCLUDE) prog.cxx -o prog.o

#############################################

prog2: prog2.o
	$(CPP) $(FLAGS) $(INCLUDE) prog2.o $(LIBS) -o bin/prog2

prog2.o: prog2.cxx
	$(CPP) -c  $(FLAGS) $(INCLUDE) prog2.cxx -o prog2.o

#############################################


clean:
	rm *.o







