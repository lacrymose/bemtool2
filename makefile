INCLUDE = -I/usr/include/gmm -I/usr/include/boost/math/special_functions/
LIBS    = lib/quad2D.o -lgsl -lblas
CPP     = clang++
FLAGS   = -g

all: clear prog3 clean

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

prog3: prog3.o
	$(CPP) $(FLAGS) $(INCLUDE) prog3.o $(LIBS) -o bin/prog3

prog3.o: prog3.cxx
	$(CPP) -c  $(FLAGS) $(INCLUDE) prog3.cxx -o prog3.o

#############################################

clean:
	rm *.o

clear:
	clear







