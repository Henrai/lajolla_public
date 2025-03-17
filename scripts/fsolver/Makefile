CXXFLAGS=-g -Wall -DNDEBUG -fopenmp -ffast-math  -I/usr/local/opt/freeglut/include -L/usr/local/opt/freeglut/lib -I/usr/local/opt/libpng/include -L/usr/local/opt/libpng/lib -lpng -I/usr/local/opt/boost/include
CXX=/usr/local/Cellar/llvm/19.1.7_1/bin/clang++
 
all: fsolver
clean:
	rm -f *.o fsolver

display.o: display.cpp fsolver.h
fsolver.o: fsolver.cpp fsolver.h

fsolver: display.o fsolver.o
	$(CXX) -fopenmp display.o fsolver.o -o fsolver 
