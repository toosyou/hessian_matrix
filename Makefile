CC=gcc-7
CXX=g++-7
INCLUDE = progressbar/include/progressbar
CXXFLAGS += -fopenmp
CXXFLAGS += -lncurses -I$(INCLUDE) -Lprogressbar/ -lprogressbar -lgsl -lgslcblas

all: test

test: hessian_matrix.o main.o
	$(CXX) hessian_matrix.o main.cpp $(CXXFLAGS) -o test

hessian_matrix.o:hessian_matrix.cpp hessian_matrix.h progressbar/libprogressbar.so
	$(CXX) -c hessian_matrix.cpp $(CXXFLAGS) -o hessian_matrix.o

main.o: main.cpp hessian_matrix.h
	$(CXX) -c main.cpp $(CXXFLAGS) -o main.o

progressbar/libprogressbar.so:progressbar/MakeFile
	cd progressbar && make && cd ..;

progressbar/MakeFile:
	git submodule update --init --recursive

clean:
	rm -f test hessian_matrix.o main.o && cd progressbar && make clean;
