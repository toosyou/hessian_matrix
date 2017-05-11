CC=gcc-5
CXX=g++-5
INCLUDE=progressbar/include/
CXXFLAGS=-ltiff -fopenmp -lncurses -I$(INCLUDE) -Lprogressbar/ -lprogressbar -lgsl -lgslcblas

all: neuron_detection_in_tiff

neuron_detection_in_tiff:tomo_tiff.o main.o progressbar/libprogressbar.so
	$(CXX) tomo_tiff.o main.o $(CXXFLAGS) -o neuron_detection_in_tiff

tomo_tiff.o:tomo_tiff.cpp tomo_tiff.h progressbar/libprogressbar.so
	$(CXX) $(CXXFLAGS) -c tomo_tiff.cpp -o tomo_tiff.o

main.o:main.cpp tomo_tiff.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

progressbar/libprogressbar.so:progressbar/MakeFile
	cd progressbar && make && cd ..;

progressbar/MakeFile:
	git submodule update --init --recursive

clean:
	rm -f neuron_detection_in_tiff tomo_tiff.o main.o && cd progressbar && make clean;
