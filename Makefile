CXX = g++ -pedantic -Wall -std=c++11

MXX = mpicxx -pedantic -Wall -std=c++11

all: debug fast parallel mpi

debug: serial_simple.cpp
	$(CXX) -g serial_simple.cpp -o serial_dbg

fast: serial_simple.cpp
	$(CXX) -O3 serial_simple.cpp -o serial

parallel: serial_simple.cpp
	$(CXX) -fopenmp -O3 serial_simple.cpp -o parallel

mpi: mpi.cpp
	$(MXX) -fopenmp -g mpi.cpp -o mpi


clean:
	rm -f serial
	rm -f serial_dbg
	rm -f parallel
	rm -f mpi

.PHONY: all clean
