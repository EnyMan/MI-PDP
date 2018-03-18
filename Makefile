CXX = g++ -pedantic -Wall -std=c++11

FXX = g++ -pedantic -Wall -std=c++11 -fopenmp

all: debug fast

debug: serial_simple.cpp
	$(CXX) -g serial_simple.cpp -o serial_dbg

fast: serial_simple.cpp
	$(CXX) -O3 serial_simple.cpp -o serial

task: task.cpp
	$(FXX) -O3 task.cpp -o task

clean:
	rm -f serial
	rm -f serial_dbg

.PHONY: all clean
