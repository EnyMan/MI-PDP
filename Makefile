CXX = g++ -pedantic -Wall -std=c++11

all: debug fast

debug: serial_simple.cpp
	$(CXX) -g serial_simple.cpp -o serial_dbg

fast: serial_simple.cpp
	$(CXX) -O3 serial_simple.cpp -o serial

clean:
	rm -f serial
	rm -f serial_dbg

.PHONY: all clean
