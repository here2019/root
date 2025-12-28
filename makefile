CXX=g++
CXXFLAGS=-O3 -march=native -DNDEBUG -fopenmp

# g++ 9: C++20 often via -std=c++2a
STD=-std=c++2a

all: wheel105_pi

wheel105_pi: wheel105_pi_omp.cpp
	$(CXX) $(CXXFLAGS) $(STD) $< -o $@

run:
	./wheel105_pi 100000000 16 4 5

clean:
	rm -f wheel105_pi
