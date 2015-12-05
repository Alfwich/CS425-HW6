CC = g++
CFLAGS = -std=c++11 -Wall -Wextra

all: executable.x 

executable.x: main.cpp tspsolver.o
	$(CC) $(CFLAGS) tspsolver.o main.cpp -o executable.x

tspsolver.o: tspsolver.cpp
	$(CC) $(CFLAGS) tspsolver.cpp -o tspsolver.o
