CC = g++
CFLAGS = -std=c++11 -g
#-Wall -Wextra

all: executable.x 

executable.x: main.cpp tspsolver.o
	$(CC) $(CFLAGS) tspsolver.o main.cpp -o executable.x

tspsolver.o: tspsolver.cpp
	$(CC) $(CFLAGS) tspsolver.cpp -c -o tspsolver.o

clean:
	rm *.o *.x
