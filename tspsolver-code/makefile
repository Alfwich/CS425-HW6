CC = g++
CFLAGS = -pthread -std=c++11 -fopenmp -Wl,--no-as-needed
OPTS = -o3
TARGET = executable.x

#-Wall -Wextra

all: $(TARGET)

prof: CFLAGS += -pg
prof: TARGET = executable.prof.x
prof: $(TARGET)
	
$(TARGET): main.cpp tspsolver.o
	$(CC) $(CFLAGS) $(OPTS) tspsolver.o main.cpp -o $(TARGET)

tspsolver.o: tspsolver.cpp
	$(CC) $(CFLAGS) $(OPTS) tspsolver.cpp -c -o tspsolver.o

clean:
	rm *.o *.x

.PHONY: prof
