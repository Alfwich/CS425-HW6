CC = g++
CFLAGS = -std=c++11 
OPTS = -o3
TARGET = executable.x

#-Wall -Wextra

all: $(TARGET)

prof: CFLAGS += -pg
prof: TARGET = executable.prof.x
prof: $(TARGET)

$(TARGET): main.cpp tspsolver.o utils.o
	$(CC) $(CFLAGS) $(OPTS) tspsolver.o utils.o main.cpp -o $(TARGET)

tspsolver.o: tspsolver.cpp
	$(CC) $(CFLAGS) $(OPTS) tspsolver.cpp -c -o tspsolver.o

utils.o: utils.cpp
	$(CC) $(CFLAGS) $(OPTS) utils.cpp -c -o utils.o

clean:
	rm *.o *.x

.PHONY: prof
