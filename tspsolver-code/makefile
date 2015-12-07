CC = g++
CFLAGS = -std=c++11 -o3
TARGET = executable.x

#-Wall -Wextra

all: $(TARGET)

prof: CFLAGS += -pg
prof: TARGET = executable.prof.x
prof: $(TARGET)
	
$(TARGET): main.cpp tspsolver.o
	$(CC) $(CFLAGS) tspsolver.o main.cpp -o $(TARGET)

tspsolver.o: tspsolver.cpp
	$(CC) $(CFLAGS) tspsolver.cpp -c -o tspsolver.o

clean:
	rm *.o *.x

.PHONY: prof
