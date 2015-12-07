CC = g++
CFLAGS = -std=c++11
TARGET = executable.x

#-Wall -Wextra

all: $(TARGET)

debug: CFLAGS += -g
debug: TARGET = executable.debug.x
debug: $(TARGET)
	
$(TARGET): main.cpp tspsolver.o
	$(CC) $(CFLAGS) tspsolver.o main.cpp -o $(TARGET)

tspsolver.o: tspsolver.cpp
	$(CC) $(CFLAGS) tspsolver.cpp -c -o tspsolver.o

clean:
	rm *.o *.x

.PHONY: debug
