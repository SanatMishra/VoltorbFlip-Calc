sources := $(wildcard src/*.cpp)
headers := $(wildcard hdr/*.h)
objects := $(sources:src/%.cpp=bin/%.o)

DEBUG := 0

cflags := -Ihdr -std=c++17

ifeq ($(DEBUG), 1)
	cflags += -g -Og -pg -fprofile-arcs -ftest-coverage
else
	cflags += -O3
endif

calc.exe : $(objects)
	g++ $(cflags) -o calc.exe $^

bin/main.o: src/main.cpp hdr/util.h hdr/VFCalc.h
	g++ $(cflags) -o $@ -c $<

bin/VFCalc.o: src/VFCalc.cpp hdr/util.h hdr/VFCalc.h
	g++ $(cflags) -o $@ -c $<

bin/getAllBoards.o: src/getAllBoards.cpp hdr/util.h hdr/VFCalc.h
	g++ $(cflags) -o $@ -c $<

bin/util.o: src/util.cpp hdr/util.h
	g++ $(cflags) -o $@ -c $<
