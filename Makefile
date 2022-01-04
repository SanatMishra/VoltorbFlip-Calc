sources := $(wildcard src/*.cpp)
headers := $(wildcard hdr/*.h)
objects := $(sources: %.cpp = bin/%.cpp)

DEBUG := 0

ifeq ($(DEBUG), 1)
	cflags += -g -Og -pg -fprofile-arcs -ftest-coverage
else
	cflags += -O4
endif

a.exe : $(objects)
	g++ $(cflags) a.exe $(objects)

src/calc.cpp: hdr/calc.h
	g++ -o calc.o -c calc.h calc.cpp

%.o : %.cpp
	g++ $(cflags)
