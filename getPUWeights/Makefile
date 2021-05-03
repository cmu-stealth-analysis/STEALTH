# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cpp src
vpath %.h include
vpath % bin

all: makePUWeights

makePUWeights: makePUWeights.o
	g++ -o bin/makePUWeights obj/makePUWeights.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

makePUWeights.o: PUWeights.cpp PUWeights.h
	g++ -c -g -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/makePUWeights.o src/PUWeights.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
