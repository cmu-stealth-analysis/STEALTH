# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cc src
vpath %.h include
vpath % bin

all: runFits

runFits: fitScripts.o
	g++ -o bin/runFits obj/fitScripts.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --ldflags --glibs` -lTreePlayer -lRooFit -lRooFitCore -lMinuit

fitScripts.o: fitScripts.cc fitScripts.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/fitScripts.o src/fitScripts.cc -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags`

clean:
	rm -rf bin/*
	rm -rf obj/*
