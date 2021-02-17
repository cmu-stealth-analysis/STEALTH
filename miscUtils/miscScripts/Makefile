# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cc src
vpath %.h include
vpath % bin

all: makeSTHistogramsUnifiedSelectionQCD

makeSTHistogramsUnifiedSelectionQCD: makeSTHistogramsUnifiedSelectionQCD.o
	g++ -o bin/makeSTHistogramsUnifiedSelectionQCD obj/makeSTHistogramsUnifiedSelectionQCD.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

makeSTHistogramsUnifiedSelectionQCD.o: makeSTHistogramsUnifiedSelectionQCD.cc makeSTHistogramsUnifiedSelectionQCD.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/makeSTHistogramsUnifiedSelectionQCD.o src/makeSTHistogramsUnifiedSelectionQCD.cc -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
