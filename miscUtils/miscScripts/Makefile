# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cc src
vpath %.h include
vpath % bin

all: getSumMCWeights makeHistogramsLowST makeSTHistogramsUnifiedSelectionQCD optimizeNormBin

getSumMCWeights: getSumMCWeights.o
	g++ -o bin/getSumMCWeights obj/getSumMCWeights.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --ldflags --glibs`

getSumMCWeights.o: getSumMCWeights.cc getSumMCWeights.h
	g++ -g -c -std=c++14 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/getSumMCWeights.o src/getSumMCWeights.cc -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags`

makeHistogramsLowST: makeHistogramsLowST.o
	g++ -o bin/makeHistogramsLowST obj/makeHistogramsLowST.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --ldflags --glibs` -lTreePlayer

makeHistogramsLowST.o: makeHistogramsLowST.cc makeHistogramsLowST.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/makeHistogramsLowST.o src/makeHistogramsLowST.cc -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags`

makeSTHistogramsUnifiedSelectionQCD: makeSTHistogramsUnifiedSelectionQCD.o
	g++ -o bin/makeSTHistogramsUnifiedSelectionQCD obj/makeSTHistogramsUnifiedSelectionQCD.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --ldflags --glibs` -lTreePlayer

makeSTHistogramsUnifiedSelectionQCD.o: makeSTHistogramsUnifiedSelectionQCD.cc makeSTHistogramsUnifiedSelectionQCD.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/makeSTHistogramsUnifiedSelectionQCD.o src/makeSTHistogramsUnifiedSelectionQCD.cc -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags`

optimizeNormBin: optimizeNormBin.o
	g++ -o bin/optimizeNormBin obj/optimizeNormBin.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --ldflags --glibs` -lTreePlayer -lRooFit -lRooFitCore -lMinuit

optimizeNormBin.o: optimizeNormBin.cc optimizeNormBin.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/optimizeNormBin.o src/optimizeNormBin.cc -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags`

clean:
	rm -rf bin/*
	rm -rf obj/*
