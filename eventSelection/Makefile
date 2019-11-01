# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cpp src
vpath %.h include
vpath % bin

all: runEventSelection mergeStatisticsHistograms mergeEventSelections

runEventSelection: runEventSelection.o
	g++ -o bin/runEventSelection obj/runEventSelection.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

runEventSelection.o: runEventSelection.cpp runEventSelection.h constants.h miscDataStructures.h parameters.h collections.h shiftedObservablesStruct.h selectionCriteria.h objectProperties.h statisticsHistograms.h MCRegions.h selectionRegionUtils.h hltEmulation.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/runEventSelection.o src/runEventSelection.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

mergeStatisticsHistograms: mergeStatisticsHistograms.o
	g++ -o bin/mergeStatisticsHistograms obj/mergeStatisticsHistograms.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

mergeStatisticsHistograms.o: mergeStatisticsHistograms.cpp mergeStatisticsHistograms.h constants.h miscDataStructures.h parameters.h collections.h shiftedObservablesStruct.h selectionCriteria.h objectProperties.h statisticsHistograms.h MCRegions.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/mergeStatisticsHistograms.o src/mergeStatisticsHistograms.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

mergeEventSelections: mergeEventSelections.o
	g++ -o bin/mergeEventSelections obj/mergeEventSelections.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

mergeEventSelections.o: mergeEventSelections.cpp mergeEventSelections.h
	g++ -g -c -std=c++11 -Wall -Wextra -Werror -pedantic-errors -fPIC -O3 -o obj/mergeEventSelections.o src/mergeEventSelections.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
