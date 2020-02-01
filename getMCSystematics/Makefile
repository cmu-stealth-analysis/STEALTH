# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cpp src
vpath %.h include:../eventSelection/include
vpath % bin

all: makePUWeights getEventHistograms getMCUncertainties

makePUWeights: makePUWeights.o
	g++ -o bin/makePUWeights obj/makePUWeights.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

makePUWeights.o: makePUWeights.cpp makePUWeights.h
	g++ -c -g -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/makePUWeights.o src/makePUWeights.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

getEventHistograms: getEventHistograms.o
	g++ -o bin/getEventHistograms obj/getEventHistograms.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

getEventHistograms.o: getEventHistograms.cpp getEventHistograms.h STRegionsStruct.h shiftedObservablesStruct.h
	g++ -c -g -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/getEventHistograms.o src/getEventHistograms.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

getMCUncertainties: getMCUncertainties.o
	g++ -o bin/getMCUncertainties obj/getMCUncertainties.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs`

getMCUncertainties.o: getMCUncertainties.cpp getMCUncertainties.h STRegionsStruct.h shiftedObservablesStruct.h
	g++ -c -g -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/getMCUncertainties.o src/getMCUncertainties.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
