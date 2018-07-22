# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cpp src
vpath %.h include
vpath % bin

all: getEventHistograms getJECUncertaintySystematics

getJECUncertaintySystematics: getJECUncertaintySystematics.o
	g++ -o bin/getJECUncertaintySystematics obj/getJECUncertaintySystematics.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs`

getJECUncertaintySystematics.o: getJECUncertaintySystematics.cpp getJECUncertaintySystematics.h
	g++ -c -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/getJECUncertaintySystematics.o src/getJECUncertaintySystematics.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

getEventHistograms: getEventHistograms.o
	g++ -o bin/getEventHistograms obj/getEventHistograms.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

getEventHistograms.o: getEventHistograms.cpp getEventHistograms.h
	g++ -c -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/getEventHistograms.o src/getEventHistograms.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
