# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cpp src
vpath %.h include
vpath % bin

all: runEventSelection getParametersAndOptions

runEventSelection: runEventSelection.o
	g++ -g -o bin/runEventSelection obj/runEventSelection.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

runEventSelection.o: runEventSelection.cpp runEventSelection.h
	g++ -g -c -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/runEventSelection.o src/runEventSelection.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

getParametersAndOptions: getParametersAndOptions.o
	g++ -o bin/getParametersAndOptions obj/getParametersAndOptions.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

getParametersAndOptions.o: getParametersAndOptions.cpp getParametersAndOptions.h
	g++ -c -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/getParametersAndOptions.o src/getParametersAndOptions.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
