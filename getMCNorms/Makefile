# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

vpath %.a lib
vpath %.so lib
vpath %.o obj
vpath %.cpp src
vpath %.h include
vpath % bin

all: common.o singlephoton

common.o: common.cpp common.h
	g++ -c -g -Wall -Wextra -Werror -pedantic -std=c++14 -fPIC -O3 -o obj/common.o src/common.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

singlephoton: singlephoton.o common.o
	g++ -o bin/singlephoton obj/singlephoton.o obj/common.o -L$(TMCPPUTILS)/generalUtils/lib -ltmGeneralUtils -L$(TMCPPUTILS)/ROOTUtils/lib -ltmROOTUtils `root-config --cflags --glibs` -lTreePlayer

singlephoton.o: singlephoton.cpp singlephoton.h common.h
	g++ -c -g -Wall -Wextra -Werror -pedantic -std=c++14 -fPIC -O3 -o obj/singlephoton.o src/singlephoton.cpp -I$(TMCPPUTILS)/generalUtils/include -I$(TMCPPUTILS)/ROOTUtils/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
