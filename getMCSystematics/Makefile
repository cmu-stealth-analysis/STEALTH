# Remember: order of arguments to g++ is important, the `root-config` must follow the rest of the arguments

JEC: JEC.o
	g++ -o bin/getMCSystematics_JEC obj/JEC.o -L$(TMCPPUTILS)/argumentParser/lib -ltmArgumentParser `root-config --cflags --glibs`

JEC.o: src/JEC.cpp include/JEC.h
	g++ -c -Wall -Wextra -Werror -pedantic -std=c++11 -fPIC -O3 -o obj/JEC.o src/JEC.cpp -I$(TMCPPUTILS)/argumentParser/include `root-config --cflags --glibs`

clean:
	rm -rf bin/*
	rm -rf obj/*
	rm -rf lib/*
