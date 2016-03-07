target=bin/fakems

all:$(target)

INC=-I ../mscreate/include -I /usr/local/include/casacore/
CXXFLAGS=-O3
LDFLAGS= -lcasa_casa -lcasa_ms -lcasa_measures -lcasa_tables -L ../mscreate/lib -lmscreate

bin/fakems:obj/fakems.o
	mkdir -p bin&&g++ -o $@ obj/fakems.o $(LDFLAGS) -g

obj/fakems.o:src/fakems.cc
	mkdir -p obj&&g++ -c $<  $(CXXFLAGS) $(INC) -g -o $@

clean:
	rm -f `find -iname *.o` `find -iname *~`

