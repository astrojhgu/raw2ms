target=bin/fakems bin/raw2ms bin/raw2ms_splited bin/test_date_str

all:$(target)

INC=-I ../mscreate/include -I /usr/local/include/casacore/ -I ./include
CXXFLAGS=-O3 -std=c++11
LDFLAGS=-L ../mscreate/lib -lmscreate -lcasa_casa -lcasa_ms -lcasa_measures -lcasa_tables 

bin/fakems:obj/fakems.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/fakems.o:src/fakems.cpp
	mkdir -p obj&&$(CXX) -c $<  $(CXXFLAGS) $(INC) -g -o $@

bin/raw2ms:obj/raw2ms.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms.o:src/raw2ms.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/raw2ms_splited:obj/raw2ms_splited.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms_splited.o:src/raw2ms_splited.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/test_date_str:obj/test_date_str.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/test_date_str.o: src/test_date_str.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/date_time.o:src/date_time.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

clean:
	rm -f `find -iname *.o` `find -iname *~`

