target=bin/fakems bin/raw2ms bin/raw2ms_splited bin/raw2ms_splited_lp bin/raw2ms_splited_precal bin/test_date_str bin/test_uvw_with_uvmap bin/test_uvw bin/bin2uvmap bin/ms2uvmap bin/bfraw2ms

all:$(target)

INC=-I ../mscreate/include -I /usr/local/include/casacore/ -I ./include
CXXFLAGS=-O3 -std=c++11
LDFLAGS=-L ../mscreate/lib -lmscreate -lcasa_casa -lcasa_ms -lcasa_measures -lcasa_tables -lfio -lcfitsio

bin/fakems:obj/fakems.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/fakems.o:src/fakems.cpp
	mkdir -p obj&&$(CXX) -c $<  $(CXXFLAGS) $(INC) -g -o $@

bin/raw2ms:obj/raw2ms.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms.o:src/raw2ms.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/bfraw2ms:obj/bfraw2ms.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/bfraw2ms.o:src/bfraw2ms.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/raw2ms_splited:obj/raw2ms_splited.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms_splited.o:src/raw2ms_splited.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/raw2ms_splited_lp:obj/raw2ms_splited_lp.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms_splited_lp.o:src/raw2ms_splited_lp.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/raw2ms_splited_precal:obj/raw2ms_splited_precal.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms_splited_precal.o:src/raw2ms_splited_precal.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/test_uvw_with_uvmap:obj/test_uvw_with_uvmap.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/test_uvw_with_uvmap.o:src/test_uvw_with_uvmap.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/test_uvw:obj/test_uvw.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/test_uvw.o:src/test_uvw.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/bin2uvmap:obj/bin2uvmap.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g -lfio -lcfitsio

obj/bin2uvmap.o:src/bin2uvmap.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/ms2uvmap:obj/ms2uvmap.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g -lfio -lcfitsio

obj/ms2uvmap.o:src/ms2uvmap.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

bin/test_date_str:obj/test_date_str.o obj/date_time.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/test_date_str.o: src/test_date_str.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/date_time.o:src/date_time.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

clean:
	rm -f `find -iname *.o` `find -iname *~`

