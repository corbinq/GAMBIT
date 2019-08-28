VERSION=0.3

CPPFLAGS= -Ofast -flto -pipe -I$(cdir)/src/htslib -I$(cdir)/src/cdflib -I$(cdir)/src/tabixpp -L$(cdir)/src/htslib -L$(cdir)/src/cdflib -L$(cdir)/src/tabixpp 
CXXFLAGS= -std=c++11 
FFLAGS=
LDFLAGS= -lz -lm

cppsrc = $(wildcard src/*.cpp) $(wildcard src/ROOT_Math/*.cpp) $(wildcard src/cdflib/*.cpp)  $(wildcard src/tabixpp/*.cpp)
csrc = $(wildcard src/*.c) src/htslib/bgzf.c src/htslib/kstring.c src/htslib/knetfile.c src/htslib/index.c

objs = $(cppsrc:.cpp=.o) $(csrc:.c=.o) 

cdir = ${CURDIR}

bin/GAMBIT: $(objs) src/libMvtnorm/libMvtnorm.a
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

src/libMvtnorm/libMvtnorm.a:
	cd src/libMvtnorm && $(MAKE)

.PHONY: clean

clean: 
	rm -f src/*.o src/*/*.o bin/*
	cd src/libMvtnorm && $(MAKE) clean

