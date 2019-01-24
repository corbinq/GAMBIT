VERSION=0.2

CPPFLAGS= -Ofast -flto -pipe
CXXFLAGS= -std=c++11 -DNDEBUG
FFLAGS=

LDFLAGS= -lz -lm

cppsrc = $(wildcard src/*.cpp) src/tabix_util/tabix.cpp
csrc = $(wildcard src/*.c) src/tabix_util/index.c src/tabix_util/bgzf.c

objs = $(cppsrc:.cpp=.o) $(csrc:.c=.o) 

bin/GaMBIT: $(objs) src/libMvtnorm/libMvtnorm.a
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

src/libMvtnorm/libMvtnorm.a:
	cd src/libMvtnorm && $(MAKE)

.PHONY: clean

clean: 
	rm -f src/*.o src/*/*.o bin/*
	cd src/libMvtnorm && $(MAKE) clean

