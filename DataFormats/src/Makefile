OBJS = TBunchCrossing.o TVertex.o TTrack.o TPixelHit.o TStripHit.o \
	Coord.o Hit.o Pixel.o PixelUnit.o Crossing.o dict.o

FLAGS  = -Wall -O3 `root-config --cflags`

###############################################################################
all : ../lib/libdict.a dict.so dict.o

###############################################################################
../lib/libdict.a : $(OBJS)
	ar rs $@ $^

dict.so : dict.cc
	g++ $(FLAGS) -o $@ $^ `root-config --ldflags` -shared

dict.o : dict.cc
	g++ $(FLAGS) -c -o $@ $<

dict.cc : $(wildcard ../interface/T*.h) ../interface/Coord.h \
		../interface/Hit.h ../interface/Pixel.h \
		../interface/PixelUnit.h ../interface/Crossing.h LinkDef.h
	rm -f dict.C dict.h
	rootcint -f $@ -c $^

###############################################################################
%.o : %.cc $(wildcard ../interface/*.h)
	g++ $(FLAGS) -c -o $@ $<

clean :
	rm -f ../lib/*.a *.o *.so dict.cc dict.h
