OBJS = TBunchCrossing.o TVertex.o TTrack.o \
	TChipId.o TLayer.o TSlimMeasurement.o \
	Point.o Hit.o Pixel.o Crossing.o \
	RandomGenerator.o \
	dict.o

all : ../lib/libdict.a dict.so dict.o

###############################################################################
CFLAGS  = -Wall -O3 `root-config --cflags`

###############################################################################
../lib/libdict.a : $(OBJS)
	ar rs $@ $^

dict.so : dict.cc
	g++ $(CFLAGS) -o $@ $^ `root-config --ldflags` -shared

dict.o : dict.cc
	g++ $(CFLAGS) -c -o $@ $<

dict.cc : $(wildcard ../interface/T*.h) \
		../interface/Point.h ../interface/Hit.h ../interface/Pixel.h \
		../interface/Crossing.h LinkDef.h
	rm -f dict.C dict.cc dict.h
	rootcint -v4 -f $@ -c $^

###############################################################################
%.o : %.cc $(wildcard ../interface/*.h)
	g++ $(CFLAGS) -c -o $@ $<

clean :
	rm -f ../lib/*.a *.o *.so dict.cc dict.h
