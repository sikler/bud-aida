PROGS = ../bin/simulateEvents ../bin/plotEvents

all : $(PROGS)

###############################################################################
#TREEOBJS = TBunchCrossing.o TVertex.o TTrack.o TPixelHit.o TStripHit.o \
#           Coord.o Hit.o dict.o

#DEDXOBJS = CubicSpline.o LinearProgramming.o NewtonMinimizer.o \
#           ClusterGenerator.o ElossModel.o ModelBichsel.o TouchedChannels.o \
#           MostProbable.o FitStripCluster.o

FLAGS  = -Wall -O3 `root-config --cflags`
LFLAGS = -Wl,--no-as-needed
LIBS   = -lm -L/usr/local/lib -lCLHEP `root-config --libs` -lQuadp \
	-L../../DataFormats/lib      -ldict \
	-L../../siEnergyLoss/lib     -lsiEnergyLoss \
	-L../../minBiasVertexing/lib -lminBiasVertexing

###############################################################################
../bin/plotEvents : plotEvents.o \
			../../DataFormats/lib/libdict.a
	g++ $(FLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

../bin/simulateEvents : simulateEvents.o KalmanTracking.o \
			../../DataFormats/lib/libdict.a \
			../../siEnergyLoss/lib/libsiEnergyLoss.a \
			../../minBiasVertexing/lib/libminBiasVertexing.a
	g++ $(FLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

###############################################################################
#dict.so : dict.cc
#	g++ $(FLAGS) -o $@ $^ `root-config --ldflags` -shared
#
#dict.o : dict.cc
#	g++ $(FLAGS) -c -o $@ $<
#
#dict.cc : $(wildcard ../interface/T*.h ../interface/Coord.h ../interface/Hit.h) LinkDef.h
#	rm -f dict.C dict.h
#	rootcint -f $@ -c $^
#
###############################################################################
%.o : %.cc $(wildcard ../interface/*.h)
	g++ $(FLAGS) -c -o $@ $<

clean :
	rm -f $(PROGS) *.o
