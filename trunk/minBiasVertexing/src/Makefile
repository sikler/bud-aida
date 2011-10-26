FFLAGS = -Wall -O3
CFLAGS = -Wall -O3 `root-config --cflags`

LIBS = -lm `root-config --libs`

all : ../bin/vertexFinder

../bin/vertexFinder : vertexFinder.o \
                      PairGroupMethod.o FastPairwiseNearestNeighbor.o \
                      NeighborJoining.o \
                      GaussianMixture.o KMeansMethod.o OptimalTree.o
	g++ $(CFLAGS) -o $@ $^ $(LIBS)

vertexFinder.o : vertexFinder.cc
	g++ $(CFLAGS) -c -o $@ $<

#DivisiveMethod.o : DivisiveMethod.cc ../interface/DivisiveMethod.h \
#                   $(wildcard CommonTools/Clustering1D/*/*) \
#                   $(wildcard DataFormats/GeometryCommonDetAlgo/*/*) \
#                   $(wildcard RecoPixelVertexing/PixelVertexFinding/*/*)
#	g++ $(CFLASS) -c -o $@ $<

PairGroupMethod.o : PairGroupMethod.cc ../interface/PairGroupMethod.h
	g++ $(CFLAGS) -c -o $@ $<

FastPairwiseNearestNeighbor.o : FastPairwiseNearestNeighbor.cc \
                   ../interface/FastPairwiseNearestNeighbor.h
	g++ $(CFLAGS) -c -o $@ $<

NeighborJoining.o : NeighborJoining.cc ../interface/NeighborJoining.h
	g++ $(CFLAGS) -c -o $@ $<

GaussianMixture.o : GaussianMixture.cc ../interface/GaussianMixture.h
	g++ $(CFLAGS) -c -o $@ $<

KMeansMethod.o : KMeansMethod.cc ../interface/KMeansMethod.h
	g++ $(CFLAGS) -c -o $@ $<

OptimalTree.o : OptimalTree.cc ../interface/OptimalTree.h
	g++ $(CFLAGS) -c -o $@ $<

#../bin/minimumBias : minimumBias.o pythia-6.4.20.o
#	gfortran $(FFLAGS) -o $@ $^
#
#minimumBias.o : minimumBias.f
#	gfortran $(FFLAGS) -c -o $@ $<
#
#pythia-6.4.20.o : pythia-6.4.20.f
#	gfortran $(FFLAGS) -c -o $@ $<

clean :
	rm -f *.o ../bin/vertexFinder ../bin/minimumBias
