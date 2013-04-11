OBJS = PairGroupMethod.o FastPairwiseNearestNeighbor.o \
       GaussianMixture.o KMeansMethod.o OptimalTree.o VertexFinder.o

all : ../lib/libminBiasVertexing.a

###############################################################################
CFLAGS = -Wall -O3 `root-config --cflags`
LIBS   = -lm `root-config --libs`

###############################################################################
../lib/libminBiasVertexing.so : $(OBJS)
	g++ -fPIC -o $@ $^ $(LIBS) -shared

../lib/libminBiasVertexing.a : $(OBJS)
	ar rs $@ $^

###############################################################################
%.o : %.cc $(wildcard ../../*/interface/*.h)
	g++ $(CFLAGS) -c -o $@ $<

clean :
	rm -f $(OBJS) ../lib/*.a ../lib/*.so
