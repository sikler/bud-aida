OBJS = CubicSpline.o NewtonMinimizer.o \
       ClusterGenerator.o ElossModel.o ModelBichsel.o TouchedChannels.o \
       MostProbable.o FitStripCluster.o ClusterReco.o ElossEstimator.o

all : ../lib/libsiEnergyLoss.a

###############################################################################
CFLAGS = -Wall -O3 `root-config --cflags`
LIBS   = -lm `root-config --libs`

###############################################################################
../lib/libsiEnergyLoss.so : $(OBJS)
	g++ -fPIC -o $@ $^ $(LIBS) -shared

../lib/libsiEnergyLoss.a : $(OBJS)
	ar rs $@ $^
 
###############################################################################
%.o : %.cc $(wildcard ../../*/interface/*.h)
	g++ $(CFLAGS) -c -o $@ $<

clean :
	rm -f $(OBJS) ../lib/*.a ../lib/*.so

