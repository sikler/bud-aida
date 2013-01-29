CFLAGS = -Wall -O3 `root-config --cflags`

LIBS = -lm -L/usr/local/lib -lCLHEP `root-config --libs` -lQuadp

PROG = ../bin/clusterFit ../bin/generate ../bin/deposit \
       ../bin/fisher     ../bin/realLife ../bin/demo

all : ../lib/libsiEnergyLoss.a

OBJS = CubicSpline.o LinearProgramming.o NewtonMinimizer.o \
       ClusterGenerator.o ElossModel.o ModelBichsel.o TouchedChannels.o \
       MostProbable.o FitStripCluster.o RecoClusters.o ElossEstimator.o

../bin/demo : demo.o $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LIBS)

../bin/realLife : realLife.o $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LIBS)

../bin/clusterFit : clusterFit.o $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LIBS)

../bin/generate : generate.o $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LIBS)

../bin/deposit : deposit.o $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LIBS)

../bin/fisher : fisher.o $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LIBS)

#../bin/mostProbable : mostProbable.o $(OBJS)
#	g++ $(CFLAGS) -o $@ $^ $(LIBS)
#
#../bin/distribution : distribution.o
#	g++ $(CFLAGS) -o $@ $^ $(LIBS)

../lib/libsiEnergyLoss.so : $(OBJS)
	g++ -fPIC -o $@ $^ $(LIBS) -shared

../lib/libsiEnergyLoss.a : $(OBJS)
	ar rs $@ $^
 
%.o : %.cc $(wildcard ../interface/*.h)
	g++ $(CFLAGS) -c -o $@ $<

clean :
	rm -f *.o $(PROG)

