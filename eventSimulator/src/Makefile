PROGS = ../bin/simulateEvents \
	../bin/reconstructEvents \
	../bin/calibrateGain \
	../bin/plotEvents

all : $(PROGS)

###############################################################################
CFLAGS = -Wall -O3 `root-config --cflags`
LFLAGS = -Wl,--no-as-needed
LIBS   = -lm `root-config --libs` \
	-L../../siEnergyLoss/lib     -lsiEnergyLoss \
	-L../../minBiasVertexing/lib -lminBiasVertexing \
	-L../../DataFormats/lib      -ldict

###############################################################################
../bin/plotEvents : plotEvents.o \
			../../DataFormats/lib/libdict.a
	g++ $(CFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

../bin/calibrateGain : calibrateGain.o \
			../../DataFormats/lib/libdict.a \
			../../siEnergyLoss/lib/libsiEnergyLoss.a
	g++ $(CFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

../bin/reconstructEvents : reconstructEvents.o KalmanTracking.o \
			../../DataFormats/lib/libdict.a \
			../../siEnergyLoss/lib/libsiEnergyLoss.a \
			../../minBiasVertexing/lib/libminBiasVertexing.a
	g++ -Wno-unused-variable $(CFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

../bin/simulateEvents : simulateEvents.o KalmanTracking.o gzstream.o \
			../../DataFormats/lib/libdict.a \
			../../siEnergyLoss/lib/libsiEnergyLoss.a
	g++ $(CFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

###############################################################################
%.o : %.cc $(wildcard ../../*/interface/*.h)
	g++ $(CFLAGS) -c -o $@ $<

clean :
	rm -f $(PROGS) *.o ../out/hits.bin ../out/*.root
