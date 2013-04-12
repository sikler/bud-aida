all : process view

process : ../out/bunchCrossings_reco.root
view : ../view/bunchCrossing.zip ../view/clusters.zip

#PRE = valgrind --leak-check=full
PRE =

NEVENTS = 500

BASEURL = http://localhost:8000/view
#BASEURL = http://sikler.web.cern.ch/sikler/bud-aida

###############################################################################
../view/%.zip : ../out/%.m
	cd ../out ; zip -q $@ $*.m

../out/bunchCrossing.m : plotEvents \
				../out/bunchCrossings_reco.root
	$(PRE) ./$< -event 6 -baseUrl $(BASEURL) \
	 -i ../out/bunchCrossings_reco.root

###############################################################################
#../out/gains.dat : calibrateGain \
#				../out/hits.bin
#	$(PRE) ./$< -i ../out/hits.bin -o $@
#
###############################################################################
#../out/hits.bin :		../out/bunchCrossings_reco.root
#	touch $@

../out/bunchCrossings_reco.root : reconstructEvents calibrateGain
#       				../out/bunchCrossings_simu.root 
	$(PRE) ./reconstructEvents -collectHits asPion \
	 -i ../out/bunchCrossings_simu.root -o $@
	#
	$(PRE) ./calibrateGain \
	 -i ../out/hits.bin -o ../out/gains_asPion.dat
	#
	$(PRE) ./reconstructEvents -collectHits ifIdentified \
	 -i ../out/bunchCrossings_simu.root -o $@
	#
	$(PRE) ./calibrateGain \
	 -i ../out/hits.bin -o ../out/gains_ifIdentified.dat \
	         -previousGain ../out/gains_asPion.dat
	#
	$(PRE) ./reconstructEvents -collectHits notNeeded \
	 -i ../out/bunchCrossings_simu.root -o $@

###############################################################################
../out/bunchCrossings_simu.root : simulateEvents \
				../data/pp_minBias.dat.gz ../data/cms.dat
	$(PRE) ./$< -numberOfBunchCrossings $(NEVENTS) \
	 -numberOfCollisionsPoissonMean 5.0 \
	 -interactionRegionSigmaZ       5.0 \
	 -randomGainHalfWidth           0.5 \
	 -i ../data/pp_minBias.dat.gz \
	 -o $@
