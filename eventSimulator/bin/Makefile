all : ../view/bunchCrossing.zip \
      ../view/clusters.zip

../view/bunchCrossing.zip : ../out/bunchCrossing.m
	cd ../out ; zip -q ../view/bunchCrossing.zip bunchCrossing.m

../view/clusters.zip : ../out/clusters.m
	cd ../out ; zip -q ../view/clusters.zip clusters.m

../out/bunchCrossing.m : plotEvents ../out/bunchCrossings.root
	./$<

../out/bunchCrossings.root : simulateEvents ../data/minBias.dat
	./$<

# setenv LD_LIBRARY_PATH
# /usr/local/lib:/usr/local/lib/root:/home/sikler/aida/bud-aida/siEnergyLoss/lib

#%.zip : %.m
#	zip -q -m $@ $<
