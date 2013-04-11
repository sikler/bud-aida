all :
	cd DataFormats/src      ; make
	cd minBiasVertexing/src ; make
	cd siEnergyLoss/src     ; make
	cd eventSimulator       ; make

force :
	touch eventSimulator/data/cms.dat
	make all

plots :
	cd eventSimulator/gnu ; make
	cd eventSimulator/eps ; make
	cd eventSimulator/tex ; make

tex :
	cd eventSimulator/eps ; make
	cd eventSimulator/tex ; make

src :
	cd DataFormats/src      ; make
	cd minBiasVertexing/src ; make
	cd siEnergyLoss/src     ; make
	cd eventSimulator/src   ; make

clean ::
	cd DataFormats/src      ; make clean
	cd minBiasVertexing/src ; make clean
	cd siEnergyLoss/src     ; make clean
	cd eventSimulator/src   ; make clean

tar ::
	tar --exclude=\*svn\* --exclude=\*.o --exclude=\*.a --exclude=\*.m \
            --exclude=\*.dat  --exclude=\*.eps --exclude=\*.zip \
            --exclude=\*.root \
            -cvzf eventSimulator.tgz eventSimulator
	tar --exclude=\*svn\* --exclude=\*.o --exclude=\*.a --exclude=\*.m \
            --exclude=\*.dat  --exclude=\*.eps --exclude=\*.zip \
            --exclude=\*.root \
            -cvzf DataFormats.tgz DataFormats
	tar --exclude=\*svn\* --exclude=\*.o --exclude=\*.a --exclude=\*.m \
            --exclude=\*.dat  --exclude=\*.eps --exclude=\*.zip \
            --exclude=\*.root \
            -cvzf minBiasVertexing.tgz minBiasVertexing
	tar --exclude=\*svn\* --exclude=\*.o --exclude=\*.a --exclude=\*.m \
            --exclude=\*.dat  --exclude=\*.eps --exclude=\*.zip \
            --exclude=\*.root \
            -cvzf siEnergyLoss.tgz siEnergyLoss
	mv *.tgz tgz/

cc ::
	find . -name \*.cc -print

h ::
	find . -name \*.h -print
