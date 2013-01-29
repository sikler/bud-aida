all :
	cd DataFormats/src      ; make
	cd minBiasVertexing/src ; make
	cd siEnergyLoss/src     ; make
	cd eventSimulator/src   ; make
	cd eventSimulator/bin   ; ./simulateEvents ../data/minBias.dat
#	cd eventSimulator/bin   ; make

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
            -cvzf eventSimulator.tgz eventSimulator
	tar --exclude=\*svn\* --exclude=\*.o --exclude=\*.a --exclude=\*.m \
            --exclude=\*.dat  --exclude=\*.eps --exclude=\*.zip \
            -cvzf DataFormats.tgz DataFormats
	tar --exclude=\*svn\* --exclude=\*.o --exclude=\*.a --exclude=\*.m \
            --exclude=\*.dat  --exclude=\*.eps --exclude=\*.zip \
            -cvzf minBiasVertexing.tgz minBiasVertexing
	tar --exclude=\*svn\* --exclude=\*.o --exclude=\*.a --exclude=\*.m \
            --exclude=\*.dat  --exclude=\*.eps --exclude=\*.zip \
            -cvzf siEnergyLoss.tgz siEnergyLoss
	mv *.tgz tgz/
