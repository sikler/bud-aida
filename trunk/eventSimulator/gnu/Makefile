all : ../eps/energyLoss.eps ../eps/gain.eps \
	../eps/residual_x.eps \
	../eps/resolution_x_pixels.eps \
	.charge .chi \
	.fitPlot .vertex 

../eps/energyLoss.eps : energyLoss.gnu epsilon.gnu myStyleEps.gnu \
		../../siEnergyLoss/data/density_Sternheimer84.par \
		../out/estimates.dat
	gnuplot $<

../eps/gain.eps : gain.gnu myStyleEps.gnu compare.awk \
		../out/gains.dat
	gnuplot $<

../eps/residual_x.eps : result.gnu myStyleEps.gnu \
		../out/result_pixels.dat
	gnuplot $<

../eps/resolution_x_pixels.eps : resolution.gnu myStyleEps.gnu  \
		../out/result_pixels.dat
	gnuplot $<

.charge : charge.gnu myStyleEps.gnu \
          ../out/result_pixels.dat \
          ../out/result_strips.dat
	gnuplot $<
	touch $@

.fitPlot : fitPlot.gnu myStyleEps.gnu .fitRes ellipse.awk
	gnuplot $<
	touch $@

.fitRes : fitRes.gnu \
          fitRes_pixel_.gnu myStyleEps.gnu \
          ../out/result_pixels.dat
	gnuplot $<
	touch $@

.mass : mass.gnu mass_.gnu myStyleEps.gnu \
          $(wildcard ../out/*_*.spec)
	gnuplot $<
	touch $@

.chi : chi.gnu myStyleEps.gnu \
          $(wildcard ../out/chi.dat)
	gnuplot $<
	touch $@

.vertex : vertex.gnu ../out/result.dat myStyleEps.gnu \
        fitMicro.gnu fitMacro.gnu \
        ../out/multi.dat ../out/macro.dat
	gnuplot $<
	touch $@
