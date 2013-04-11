f(x) = (abs(x) < 10*s ? abs(a) * exp(-0.5*(x/s)**2) : 0)
g(x) = (abs(x) < 10*s ? abs(a) * exp(-0.5*((x-m)/s)**2) : 0)

set fit errorvariables
set fit quiet

set print "../out/fitRes_pixel.dat"
call "fitRes_pixel_.gnu"  2
call "fitRes_pixel_.gnu"  3
call "fitRes_pixel_.gnu"  4
call "fitRes_pixel_.gnu"  5
call "fitRes_pixel_.gnu"  6
call "fitRes_pixel_.gnu"  7
call "fitRes_pixel_.gnu"  8
call "fitRes_pixel_.gnu"  9
call "fitRes_pixel_.gnu" 10

set print "../out/fitRes_strip.dat"
call "fitRes_strip_.gnu"  2
call "fitRes_strip_.gnu"  3
call "fitRes_strip_.gnu"  4
call "fitRes_strip_.gnu"  5
call "fitRes_strip_.gnu"  6
call "fitRes_strip_.gnu"  7
call "fitRes_strip_.gnu"  8
call "fitRes_strip_.gnu"  9
call "fitRes_strip_.gnu" 10

