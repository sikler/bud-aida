#!/bin/tcsh -f

if($#argv == 5) then
awk '{ \
 if($v>min && $v<max) \
 { i=int(($v-min)/(max-min)*bin); a[i]++; sum++ } \
} END { \
 for(i=0; i<bin; i++) print min+(i+0.5)/bin*(max-min),a[i]+0,sum \
}' v=$1 min=$2 max=$3 bin=$4 $5
endif

if($#argv == 9) then
awk '{ \
 if($x>xmin && $x<xmax && $y>ymin && $y<ymax) \
 { i=int(($x-xmin)/(xmax-xmin)*xbin); \
   j=int(($y-ymin)/(ymax-ymin)*ybin); \
   a[i,j]++; sum++ } \
} END { \
 for(i=0; i<xbin+1; i++) \
 { \
 for(j=0; j<ybin+1; j++) \
   print xmin+(i+0.5)/xbin*(xmax-xmin), \
         ymin+(j+0.5)/ybin*(ymax-ymin), a[i,j]+0,sum \
 print ""; \
 } \
}' x=$1 xmin=$2 xmax=$3 xbin=$4 y=$5 ymin=$6 ymax=$7 ybin=$8 $9
endif
