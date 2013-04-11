#!/usr/bin/awk -f

{
  if(ARGIND == 1) # original
  {
    o[$1","$2","$3] = $4

                               l[$1","$2","$3] = "strip"
    if($1<5 || $1==7 || $1==8) l[$1","$2","$3] = "pixel"
  }

  if(ARGIND == 2) # estimated
  {
    v[$1","$2","$3] = $5
    e[$1","$2","$3] = $6
  }
} END {
  for(x in o) print 1/o[x], v[x]+0, e[x]+0, l[x]
}
