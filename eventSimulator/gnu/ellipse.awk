#!/usr/bin/awk -f

{
  dx = 0.1
  fwhm = 2 * sqrt(2*log(2))

  # Sum
  if($1 >= 2 && $1 <= 10)
  {
  if(i==0)
   printf("set object %d ellipse at %f,%f size 0.5,%f back fc ls 1 fs transparent solid 0.2; ",
           100*i+$1, $1-dx,$14,fwhm*$18)

  # Fitter
  if(i==1)
   printf("set object %d ellipse at %f,%f size 0.5,%f back fc ls 3 fs transparent pattern 5; ",
           100*i+$1, $1+dx,$16,fwhm*$20)
  }
}
