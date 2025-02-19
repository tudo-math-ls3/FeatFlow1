set key 8.0,1.5
set terminal postscript landscape
set output "uall.ps"
set xlabel "time"
set ylabel "u-values"
set size 0.9,0.9
plot [0.0:10.1] [-1.2:2.2] "u1_0" title 'u1' w l 1 , "u2_0" title 'u2' w l 2,"u3_0" title 'u3' w l 0,"u4_0" title 'u4' w l 3


!pageview -right -Ws 920 750 all.ps
