set key 8.0,1.5
set terminal postscript landscape
set output "x.ps"
set xlabel "time"
set ylabel "value"
set size 0.9,0.9
plot [0.0:10.1] [1.2:1.6] "p1_0" title 'p1' w l 1


!pageview -right -Ws 920 750 x.ps
