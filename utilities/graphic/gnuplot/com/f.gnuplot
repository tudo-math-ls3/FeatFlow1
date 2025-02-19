set key 8.0,1.5
set terminal postscript landscape
set output "f.ps"
set xlabel "time"
set ylabel "flux - values"
set size 0.9,0.9
plot [0.0:10.1] [0.1:0.3] "f1_0" title 'f1' w l 1


!pageview -right -Ws 920 750 f.ps
