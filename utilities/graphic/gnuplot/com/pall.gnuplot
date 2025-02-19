set key 8.0,1.5
set terminal postscript landscape
set output "pall.ps"
set xlabel "time"
set ylabel "p-values"
set size 0.9,0.9
plot [0.0:10.1] [-1.2:2.2] "p1_0" title 'p1' w l 1 , "p2_0" title 'p2' w d ,"p3_0" title 'p3' w l 0 , "p5_0" title 'p5' w l 2 ,"p7_0" title 'p7' w l 3


!pageview -right -Ws 920 750 pall.ps
