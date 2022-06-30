a = 0
b = 10000
reset
#set terminal png
set terminal postscript enhanced eps  font 'Helvetica,16'
#set style data lines
set grid
set key left
#set xrange[0.04:0.5]
#set yrange[0.00:0.5]
v= sprintf("%g",a+1)
set output './image/blastE'.v.'.eps' #'roll1D'.v.'.eps'
set xlabel 'x'
set ylabel 'internal energy'
plot  './resu/blast1600.out' u 1:5 w p pt 1 title '1600 cells' , './resu/blast3200.out' u 1:5  w p pt 8 title '3200 cells', './resu/blast6400.out' u 1:5  w p pt 6  title '6400 cells'
pause -1
a=a+1
if(a<b) reread














plot  './resu/sod.out' u 1:2 w p pt 1 title 'num, 400 cells' , './resu/TESTsod.out' u 1:2  with lines lt -1 lw 4  title 'exact'  
pause -1

set ylabel 'velocity'
plot  './resu/sod.out' u 1:3 w p pt 1 title 'num, 400 cells' , './resu/TESTsod.out' u 1:3  with lines lt -1 lw 4  title 'exact'  
set output './image/sodu'.v.'.eps' #'roll1D'.v.'.eps'
pause -1
set ylabel 'pressure'
plot  './resu/sod.out' u 1:4 w p pt 1 title 'num, 400 cells' , './resu/TESTsod.out' u 1:4  with lines lt -1 lw 4  title 'exact'  
set output './image/sodp'.v.'.eps' #'roll1D'.v.'.eps'
pause -1
set ylabel 'energy'
plot  './resu/sod.out' u 1:5 w p pt 1 title 'num, 400 cells' , './resu/TESTsod.out' u 1:5  with lines lt -1 lw 4  title 'exact'  
pause -1
