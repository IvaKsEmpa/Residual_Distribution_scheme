reset
set style data lines    #p
set xlabel 'x' 
set ylabel 'ro'
plot  './resu/blast1600.out' u 1:2 w p pt 1 title '1600 cells' , './resu/blast3200.out' u 1:2  w p pt 2 title '3200 cells', './resu/blast6400.out' u 1:2  w p pt 3  title '6400 cells'  #, './resu/TEST123Exact.out' u 1:2 w l title 'exact'

pause(-1)
set xlabel 'x' 
set ylabel 'u'
plot  './resu/blast1600.out' u 1:3 w p pt 1 title '1600 cells' , './resu/blast3200.out' u 1:3  w p pt 2 title '3200 cells', './resu/blast6400.out' u 1:3  w p pt 3  title '6400 cells'  #, './resu/TEST123Exact.out' u 1:3 w l title 'exact'

pause(-1)
set xlabel 'x' 
set ylabel 'ein'
plot  './resu/blast1600.out' u 1:4 w p pt 1 title '1600 cells' , './resu/blast3200.out' u 1:4  w p pt 2 title '3200 cells', './resu/blast6400.out' u 1:4  w p pt 3  title '6400 cells'  #, './resu/TEST123Exact.out' u 1:4 w l title 'exact'


pause(-1)
set xlabel 'x' 
set ylabel 'pressure'
plot  './resu/blast1600.out' u 1:5 w p pt 1 title '1600 cells' , './resu/blast3200.out' u 1:5  w p pt 2 title '3200 cells', './resu/blast6400.out' u 1:5  w p pt 3  title '6400 cells'  #, './resu/TEST123Exact.out' u 1:5 w l title 'exact'




