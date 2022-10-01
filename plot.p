
plot 'test.dat' with line title 'numerically','theoritical_results.dat' with line title 'analytically'
set grid
set xlabel 'time [sec]'
set ylabel 'Temperature [C]'
set title 'Temperature at point x = L/2'
set term png size 1500,1000
set output 'TransientConduction.png'
replot
