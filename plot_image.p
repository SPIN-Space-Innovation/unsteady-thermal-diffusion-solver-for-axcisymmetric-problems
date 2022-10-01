plot 'results.txt' with image
set size ratio 0.125
set xlabel 'x [m]'
set ylabel 'r [m]'
set title 'Temperature Field'
set term png size 3000,2000
set output 'Temperature_field.png'
replot