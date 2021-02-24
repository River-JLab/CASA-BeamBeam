set term wxt size 1300,192
set object 1 rectangle from screen 0,0 to screen 1300,192 fc rgb "black" behind

set xrange [-0.004:0.004]
set yrange [-0.001:0.001]

set xlabel "Z" textcolor rgb "white"
set ylabel "X" textcolor rgb "white"

set tics textcolor rgb "white"
set border lc rgb "white"
set grid lc rgb "white"


do for [i=0:114]{
  datafile=sprintf("%0d_2.casa",i) # CASA_data input
  plot datafile u 5:1 w p pt 7 lc 2 lw 1 ps 0.2 title '' 
}

pause mouse any
