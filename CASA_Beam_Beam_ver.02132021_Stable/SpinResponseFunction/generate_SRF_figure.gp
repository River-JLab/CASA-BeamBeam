set term png
set output "Files_of_Output/SpinResponseFunction.png"



set grid
# set multiplot
set size 1.0,1.0
# set xrange [0:2256.614542]
# set yrange [-100:500]
set xlabel "z [m]"
set ylabel "|F|"
# set y2label "|F|"
# set y2range [-0.5:2.5]
# set y2tics 0.5 nomirror
# set y2tics format "%2.1f"
# set key box lw 0.5 width 1.9 height 1.15

# set lmargin 8
# set rmargin 9
# set tmargin 2.5
plot 'Files_of_Output/SpinFunctionOutput.casa' using 1:2 w l lc rgb 'red' lw 2.5 title 'Spin Responce'
