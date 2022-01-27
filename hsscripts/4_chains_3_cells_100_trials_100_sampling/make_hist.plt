set term pngcairo font "Arial,14" 
set colorsequence classic 
# 
data = "nw_hist.dat" 
set output "histgram.png "
#
set size square
# set xrange [0:1.0]
#set yrange [0:100]
#
set xlabel "Arg. Con."
set ylabel "Freq."

set style fill solid 0.5
set boxwidth 0.012859999999999983
#
plot data w boxes noti