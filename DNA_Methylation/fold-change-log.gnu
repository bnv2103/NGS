# Need to check the max/min options for this one

reset
n=100	#number of intervals
max=10.0	#max value
min=-10.0	#min value
width=(max-min)/n	#interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set term png	#output terminal and file
set output "fold-change-log.png"
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5	#fillstyle
set tics out nomirror
set xlabel "Log fold change between group A and B"
set ylabel "Frequency"
#count and plot
plot	"stat.txt" u (hist($6,width)):(1.0) smooth freq w boxes lc rgb"red" notitle

