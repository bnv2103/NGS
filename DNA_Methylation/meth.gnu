# Need to check the max/min options for this one

reset
n=100	#number of intervals
max=10.0	#max value
min=0.0	#min value
width=(max-min)/n	#interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set term png	#output terminal and file
set output "meth.png"
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
#set style fill solid 0.5	#fillstyle
set tics out nomirror
set xlabel "log(2)+1: Methylation rate[A/B]"
set ylabel "Frequency"
#count and plot
plot	"meth.txt" u (hist($2,width)):(1.0) smooth freq w lines lc rgb"blue" title"Methylation rate A", \
	"meth.txt" u (hist($3,width)):(1.0) smooth freq w lines lc rgb"green" title"Methylation rate B"

