# Need to check the max/min options for this one

reset
n=100	#number of intervals
max=10.0	#max value
min=0.0	#min value
width=(max-min)/n	#interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set term png	#output terminal and file
set output "stat-rpkm-log.png"
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
#set boxwidth width*0.9
set style fill solid 0.5	#fillstyle
set tics out nomirror
set xlabel "Log2 of RPKM+1"
set ylabel "Frequency"
#count and plot
plot	"stat-Sample_RK1.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 1", \
	"stat-Sample_RK2.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 2", \
	"stat-Sample_RK3.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 3", \
	"stat-Sample_RK4.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 4", \
	"stat-Sample_RK5.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 5", \
	"stat-Sample_RK6.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 6", \
	"stat-Sample_RK7.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 7", \
	"stat-Sample_RK8.txt" u (hist($4,width)):(1.0) smooth freq w lines title"Sample 8"

