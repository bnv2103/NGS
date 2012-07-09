args <- commandArgs(TRUE)

numReads = as.numeric(args[1])
outfile = args[2]
a = read.table(file="ACTG.txt", header=F)

b = a[,2:102]
A = b[1,]
T = b[2,]
C = b[3,]
G = b[4,]

A = A / numReads
T = T / numReads
C = C / numReads
G = G / numReads


x = 1:101
pdf(file = outfile)

plot(0, xlim=c(1,101), ylim=c(0,1), xlab = "Cycle", ylab = "Frequency", main="ACGT Distribution",  type='n')
lines(x, A, lty=1, col="blue")
lines(x, C, lty=1, col="green")
lines(x, T, lty=1, col="red")
lines(x, G, lty=1, col="black")
legend(80, 0.8, lty=c(1,1,1,1), col=c("blue", "green", "red", "black"),legend=c("A", "C", "T", "G"))

qual1 = read.table(file="Qual.txt", header=F)
qual1 = as.numeric(qual1)
qual1 = qual1/numReads
plot(0, xlim=c(1,101), ylim=c(20,80), xlab = "Cycle", ylab = "Phred Quality", main="Per Cycle Read Quality", type='n')
lines(x, qual1, lty =1, col="darkgreen")
dev.off()