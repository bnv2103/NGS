args <- commandArgs(TRUE)

numReads = as.numeric(args[1])
outfile = args[2]
actgfile = args[3]
qualfile = args[4]
qscore_f = args[5]
qscore_r = args[6]

a = read.table(file=actgfile, header=F)
readLength = length(a[1,])-1

b = a[,2:(readLength+1)]
A = b[1,]
T = b[2,]
C = b[3,]
G = b[4,]
N = b[5,]

A = A / numReads
T = T / numReads
C = C / numReads
G = G / numReads
N = N / numReads

x = 1:readLength
pdf(file = outfile)

plot(0, xlim=c(1,readLength), ylim=c(0,1), xlab = "Cycle", ylab = "Frequency", main="ACGT Distribution",  type='n')
lines(x, A, lty=1, col="blue")
lines(x, C, lty=1, col="green")
lines(x, T, lty=1, col="red")
lines(x, G, lty=1, col="black")
lines(x, N, lty=1, col="yellow")
legend(80, 0.8, lty=c(1,1,1,1,1), col=c("blue", "green", "red", "black", "yellow" ),legend=c("A", "C", "T", "G", "N"))

qual1 = read.table(file=qualfile, header=F)
qual1 = as.numeric(qual1)
qual1 = qual1/numReads
plot(0, xlim=c(1,readLength), ylim=c(0,80), xlab = "Cycle", ylab = "Phred Quality", main="Mean Quality by Cycle", type='n')
lines(x, qual1, lty =1, col="darkgreen")


qscore1 = read.table(file=qscore_f, header=F)
maxX = nrow(qscore1)
qscore1 = as.numeric(t(qscore1))
temp= 1:maxX
maxY = max(qscore1)
plot(0, xlim=c(1,45), ylim=c(0,maxY),  xlab = "QScore", ylab = "Frequency", main="QScore Histogram Read1",  type='n')
lines(temp, qscore1, lty=1, col="blue")

qscore2 = read.table(file=qscore_r, header=F)
qscore2 = as.numeric(t(qscore2))
maxY = max(qscore2)
plot(0, xlim=c(1,45), ylim=c(0,maxY),  xlab = "QScore", ylab = "Frequency", main="QScore Histogram Read3",  type='n')
lines(temp, qscore2, lty=1, col="green")

dev.off()

