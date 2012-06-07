args <- commandArgs(TRUE)

fileIn = args[1]
fileOut = args[2]

a = read.table(fileIn, header=F)
a.s = a[order(-a$V1),]

fileOut = paste(fileOut, ".csv", sep="")
write.csv(a.s, fileOut)