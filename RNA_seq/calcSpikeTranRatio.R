args <- commandArgs(TRUE)

# file1 = args[1]
# file2 = args[2]

fileList = list.files(pattern="*.csv")
pdf(file = "outfile.pdf")
for (i in 1:length(fileList) ){
    for (j in 1:(length(fileList)) ){
    	if (i > j){

a = read.csv(fileList[i], head=T)
b = read.csv(fileList[j], head=T)

a$ratio = a$V1/sum(a$V1)
b$ratio = b$V1/sum(b$V1)

rownames(a) = a$V2
rownames(b) = b$V2


commID = intersect(a$V2, b$V2)

d = cbind(a[commID,], b[commID,])

# pdf(file = "outfile.pdf")

x = strsplit(fileList[i], "_")
y = strsplit(fileList[j], "_")

plot(d[,4], d[,8],  xlab=x[[1]][6], ylab=y[[1]][6],  main="Ratio", pch=19)
lines(lowess(d[,4], d[,8]), col="blue")
sum = summary(lm(d[,4] ~ d[,8]))
r2 = round(sum$r.squared, digits=3)
r2Test = substitute(R^2 == r2, list(r2=r2))
text(0.25, 0.25, r2Test, col="blue")
}
}
}
dev.off()