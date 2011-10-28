args <- commandArgs(TRUE)
outPDF = args[1]
isMA = args[2]
fileList = list.files(pattern="*isoforms*")
pdf(file = outPDF)
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    # tempName = paste(fileList[i], ".csv", sep="")
    # write.csv(a, tempName)
    B = strsplit(fileList[i], "_")
    title = B[[1]][4]
    brNum = 500
    
    hist(log10(a$FPKM), br=brNum, xlim=range(c(-10, 10)),xlab="log10(FPKM)", main = paste("Histogram of" , title))
    a.s = a[order(-a$FPKM*a$length), ]
    c = a.s[1:100, ]
    plot(1:100, log(c$FPKM*c$length), ylab="Normalized # Reads: log(FPKM*length)", main=paste("Plot of Top 100 Isoforms", title))
    tempName = paste(fileList[i], ".csv", sep="")
    write.csv(a, tempName)

}
# print MA plot
if (isMA ==1){
for (i in 1:(length(fileList)) ){
    for (j in i:(length(fileList)) ){
    	if (i != j){
	t1 = strsplit(fileList[i], "_")
        t2 = strsplit(fileList[j], "_")
        title1 = t1[[1]][4]
        title2 = t2[[1]][4]
	title_sub1 = strsplit(title1, "-")
	title_sub2 = strsplit(title2, "-")
	if (title_sub1[[1]][2] == title_sub2[[1]][2]){

    	   a = read.table(fileList[i], header=T)
    	   b = read.table(fileList[j], header=T)
    	   M = log2(a$FPKM) - log2(b$FPKM)
    	   A = (1/2) * (log2(a$FPKM) + log2(b$FPKM))
	
	   mainText = paste("M-A plot:", title1, "~", title2,  sep=" ")
           plot(A,M, xlab="A = (1/2) * (log2(FPKM) + log2(FPKM))", ylab="M = log2(FPKM) - log2(FPKM)", main=mainText, ylim=c(-20,20), xlim=c(-10,20))
	  }
	}
    }
 }
}
dev.off()
fileList = list.files(pattern="*genes*")
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    tempName = paste(fileList[i], ".csv", sep="")
    
    write.csv(a, tempName)

}
