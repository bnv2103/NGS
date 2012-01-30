args <- commandArgs(TRUE)

# 1. create csv files for isoforms and genes
# 2. print histogram for each sample
# 3. print MA plot for replicates (optional)

outPDF = args[1]
isMA = args[2]
fileList = list.files(pattern="*isoforms")
pdf(file = outPDF)
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    # tempName = paste(fileList[i], ".csv", sep="")
    # write.csv(a, tempName)
    B = strsplit(fileList[i], "_")
    # get the smaple name
    title = B[[1]][6]
    brNum = 500
    
    hist(log10(a$FPKM), br=brNum, xlim=range(c(-10, 10)),xlab="log10(FPKM)", main = paste("Histogram of" , title))
    a.s = a[order(-a$FPKM*a$length), ]
    c = a.s[1:100, ]
    plot(1:100, log(c$FPKM*c$length), ylab="Normalized # Reads: log(FPKM*length)", main=paste("Plot of Top 100 Isoforms", title))
    
    tempName = paste(B[[1]][5], B[[1]][6], sep="_")
    tempName = paste(tempName, B[[1]][7], sep="_")
    tempName = paste(tempName, B[[1]][8], sep="_")
    
    #tempName = paste(fileList[i], ".csv", sep="")
    tempName = paste(tempName, ".csv", sep="")
    write.csv(a, tempName)

}
# print MA plot
if (isMA ==1){
for (i in 1:(length(fileList)) ){
    for (j in i:(length(fileList)) ){
    	if (i != j){
	t1 = strsplit(fileList[i], "_")
        t2 = strsplit(fileList[j], "_")
        title1 = t1[[1]][6]
        title2 = t2[[1]][6]
	title_sub1 = strsplit(title1, "-")
	title_sub2 = strsplit(title2, "-")
	# if (title_sub1[[1]][2] == title_sub2[[1]][2]){
    	   a = read.table(fileList[i], header=T)
    	   b = read.table(fileList[j], header=T)
    	   M = log2(a$FPKM) - log2(b$FPKM)
    	   A = (1/2) * (log2(a$FPKM) + log2(b$FPKM))
	   #corrValue = cor(a$FPKM, b$FPKM)	
	   mainText = paste("M-A plot:", title1, "~", title2, sep=" ")
           plot(A,M, xlab="A = (1/2) * (log2(FPKM) + log2(FPKM))", ylab="M = log2(FPKM) - log2(FPKM)", main=mainText, ylim=c(-20,20), xlim=c(-10,20))
	  }
	# }
    }
 }
}
dev.off()
fileList = list.files(pattern="*genes")
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    B = strsplit(fileList[i], "_")
    tempName = paste(B[[1]][5], B[[1]][6], sep="_")
    tempName = paste(tempName, B[[1]][7], sep="_")
    tempName = paste(tempName, B[[1]][8], sep="_")
    tempName = paste(tempName, ".csv", sep="")

    # tempName = paste(fileList[i], ".csv", sep="")
    
    write.csv(a, tempName)
}
