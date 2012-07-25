 library(edgeR)
 library(limma)
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
    # par(mfrow = c(2,1))
    hist(log10(a$FPKM), br=brNum, xlim=range(c(-10, 10)),xlab="log10(FPKM)", main = paste("Histogram of" , title))
    # plot(density(log10(a$FPKM)), xlim=range(c(-10, 10)),xlab="log10(FPKM)")
    a.s = a[order(-a$FPKM*a$length), ]
    c = a.s[1:100, ]
    plot(1:100, log(c$FPKM*c$length), ylab="Normalized # Reads: log(FPKM*length)", main=paste("Plot of Top 100 Isoforms", title))
    
   # tempName = paste(B[[1]][5], B[[1]][6], sep="_")
   # tempName = paste(tempName, B[[1]][7], sep="_")
    tempName = paste( B[[1]][6], B[[1]][9], sep="_")
    
    #tempName = paste(fileList[i], ".csv", sep="")
    tempName = paste(tempName, ".csv", sep="")
    write.csv(a, tempName)

}

# qplot(a$FPKM, b$FPKM, ylim=c(-10, 2000), xlim=c(-10,2000))
# print MA plot
if (isMA ==1){
   # if (length(fileList) <=3 ){
   #   par(mfrow = c(length(fileList),length(fileList)))
   # }
   par(mfrow = c(2,1))
for (i in 1:(length(fileList)) ){
    for (j in 1:(length(fileList)) ){
        # par(mfrow= c(length(fileList),length(fileList)))
        if (i < j){
           t1 = strsplit(fileList[i], "_")
           t2 = strsplit(fileList[j], "_")
           title1 = t1[[1]][6]
           title2 = t2[[1]][6]
           title_sub1 = strsplit(title1, "-")
           title_sub2 = strsplit(title2, "-")
        # if (title_sub1[[1]][2] == title_sub2[[1]][2]){
           a = read.table(fileList[i], header=T)
           b = read.table(fileList[j], header=T)
	   mainText = paste("Scatter plot:", title1, "~", title2, sep=" ")
           plot(a$FPKM, b$FPKM,,xlab=title1, ylab=title2, main=mainText, ylim=c(0,2000), xlim=c(0,2000))
           lines(lowess(a$FPKM, b$FPKM), col=2)
           sum = summary(lm(a$FPKM ~ b$FPKM))
           r2 = round(sum$r.squared, digits=3)
           r2Test = substitute(R^2 == r2, list(r2=r2))
           text(800, 800, r2Test, cex=1.5, pos=4, col="blue")
	   }
	   }
}
}

if (isMA ==1){
    if (length(fileList) <=3 ){
      par(mfrow = c(length(fileList),length(fileList)))
      }
      else { 
      par(mfrow = c(2,1))
   }

for (i in 1:(length(fileList)) ){
    for (j in 1:(length(fileList)) ){
    	# par(mfrow= c(length(fileList),length(fileList)))
    	if (i < j){
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

	   tmp=maPlot(a$FPKM, b$FPKM, main=mainText, ylim=c(-20, 20), xlim=c(-10,20))
           lines(lowess(tmp[["A"]],tmp[["M"]]),col=2)
           abline(h=0,col="grey")


           # plot(A,M, xlab="A = (1/2) * (log2(FPKM) + log2(FPKM))", ylab="M = log2(FPKM) - log2(FPKM)", main=mainText, ylim=c(-20,20), xlim=c(-10,20))
	   # lines(lowess(!is.na(A), !is.na(M)), col=2)
	   # abline(h=0,col="grey")
		
    }
 }
}
}

dev.off()


fileList = list.files(pattern="*genes")
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    B = strsplit(fileList[i], "_")
    # tempName = paste(B[[1]][5], B[[1]][6], sep="_")
    # tempName = paste(tempName, B[[1]][7], sep="_")
    # tempName = paste(tempName, B[[1]][8], sep="_")
    tempName = paste( B[[1]][6], B[[1]][9], sep="_")
    tempName = paste(tempName, ".csv", sep="")

    # tempName = paste(fileList[i], ".csv", sep="")
    
    write.csv(a, tempName)
}


