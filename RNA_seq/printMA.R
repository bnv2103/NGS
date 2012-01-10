args <- commandArgs(TRUE)
outPDF = args[1]

# print MA plot in pdf format

fileList = list.files(pattern="*")

pdf(file = outPDF)

# print MA plot
for (i in 1:(length(fileList)) ){
    for (j in i:(length(fileList)) ){
    	if (i != j){
	t1 = strsplit(fileList[i], "_")
        t2 = strsplit(fileList[j], "_")
   
	
	 title1 = t1[[1]][4]
         title2 = t2[[1]][4]
	# title_sub1 = strsplit(title1, "-")
	# title_sub2 = strsplit(title2, "-")
	# if (title_sub1[[1]][2] == title_sub2[[1]][2]){
    	   a = read.table(fileList[i], header=T)
    	   b = read.table(fileList[j], header=T)
    	   M = log2(a$FPKM) - log2(b$FPKM)
    	   A = (1/2) * (log2(a$FPKM) + log2(b$FPKM))
	   #corrValue = cor(a$FPKM, b$FPKM)	
	   mainText = paste("M-A plot:", title1, "~", title2, sep=" ")
           plot(A,M, xlab="A = (1/2) * (log2(FPKM) + log2(FPKM))", ylab="M = log2(FPKM) - log2(FPKM)", main=mainText, ylim=c(-20,20), xlim=c(-10,20))
	  }
	#}
    }
 }
dev.off()
