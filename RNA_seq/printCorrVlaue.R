args <- commandArgs(TRUE)

fileList = list.files(pattern="*genes")

# calculate coeff. 

fileName = "corrValue.txt"

title = c(1:length(fileList))
write("", fileName, sep="\t")
for (i in 1:(length(fileList)) ){
    t1 = strsplit(fileList[i], "_")
    title[i] = t1[[1]][4]
    #write(title1, fileName, append = TRUE, sep="\t")
}
write.table(title, fileName, append = FALSE, sep="\t")
#write("", fileName, append = TRUE, sep="\n")
#for (i in 1:1 ){
x = array(dim=c(length(fileList), length(fileList)))
 for (i in 1:(length(fileList)) ){
    #x = c(1:length(fileList), 1:length(fileList))
    for (j in 1:(length(fileList)) ){
        t2 = strsplit(fileList[j], "_")
        title2 = t2[[1]][4]
           a = read.table(fileList[i], header=T)
           b = read.table(fileList[j], header=T)
	   #corrValue = cor(a$FPKM, b$FPKM)
	   a1 = merge(a, b, by="tracking_id")
	   a2 = a1[a1$FPKM.x > 0,]
	   a3 = a2[a2$FPKM.y > 0,]
	   #corrValue = cor(log(a3$FPKM.x), log(a3$FPKM.y))
	   corrValue = cor((a3$FPKM.x), (a3$FPKM.y))
	   x[i,j] = corrValue

    }
    #write.table(x, fileName, append=TRUE, sep="\t")
 }
write.table(x, fileName, append=TRUE, sep="\t")
