args <- commandArgs(TRUE)


fileList = list.files(pattern="*genes")
a = read.table(fileList[1], header=T)
a2 = read.table(fileList[2], header=T)
b = cbind(a$FPKM, a2$FPKM)

for (i in 3:length(fileList) ){
    
    a = read.table(fileList[i], header=T)
    b = cbind(b, a$FPKM)


}
write.csv(b, "result.csv")
