args <- commandArgs(TRUE)



fileList = list.files(pattern="*genes")
a = read.csv(fileList[1], header=T)
a2 = read.csv(fileList[2], header=T)
b = cbind(a$FPKM, a2$FPKM)

for (i in 3:length(fileList) ){
    
    a = read.csv(fileList[i], header=T)
    b = cbind(b, a$FPKM)


}

rownames(b) = a$tracking_id
colnames(b) = fileList
write.csv(b, "FPKM.csv")
