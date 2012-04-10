args <- commandArgs(TRUE)
fileList = list.files(pattern="*genes_nonRef")
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    B = strsplit(fileList[i], "_")
    tempName = paste( B[[1]][6], B[[1]][9], sep="_")
    tempName = paste(tempName, "_nonRef.csv", sep="")


    write.csv(a, tempName)
}
fileList = list.files(pattern="*isoforms_nonRef")
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    B = strsplit(fileList[i], "_")
    tempName = paste( B[[1]][6], B[[1]][9], sep="_")
    tempName = paste(tempName, "_nonRef.csv", sep="")


    write.csv(a, tempName)
}
