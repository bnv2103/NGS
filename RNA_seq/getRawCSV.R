args <- commandArgs(TRUE)
fileList = list.files(pattern="*fpkm_tracking")
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    tempName = paste(fileList[i], ".csv", sep="")
    write.csv(a, tempName)

}
