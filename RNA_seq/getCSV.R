args <- commandArgs(TRUE)
fileList = list.files(pattern="*gene_exp.diff")
# fileList = list.files(pattern="*fpkm_tracking")
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=T)
    tempName = paste(fileList[i], ".csv", sep="")
    b = a[a$significant == "yes", ]    
    c = a[order(a$q_value),]
    b$Notes = "T"
    for (j in 1:(dim(b)[1]) ){
    	if (b$value_1[j] < 1 || b$value_2[j] < 1){
   	   b$Notes[j] = "F"
	}
    }   
    write.csv(b, tempName)

}
