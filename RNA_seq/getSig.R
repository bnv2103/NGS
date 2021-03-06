args <- commandArgs(TRUE)
# fileList = list.files(pattern="*gene_exp.diff")
dataf = args[1]
CompName = args[2]
# for (i in 1:length(fileList) ){
    a = read.table(dataf, header=T)
    tempName = paste(CompName, dataf, sep="_")
    tempName = paste(tempName, ".csv", sep="")
    b = a[a$significant == "yes", ]
    c = a[order(a$q_value),]
    b$Notes = "T"
    for (j in 1:(dim(b)[1]) ){
        if (b$value_1[j] < 1 || b$value_2[j] < 1){
           b$Notes[j] = "F"
        }
    }
    write.csv(b, tempName)

# }
tempName = paste(CompName, ".pdf", sep="")
pdf(file = tempName)
title= paste(CompName, "p_value", sep="-")
hist(a$p_value, xlab="p_value", main=title)
title= paste(CompName, "q_value", sep="-")
hist(a$q_value, xlab="q_value", main=title)
dev.off()