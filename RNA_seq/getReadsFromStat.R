args <- commandArgs(TRUE)

fileList = list.files(pattern="*stats")
x = 0
for (i in 1:length(fileList) ){
    a = read.table(fileList[i], header=F)
    B = strsplit(fileList[i], "_")
    # get the smaple name
    title = B[[1]][6]
    b = t(a)
    write.table(b, "ReadsSummary.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE, sep="\t")
    for (h in 2:dim(b)[1]){
    	if (as.integer(b[h,2])<14500000){
	   x = x + 1
	}
    }
}
 write.table(x, "ReadsSummary.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE, sep="\t")
