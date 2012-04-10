args <- commandArgs(TRUE)
## # dir <- args[1]
# dataf = args[1]
a = read.table("readsName", header=F)


temp = unique(a$V1)
numLines = length(temp)
print( numLines )