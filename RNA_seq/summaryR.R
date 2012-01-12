#

args <- commandArgs(TRUE)
## # dir <- args[1]
# read isoforms and genes information, 
# print median of FPKM, mean of FPKM and number of isoforms/genes which FPKM > 1 and 0.1

dataf = args[1]
dataf2 = args[2]
outfile = args[3]
a = read.table(dataf, header=T)
b = read.table(dataf2, header=T)


result = c(0,0,0,0,0,0)

result[1] = signif(median(a$FPKM), 3)
result[2] = signif(mean(a$FPKM), 3)

temp = dim(a[a$FPKM > 1, ])
result[3] = temp[1]
temp =  dim(a[a$FPKM > 0.1, ])
result[4] = temp[1]
temp = dim(b[b$FPKM > 1, ])
result[5] = temp[1]
temp =  dim(b[b$FPKM > 0.1, ])
result[6] = temp[1]


print(result, quote=FALSE, row.names=FALSE, col.names=FALSE)
