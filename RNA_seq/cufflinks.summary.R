# 

args <- commandArgs(TRUE)
## # dir <- args[1]
dataf = args[1]

a = read.table(dataf, header=T)

# b = a$FPKM + 1e-200

result = c(0,0,0,0,0)

result[1] = signif(median(a$FPKM), 3)
result[2] = signif(mean(a$FPKM), 3)


a.x = sort(a$FPKM, decreasing=T)

num = dim(a)[1]
q1 = round(num * 0.25)
q3 = round(num * 0.75)

result[3] = signif(dim(a[a$FPKM>0.01,])[1] / num, 3)
result[4] = var(a.x[q1:q3]) ** 0.5 
result[5] = num
#result[6] = signif(dim(a[a$FPKM>0.01,])[1] / num, 3)
print(result)

## var(d.x[1000:27000])
# [1] 119.5404
# > dim(a[a$FPKM > 0.01, ])
# [1] 28738    13
# > dim(b[b$FPKM > 0.01, ])
