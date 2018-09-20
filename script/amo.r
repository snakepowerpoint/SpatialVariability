setwd("C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\data\\raw")

amo = read.csv("amon.csv", header=F, sep="")
colnames(amo) = c("Year", paste0("M", 1:12))

# compute mean AMO on quarter 1 (January & February) and 3 (August & September)
amo$Q1 = apply(amo[,c(2,3)], 1, mean)
amo$Q3 = apply(amo[,c(9,10)], 1, mean)

# select amo from 1991 to 2015
amo = subset(amo, subset=Year%in%c(1991:2015))
amo = subset(amo, select=c(Year, Q1, Q3))

# convert amo into vector
n = dim(amo)[1]
amo = data.frame(Year=rep(amo$Year,2), Quarter=c(rep(1,n),rep(3,n)), AMO=c(amo$Q1, amo$Q3))
amo = amo[with(amo, order(Year,Quarter)), ]
rownames(amo) = c()

save.image(file="C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\data\\amo.RData")

