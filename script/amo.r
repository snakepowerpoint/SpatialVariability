wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data\\raw\\AMO"))

amo = read.csv("amon.csv", header=F, sep="")
colnames(amo) = c("Year", paste0("M", 1:12))
year = as.character(amo$Year)
amo = data.frame(Year=year, apply(amo[, -1], 2, as.numeric))

# compute mean AMO on quarter 1 (January & February) and 3 (August & September)
amo$Q1 = apply(amo[,c(2,3)], 1, mean)
amo$Q3 = apply(amo[,c(9,10)], 1, mean)

# select amo within study period
amo = subset(amo, subset=Year%in%c(1991:2015))
amo = subset(amo, select=c(Year, Q1, Q3))

# convert amo into vector
n = dim(amo)[1]
amo = data.frame(Year=rep(amo$Year,2), Quarter=c(rep(1,n),rep(3,n)), AMO=c(amo$Q1, amo$Q3))
amo = amo[with(amo, order(Year,Quarter)), ]
rownames(amo) = c()

write.csv(amo, file=paste0(wd, "output\\amo.csv"), row.names=FALSE)

