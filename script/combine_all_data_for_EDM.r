wd = "C:\\Users\\b9930\\Google ���ݵw��\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data\\"))



### CPUE data
cpue = read.csv("compiled_cpue_data.csv", header = TRUE)
cpue = split(cpue, f = cpue$Species)

# compute total CPUE and mean and variance of CPUE
cpue <- lapply(cpue, FUN=function(y){
    data = data.frame(y)
    cpue.var = with(data, aggregate(CPUE, by=list(Year=Year, Quarter=Quarter), FUN=function(z){var(z)}))
    cpue.mean = with(data, aggregate(CPUE, by=list(Year=Year, Quarter=Quarter), FUN=function(z){mean(z)}))
    cpue.sum = with(data, aggregate(CPUE, by=list(Year=Year, Quarter=Quarter), FUN=sum))
    
    out = cbind(cpue.var, cpue.mean$x, cpue.sum$x)
    
    # change the column names and sort the data
    colnames(out) = c("Year", "Quarter", "Var.CPUE", "Mean.CPUE", "Total.CPUE")
    out[with(out, order(Year, Quarter)), ]
})


# spatial coefficient of variation (CV)
CV <- lapply(cpue, FUN=function(y){
    data = data.frame(y)
    data$CV.CPUE = sqrt(data$Var.CPUE)/data$Mean.CPUE
    data
})


### Age data
age = read.csv("compiled_age_data.csv", header = TRUE)
age = split(age, f = age$Species)

# sum up CPUE of each age class according to year and quarter
age = lapply(age, FUN=function(x){
    data = data.frame(x)
    output = with(data, aggregate(cbind(Age_0, Age_1, Age_2, Age_3, Age_4, Age_5, Age_6, Age_7, Age_8, Age_9, Age_10),
                                  by=list(Year=Year, Quarter=Quarter), FUN=sum))
    output = with(output, output[order(Year, Quarter), ])
    
    return(output)
})

# age diversity
age = lapply(age, FUN=function(x){
    data = data.frame(x)
    Shannon.age = apply(data[, 3:13], 1, FUN=function(x){
        p = x/sum(x, na.rm=T)
        sum(-p*log(p), na.rm=T)
    })
    cbind(data, Shannon.age)
})

# merge cpue data and age data according to species
result.ch = merge(x=CV[[1]], y=age[[1]], by=c("Year", "Quarter"))
result.gm = merge(x=CV[[2]], y=age[[2]], by=c("Year", "Quarter"))
result.ma = merge(x=CV[[3]], y=age[[3]], by=c("Year", "Quarter"))
result.mm = merge(x=CV[[4]], y=age[[4]], by=c("Year", "Quarter"))
result.pp = merge(x=CV[[5]], y=age[[5]], by=c("Year", "Quarter"))
result.pv = merge(x=CV[[6]], y=age[[6]], by=c("Year", "Quarter"))
result.ss = merge(x=CV[[7]], y=age[[7]], by=c("Year", "Quarter"))
result.ssp = merge(x=CV[[8]], y=age[[8]], by=c("Year", "Quarter"))
result.te = merge(x=CV[[9]], y=age[[9]], by=c("Year", "Quarter"))

# keep key variables
result.ch = subset(result.ch, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.gm = subset(result.gm, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.ma = subset(result.ma, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.mm = subset(result.mm, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.pp = subset(result.pp, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.pv = subset(result.pv, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.ss = subset(result.ss, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.ssp = subset(result.ssp, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))
result.te = subset(result.te, select=c(Year, Quarter, CV.CPUE, Shannon.age, Total.CPUE))


### AMO
setwd(paste0(wd, "output\\"))

amo = read.csv("amo.csv", header = TRUE)

# merge
result.ch = merge(result.ch, amo, by=c("Year", "Quarter"))
result.gm = merge(result.gm, amo, by=c("Year", "Quarter"))
result.ma = merge(result.ma, amo, by=c("Year", "Quarter"))
result.mm = merge(result.mm, amo, by=c("Year", "Quarter"))
result.pp = merge(result.pp, amo, by=c("Year", "Quarter"))
result.pv = merge(result.pv, amo, by=c("Year", "Quarter"))
result.ss = merge(result.ss, amo, by=c("Year", "Quarter"))
result.ssp = merge(result.ssp, amo, by=c("Year", "Quarter"))
result.te = merge(result.te, amo, by=c("Year", "Quarter"))


### Sea bottom temperature
setwd(paste0(wd, "output\\"))

bottomT = read.csv("bottomT.csv", header = TRUE)

# merge
result.ch = merge(result.ch, bottomT, by=c("Year", "Quarter"))
result.gm = merge(result.gm, bottomT, by=c("Year", "Quarter"))
result.ma = merge(result.ma, bottomT, by=c("Year", "Quarter"))
result.mm = merge(result.mm, bottomT, by=c("Year", "Quarter"))
result.pp = merge(result.pp, bottomT, by=c("Year", "Quarter"))
result.pv = merge(result.pv, bottomT, by=c("Year", "Quarter"))
result.ss = merge(result.ss, bottomT, by=c("Year", "Quarter"))
result.ssp = merge(result.ssp, bottomT, by=c("Year", "Quarter"))
result.te = merge(result.te, bottomT, by=c("Year", "Quarter"))


### Sea surface temperature
setwd(paste0(wd, "output\\"))

sst = read.csv("sst.csv", header = TRUE)

# merge
result.ch = merge(result.ch, sst, by=c("Year", "Quarter"))
result.gm = merge(result.gm, sst, by=c("Year", "Quarter"))
result.ma = merge(result.ma, sst, by=c("Year", "Quarter"))
result.mm = merge(result.mm, sst, by=c("Year", "Quarter"))
result.pp = merge(result.pp, sst, by=c("Year", "Quarter"))
result.pv = merge(result.pv, sst, by=c("Year", "Quarter"))
result.ss = merge(result.ss, sst, by=c("Year", "Quarter"))
result.ssp = merge(result.ssp, sst, by=c("Year", "Quarter"))
result.te = merge(result.te, sst, by=c("Year", "Quarter"))


# save data
species = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus",
            "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
            "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")

write.csv(result.ch, file=paste0(wd, "output\\", species[1], ".csv"), row.names=FALSE)
write.csv(result.gm, file=paste0(wd, "output\\", species[2], ".csv"), row.names=FALSE)
write.csv(result.ma, file=paste0(wd, "output\\", species[3], ".csv"), row.names=FALSE)
write.csv(result.mm, file=paste0(wd, "output\\", species[4], ".csv"), row.names=FALSE)
write.csv(result.pp, file=paste0(wd, "output\\", species[5], ".csv"), row.names=FALSE)
write.csv(result.pv, file=paste0(wd, "output\\", species[6], ".csv"), row.names=FALSE)
write.csv(result.ss, file=paste0(wd, "output\\", species[7], ".csv"), row.names=FALSE)
write.csv(result.ssp, file=paste0(wd, "output\\", species[8], ".csv"), row.names=FALSE)
write.csv(result.te, file=paste0(wd, "output\\", species[9], ".csv"), row.names=FALSE)


