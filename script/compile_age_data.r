wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data\\raw"))

age = read.csv("CPUE_per_age_per_subarea.csv", header=T)
# keep variables we want
age = subset(age, select=c(Year, Quarter, Area, Subarea, Species, 
                           Age_0, Age_1, Age_2, Age_3, Age_4, Age_5,
                           Age_6, Age_7, Age_8, Age_9, Age_10))
# use data on quarter 1 and 3
age = droplevels(subset(age, subset=Quarter%in%c(1,3)))

# split data according to species
age.species = split(age, f=age$Species) # 9 standard species

# the maximal number of subarea among quarters is ...
max(unlist(lapply(age.species, FUN=function(x){
    data = data.frame(x)
    max(with(data, table(Year, Quarter)))
}))) # 183

# therefore, we drop out quarters which subareas are less than 10 (5% of maximum)
age.species <- lapply(age.species, FUN=function(x){
    data1 = droplevels(subset(data.frame(x), subset=Quarter%in%1))
    data3 = droplevels(subset(data.frame(x), subset=Quarter%in%3))
    # for Q1
    yr.table1 = table(data1$Year) # summarize how many subareas in Q1 for each year
    yr.keep1 = names(yr.table1)[which(yr.table1>=10)] # keep years which subareas >= 10 
    data1 = droplevels(subset(data1, subset=Year%in%yr.keep1))
    # for Q3
    yr.table3 = table(data3$Year)
    yr.keep3 = names(yr.table3)[which(yr.table3>=10)] 
    data3 = droplevels(subset(data3, subset=Year%in%yr.keep3))
    # combind Q1 and Q3
    result = rbind(data1, data3)
    result[with(result, order(Year, Quarter, Subarea)), ]
})

# restrict age data in 1991-2015
age.species = lapply(age.species, FUN=function(x){
    data = data.frame(x)
    subset(x, subset = Year%in%c(1991:2015))
})

# delete subareas never large than 0
age.subarea = lapply(age.species, FUN=function(x){
    data = data.frame(x)
    with(data, aggregate(cbind(Age_0, Age_1, Age_2, Age_3, Age_4, Age_5, Age_6, Age_7, Age_8, Age_9, Age_10),
                         by = list(Subarea), FUN = sum))
})

age.subarea = lapply(age.subarea, FUN=function(x){
    data = data.frame(x)
    data.frame(Subarea = data[, 1], TotalAgeCPUE = apply(data[, -1], 1, sum, na.rm=T))
})

age.subarea.non0 = lapply(age.subarea, FUN=function(x){
    data = data.frame(x)
    droplevels(data$Subarea[data$TotalAgeCPUE > 0])
})

age.species.drop0 = list()
for (i in 1:length(age.species)){
    age.species.drop0[[i]] = subset(data.frame(age.species[[i]]), 
                                    subset = Subarea %in% age.subarea.non0[[i]])
}
names(age.species.drop0) = names(age.species)

age.species.drop0 = lapply(age.species.drop0, FUN=function(x){
    data = data.frame(x)
    data$Subarea = droplevels(data$Subarea)
    data$Area = NULL
    data
})

# save compiled age data as a csv file
age.data.final = do.call(rbind, age.species.drop0)
rownames(age.data.final) = NULL
age.data.final$Subarea = as.character(age.data.final$Subarea)
age.data.final$Subarea = format(age.data.final$Subarea, scientific = FALSE) 
age.data.final$Species = as.character(age.data.final$Species)

write.csv(age.data.final, file=paste0(wd, "data\\compiled_age_data.csv"), row.names=FALSE)






### Legacy
# sum up CPUE of each age class according to year and quarter
age.species.drop0 = lapply(age.species.drop0, FUN=function(x){
    data = data.frame(x)
    output = with(data, aggregate(cbind(Age_0, Age_1, Age_2, Age_3, Age_4, Age_5, Age_6, Age_7, Age_8, Age_9, Age_10),
                                  by=list(Year=Year, Quarter=Quarter), FUN=sum))
    output = with(output, output[order(Year, Quarter), ])
    
    return(output)
})

age.species.drop0 = lapply(age.species.drop0, FUN=function(x){
    data = data.frame(x)
    Shannon.age = apply(data[, 3:13], 1, FUN=function(x){
        p = x/sum(x, na.rm=T)
        sum(-p*log(p), na.rm=T)
    })
    cbind(data, Shannon.age)
})


rm(age, age.subarea, age.subarea.non0, age.species)
save.image(file="C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\data\\age.RData")


