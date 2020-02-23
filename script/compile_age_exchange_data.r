wd = "C:\\Users\\b9930\\Google ¶³ºİµwºĞ\\publication\\SpatialVariability\\"


### NS-IBTS
## load data
setwd(paste0(wd, "data\\survey_north_sea"))
load("NSIBTSdata.RData")

# keep variables we want
data = subset(NSIBTSdata, # NSIBTSdata
              select=c(Year, Quarter, Ship, Gear, SweepLngt, AreaType, AreaCode,
                       Age, NoAtALK, SpecCodeType, SpecCode))
#rm(NSIBTSdata)

# keep gear of "GOV"
data = subset(data, subset=data$Gear=="GOV")

# keep area code of "subarea"
data = subset(data, subset=data$AreaType=="0")

# keep year of 1991 ~ latest year
data = subset(data, subset=data$Year >= 1991)

# keep quarter of 1 and 3
data = subset(data, subset=data$Quarter %in% c(1, 3))

# convert species code to species name
setwd(paste0(wd, "data"))
spec_code = read.csv("spec_code_compiled_fill.csv")

data$SpecName = spec_code$Name[match(data$SpecCode, spec_code$SpecCode)]

# keep variable we want
data = subset(data, select=c(Year, Quarter, Ship, SweepLngt, AreaCode, Age, NoAtALK, SpecName))
data$SpecName = as.character(data$SpecName)
row.names(data) = NULL
head(data)


## pre-processing
data = data[!is.na(data$SweepLngt), ]
data$CPUE = data$NoAtALK / data$SweepLngt

# sum up catches per age class, area, quarter, year for each ship
data_agg = with(data, aggregate(
    CPUE, by=list(Year=Year, Quarter=Quarter, Ship=Ship, AreaCode=AreaCode, 
                  Age=Age, SpecName=SpecName), FUN=function(x){sum(x, na.rm=TRUE)}))

# average catches across ships to avoid effects of duplicated sampling
# (average catches per age class, area, quarter, year)
data_agg = with(data_agg, aggregate(
    x, by=list(Year=Year, Quarter=Quarter, AreaCode=AreaCode, 
               Age=Age, SpecName=SpecName), FUN=function(x){mean(x, na.rm=TRUE)}))
head(data_agg)

# remove areas which number of fish is always 0 (no matter age classes)
area_stat = with(data_agg, aggregate(
    x, by=list(AreaCode=AreaCode, SpecName=SpecName),
    FUN=function(x){sum(x, na.rm=TRUE)}
))
head(area_stat)
sum(area_stat$x==0)  # no areas which catch is 0

# To ensure that our analysis is statistical meaningful, we remove time points (quarters)
# when the number of areas is less than 10
library(tidyr)
library(plyr)

data_agg = data_agg %>% spread(Age, x)
n_var = dim(data_agg)[2]
names(data_agg) = c(names(data_agg)[1:4], paste0("Age_", names(data_agg)[5:length(data_agg)]))
head(data_agg)

num_area_per_quart = with(data_agg, aggregate(
    AreaCode, by=list(Year=Year, Quarter=Quarter, SpecName=SpecName),
    FUN=function(x){length(x)}
))
head(num_area_per_quart)
quart_to_drop = num_area_per_quart[num_area_per_quart$x >= 10, ]

head(quart_to_drop)
#sum(quart_to_drop$x)

data_agg_final = match_df(data_agg, quart_to_drop, on=c("Year", "Quarter", "SpecName"))
head(data_agg_final)

# save data frame
write.csv(data_agg_final, 
          file=paste0(wd, "data\\compiled_age_exchange_data.csv"), 
          row.names=FALSE)



### BITS
## load data
setwd(paste0(wd, "data\\survey_outside_north_sea"))
load("BITSdata.RData")

# keep variables we want
data = subset(BITSdata, 
              select=c(Year, Quarter, Ship, Gear, SweepLngt, AreaType, AreaCode,
                       Age, NoAtALK, SpecCodeType, SpecCode))
#rm(BITSdata)

# fix gear
data = subset(data, subset=data$Gear %in% c("TVS", "TVL"))

# keep area code of "subarea"
data = subset(data, subset=data$AreaType=="0")

# keep year of 2001 ~ latest year
data = subset(data, subset=data$Year >= 2001)

# keep quarter of 1 and 4
data = subset(data, subset=data$Quarter %in% c(1, 4))

# convert species code to species name
setwd(paste0(wd, "data"))
spec_code = read.csv("spec_code_compiled_fill.csv")

data$SpecName = spec_code$Name[match(data$SpecCode, spec_code$SpecCode)]

# keep variable we want
data = subset(data, select=c(Year, Quarter, Ship, SweepLngt, AreaCode, Age, NoAtALK, SpecName))
data$SpecName = as.character(data$SpecName)
row.names(data) = NULL
head(data)


## pre-processing
data = data[!is.na(data$SweepLngt), ]
data$CPUE = data$NoAtALK / data$SweepLngt

# sum up catches per age class, area, quarter, year for each ship
data_agg = with(data, aggregate(
    CPUE, by=list(Year=Year, Quarter=Quarter, Ship=Ship, AreaCode=AreaCode, 
                     Age=Age, SpecName=SpecName), FUN=function(x){sum(x, na.rm=TRUE)}))

# average catches across ships to avoid effects of duplicated sampling
# (average catches per age class, area, quarter, year)
data_agg = with(data_agg, aggregate(
    x, by=list(Year=Year, Quarter=Quarter, AreaCode=AreaCode, 
               Age=Age, SpecName=SpecName), FUN=function(x){mean(x, na.rm=TRUE)}))
head(data_agg)

# remove areas which number of fish is always 0 (no matter age classes)
area_stat = with(data_agg, aggregate(
    x, by=list(AreaCode=AreaCode, SpecName=SpecName),
    FUN=function(x){sum(x, na.rm=TRUE)}
))
head(area_stat)
sum(area_stat$x==0)  # no areas which catch is 0

# To ensure that our analysis is statistical meaningful, we remove time points (quarters)
# when the number of areas is less than 10
library(tidyr)
library(plyr)

data_agg = data_agg %>% spread(Age, x)
n_var = dim(data_agg)[2]
names(data_agg) = c(names(data_agg)[1:4], paste0("Age_", names(data_agg)[5:length(data_agg)]))
head(data_agg)

num_area_per_quart = with(data_agg, aggregate(
    AreaCode, by=list(Year=Year, Quarter=Quarter, SpecName=SpecName),
    FUN=function(x){length(x)}
))
head(num_area_per_quart)
quart_to_drop = num_area_per_quart[num_area_per_quart$x >= 10, ]

sum(quart_to_drop$x)
head(quart_to_drop)

data_agg_final = match_df(data_agg, quart_to_drop, on=c("Year", "Quarter", "SpecName"))
head(data_agg_final)

# save data frame
write.csv(data_agg_final, 
          file=paste0(wd, "data\\compiled_age_exchange_data_BITS.csv"), 
          row.names=FALSE)



### BTS
## load data
setwd(paste0(wd, "data\\survey_north_sea"))
load("BTSdata.RData")

# keep variables we want
data = subset(BTSdata, 
              select=c(Year, Quarter, Ship, Gear, SweepLngt, AreaType, AreaCode,
                       Age, NoAtALK, SpecCodeType, SpecCode))
#rm(BTSdata)

# fix gear ("BT4A", "BT4AI", "BT8")
data = subset(data, subset=data$Gear %in% c("BT8"))

# keep area code of "subarea"
data = subset(data, subset=data$AreaType=="0")

# keep quarter of 3
data = subset(data, subset=data$Quarter %in% c(3))

# convert species code to species name
setwd(paste0(wd, "data"))
spec_code = read.csv("spec_code_compiled_fill.csv")

data$SpecName = spec_code$Name[match(data$SpecCode, spec_code$SpecCode)]

# keep variable we want
data = subset(data, select=c(Year, Quarter, Ship, SweepLngt, AreaCode, Age, NoAtALK, SpecName))
data$SpecName = as.character(data$SpecName)
row.names(data) = NULL
head(data)


## pre-processing
# sum up catches per age class, area, quarter, year for each ship
data_agg = with(data, aggregate(
    NoAtALK, by=list(Year=Year, Quarter=Quarter, Ship=Ship, AreaCode=AreaCode, 
                     Age=Age, SpecName=SpecName), FUN=function(x){sum(x, na.rm=TRUE)}))

# average catches across ships to avoid effects of duplicated sampling
# (average catches per age class, area, quarter, year)
data_agg = with(data_agg, aggregate(
    x, by=list(Year=Year, Quarter=Quarter, AreaCode=AreaCode, 
               Age=Age, SpecName=SpecName), FUN=function(x){mean(x, na.rm=TRUE)}))
head(data_agg)

# remove areas which number of fish is always 0 (no matter age classes)
area_stat = with(data_agg, aggregate(
    x, by=list(AreaCode=AreaCode, SpecName=SpecName),
    FUN=function(x){sum(x, na.rm=TRUE)}
))
head(area_stat)
sum(area_stat$x==0)  # no areas which catch is 0

# To ensure that our analysis is statistical meaningful, we remove time points (quarters)
# when the number of areas is less than 10
library(tidyr)
library(plyr)

data_agg = data_agg %>% spread(Age, x)
n_var = dim(data_agg)[2]
names(data_agg) = c(names(data_agg)[1:4], paste0("Age_", names(data_agg)[5:length(data_agg)]))
head(data_agg)

num_area_per_quart = with(data_agg, aggregate(
    AreaCode, by=list(Year=Year, Quarter=Quarter, SpecName=SpecName),
    FUN=function(x){length(x)}
))
head(num_area_per_quart)
quart_to_drop = num_area_per_quart[num_area_per_quart$x >= 10, ]

sum(quart_to_drop$x)
head(quart_to_drop)

data_agg_final = match_df(data_agg, quart_to_drop, on=c("Year", "Quarter", "SpecName"))
head(data_agg_final)

# save data frame
write.csv(data_agg_final, 
          file=paste0(wd, "data\\compiled_age_exchange_data_BTS.csv"), 
          row.names=FALSE)



### Appendix
## keep species of interests
# standard species
target_spec = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus", 
                "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
                "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")
# other species
target_spec = c("Enchelyopus cimbrius", "Lophius piscatorius", "Microstomus kitt",
                "Platichthys flesus", "Trachurus trachurus")

# only Microstomus kitt has enough data
data_agg_final = subset(data_agg, subset=data_agg$SpecName %in% target_spec[5])
head(data_agg_final)
table(data_agg_final$Year, data_agg_final$Quarter)


## developing code
# NS-IBTS
target_spec = c("Glyptocephalus cynoglossus", "Lophius piscatorius", 
                "Merluccius merluccius", "Microstomus kitt", "Scophthalmus maximus")

test = subset(data, subset=data$SpecName %in% target_spec[4])
head(test)
table(test$Year, test$Quarter)

test1 = subset(test, subset=test$Year==2010 & test$Quarter==3)
unique(test1$AreaCode)

test = subset(data_agg_final, subset=data_agg_final$SpecName %in% target_spec[4])
table(test$Year, test$Quarter)


# BITS
target_spec = c("Platichthys flesus", "Scophthalmus maximus", "Limanda limanda",
                "Scophthalmus rhombus", "Psetta maxima")
target_spec = c("Platichthys flesus", "Solea solea")

test = subset(data, subset=data$SpecName %in% target_spec[1])
head(test)
table(test$Year, test$Quarter)

library(dplyr)
test %>%
    group_by(Year, Quarter) %>%
    summarise(n_distinct(AreaCode)) %>% print(n=Inf)

test = subset(data_agg_final, subset=data_agg_final$SpecName %in% target_spec[1])
table(test$Year, test$Quarter)


# BTS
target_spec = c("Pleuronectes platessa")

test = subset(data_agg_final, subset=data_agg_final$SpecName %in% target_spec[1])
table(test$Year, test$Quarter)
