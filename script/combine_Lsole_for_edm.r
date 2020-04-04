wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(wd)

age_cpue_data = read.csv("data\\compiled_age_data_mk.csv", header=T)


### data-preprocessing
## remove subareas never larger than 0
age_cols_name = grep("Age_", names(age_cpue_data), value=T)
subarea_by_age = aggregate(x=age_cpue_data[, age_cols_name],
                           by=list(age_cpue_data$SubArea), FUN=sum, na.rm=T)
subarea_cpue = data.frame(Subarea = subarea_by_age[, 1], 
                          CPUE = apply(subarea_by_age[, -1], 1, sum, na.rm=T))
which(!subarea_cpue$CPUE > 0) 
# no subarea never larger than 0


## remove quarters which number of subareas is less than 10
aggregate(SubArea ~ Year + Quarter, data=age_cpue_data, 
          FUN=function(x){length(unique(x))})
# no quarters which number of subareas is less than 10



### calculate age diversity in each time point
library(dplyr)

cal_age_diversity = function(data){
    # calculate total CPUE for each age class at each time points
    cpue_by_year_quarter = subset(data, select=-c(SubArea, Species)) %>%
        group_by(Year, Quarter) %>%
        summarise_all(sum, na.rm=TRUE)
    
    # calculate age diversity each time points
    age_only = subset(cpue_by_year_quarter, select=-c(Year, Quarter))
    shannon = apply(age_only, 1, FUN=function(x){
        prob = x/sum(x, na.rm=TRUE)
        shannon = sum(-prob * log(prob), na.rm=TRUE)
        return(shannon)
    })
    
    return(cbind(cpue_by_year_quarter[, c("Year", "Quarter")], Shannon=shannon))
}

diversity = cal_age_diversity(age_cpue_data)



### calculate total CPUE
cal_total_cpue = function(data){
    # calculate total CPUE for each subarea per year per quarter
    age_only = subset(data, select=-c(Year, Quarter, SubArea, Species))
    cpue_by_subarea = apply(age_only, 1, function(x){sum(x, na.rm=TRUE)})
    cpue_by_subarea = cbind(data[, c("Year", "Quarter", "SubArea")], TotalCPUE=cpue_by_subarea)
    
    # calculate total CPUE for each year and quarter
    total_cpue = aggregate(TotalCPUE ~ Year + Quarter, 
                           data=cpue_by_subarea, 
                           function(x){sum(x, na.rm=TRUE)})
    return(total_cpue)
}

total_cpue = cal_total_cpue(age_cpue_data)



### calculate spatial variability
cal_spatail_cv = function(data){
    # calculate total CPUE for each subarea per year per quarter
    age_only = subset(data, select=-c(Year, Quarter, SubArea, Species))
    cpue_by_subarea = apply(age_only, 1, function(x){sum(x, na.rm=TRUE)})
    cpue_by_subarea = cbind(data[, c("Year", "Quarter", "SubArea")], TotalCPUE=cpue_by_subarea)
    
    # calculate mean and variance of CPUE for each year and quarter
    cpue_var = with(cpue_by_subarea, aggregate(x=list(VarCPUE=TotalCPUE), 
                                               by=list(Year=Year, Quarter=Quarter), 
                                               FUN=function(x){var(x, na.rm=TRUE)}))
    cpue_mean = with(cpue_by_subarea, aggregate(x=list(MeanCPUE=TotalCPUE), 
                                                by=list(Year=Year, Quarter=Quarter), 
                                                FUN=function(x){mean(x, na.rm=TRUE)}))
    output_data = merge(cpue_var, cpue_mean)
    
    # calculate spatial CV of CPUE for each year and quarter
    output_data$CV.CPUE = sqrt(output_data[["VarCPUE"]]) / output_data[["MeanCPUE"]]
    
    return(subset(output_data, select=c(Year, Quarter, CV.CPUE)))
}

spatial_cv = cal_spatail_cv(age_cpue_data)



### AMO
setwd(paste0(wd, "data\\raw\\amo"))

amo = read.csv("amon.us.data_2019.csv", header=F, sep="", skip=1)
amo = amo[1:72, ]  # remove description

# re-format
colnames(amo) = c("Year", paste0("M", 1:12))
year = as.character(amo$Year)
amo = data.frame(Year=year, apply(amo[, -1], 2, as.numeric))

# compute mean AMO on quarter 1 (January & February) and 3 (August & September)
amo$Q1 = apply(amo[,c(2,3)], 1, mean)
amo$Q3 = apply(amo[,c(9,10)], 1, mean)

# select AMO within study period
amo = subset(amo, subset=Year%in%c(2006:2018))
amo = subset(amo, select=c(Year, Q1, Q3))

# convert amo into vector
n = nrow(amo)
amo = data.frame(Year=rep(amo$Year,2), Quarter=c(rep(1,n),rep(3,n)), AMO=c(amo$Q1, amo$Q3))
amo = amo[with(amo, order(Year,Quarter)), ]
rownames(amo) = NULL



### sea bottom temperature
setwd(paste0(wd, "data\\raw\\bottomT"))

files = list.files()
files_idx = grep("csv", files)

# load and combine data within study period
num_files = length(files_idx)
rawT = read.csv(files[files_idx[1]])
for(idx in 2:num_files){
    temp = read.csv(files[files_idx[idx]])
    rawT = rbind(rawT, temp)
}
rm(temp)

rawT$yyyy.mm.ddThh.mm = as.Date(rawT$yyyy.mm.ddThh.mm)  # Dtype: factor to date
rawT$Year = strftime(rawT$yyyy.mm.ddThh.mm, "%Y")  # year 
rawT$Month = strftime(rawT$yyyy.mm.ddThh.mm, "%m")  # month

rawT$Quarter = NA  
rawT$Quarter[rawT$Month %in% c('01', '02')] = 1  # quarter 1 conformed with IBST
rawT$Quarter[rawT$Month %in% c('08', '09')] = 3  # quarter 3 conformed with IBST
head(rawT)

# keep variables we want
bottomT = subset(rawT, select=c(Cruise, Station, Year, Quarter, Latitude..degrees_north.,
                                Longitude..degrees_east.,Bot..Depth..m.,TEMP..deg.C.))
colnames(bottomT) = c("Cruise", "Station", "Year", "Quarter", 
                      "Latitude", "Longitude", "Depth", "Temperature")  # rename

# keep data on quarter 1 and 3
bottomT = droplevels(subset(bottomT, subset = Quarter %in% c(1,3)))

# average duplicated sample
bottomT$Temperature = as.numeric(bottomT$Temperature)
bottomT = with(bottomT, aggregate(x=Temperature, 
                                  by=list(Cruise=Cruise, Year=Year, Quarter=Quarter,
                                          Latitude=Latitude, Longitude=Longitude, Depth=Depth), 
                                  FUN=function(x){mean(x, na.rm=TRUE)}))
bottomT = bottomT[!is.na(bottomT$x), ]
bottomT = with(bottomT, bottomT[order(Year, Quarter), ])
colnames(bottomT)[ncol(bottomT)] = "Temperature"


## compute mean bottom temperature and its spatial variability
meanT = with(bottomT, aggregate(x=list(MeanT=Temperature), 
                                by=list(Year=Year, Quarter=Quarter), 
                                FUN=function(x){mean(x, na.rm=TRUE)}))
cvT = with(bottomT, aggregate(x=list(CVT=Temperature), 
                              by=list(Year=Year, Quarter=Quarter), 
                              FUN=function(x){sqrt(var(x, na.rm=TRUE))/mean(x, na.rm=TRUE)}))

bottomT = merge(meanT, cvT)
colnames(bottomT) = c("Year", "Quarter", "MeanBT", "CVofBT")
bottomT = bottomT[order(bottomT$Year), ]



### merge all data
final_data = merge(spatial_cv, diversity)
final_data = merge(final_data, total_cpue)
final_data = merge(final_data, amo)
final_data = merge(final_data, bottomT)

write.csv(final_data, file=paste0(wd, "output\\Microstomus kitt.csv"), row.names=FALSE)
