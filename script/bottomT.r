wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data\\raw\\bottomT"))

files = list.files()

# load and combine data from 1991 to 2015
rawT = read.csv(files[1])
for(idx in 2:length(files)){
    temp = read.csv(files[idx])
    rawT = rbind(rawT, temp)
}
rm(temp)

rawT$yyyy.mm.ddThh.mm = as.Date(rawT$yyyy.mm.ddThh.mm)  # Dtype: factor to date
rawT$Year = strftime(rawT$yyyy.mm.ddThh.mm, "%Y")  # year 
rawT$Month = strftime(rawT$yyyy.mm.ddThh.mm, "%m")  # month

rawT$Quarter = NA  # quarter
rawT$Quarter[rawT$Month %in% c('01', '02')] = 1
rawT$Quarter[rawT$Month %in% c('08', '09')] = 3

# keep variables we want
bottomT = subset(rawT, select=c(Cruise, Station, Year, Quarter, Latitude..degrees_north.,
                                Longitude..degrees_east.,Bot..Depth..m.,TEMP..deg.C.))
colnames(bottomT) = c("Cruise", "Station", "Year", "Quarter", 
                      "Latitude", "Longitude", "Depth", "Temperature")  # rename

# keep data on quarter 1 and 3
bottomT = droplevels(subset(bottomT, subset = Quarter %in% c(1,3)))

## average duplicated sample
bottomT = with(bottomT, aggregate(Temperature, by=list(Cruise = Cruise, Year=Year, Quarter=Quarter,
                                                       Latitude = Latitude, Longitude = Longitude, Depth = Depth), 
                                  FUN=mean, na.rm=TRUE))
bottomT = bottomT[!is.na(bottomT$x), ]
bottomT = with(bottomT, bottomT[order(Year, Quarter), ])
colnames(bottomT)[dim(bottomT)[2]] = "Temperature"


## compute mean bottom temperature and its spatial variability
meanT = with(bottomT, aggregate(Temperature, by=list(Year=Year, Quarter=Quarter), FUN=mean, na.rm=TRUE))
cvT = with(bottomT, aggregate(Temperature, by=list(Year=Year, Quarter=Quarter), 
                              FUN=function(x){sqrt(var(x, na.rm=TRUE))/mean(x, na.rm=TRUE)}))

bottomT = cbind(meanT, cvT$x)
colnames(bottomT) = c("Year", "Quarter", "MeanBT", "CVofBT")
bottomT = bottomT[order(bottomT$Year), ]

write.csv(bottomT, file=paste0(wd, "output\\bottomT.csv"), row.names=FALSE)
