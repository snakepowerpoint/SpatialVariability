library(ncdf4)

wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
sst_location = "data\\raw\\SST"
setwd(paste0(wd, sst_location))

### COBE-SST2
nc = nc_open("sst.mon.mean.nc")  # load data information
nc$dim$time$units  #days since 1891-1-1 00:00:00

sst = ncvar_get(nc)  # get variable in the data
lon = nc$dim$lon$vals  # longitude (starts from 0.5E with interval 1)
lat = nc$dim$lat$vals  # latitude (starts from 89.5N with interval 1)
time = nc$dim$time$vals  # time (unit:days since 1891-01-01), #493 is 0, from 1850.01-2015.12

# restrict SST within 1965/01 and 2015/12
time.date = as.Date(time, origin = '1891-1-1')  #convert days to yyyy-mm-dd
time.date = substr(time.date, start=1, stop=4)  #extract year
sst = sst[, , which(time.date %in% 1965:2015)]

# restrict SST in the North Sea (IBTS: 5W-13E, 49.5N-61.5N)
# since the grid size is different, here we limit the satellite data in
# 4.5W-12.5E, 49.5N-61.5N
NS_lat = match(seq(49.5, 61.5, 1), lat)
NS_lon = match(c(seq(355.5, 359.5, 1), seq(0.5, 12.5, 1)), lon)

sst = sst[NS_lon, NS_lat, ]


### Statistics of SST data
sst.stat <- data.frame(Year = rep(1965:2015, each = 12), Month = 1:12)

sst.stat$MeanSST = apply(sst, 3, mean, na.rm=T)
sst.stat$sdSST = apply(sst, 3, sd, na.rm=T)
sst.stat$CVofSST = sst.stat$sdSST/sst.stat$MeanSST

# According to IBTS, quarter 1 refers to January and February
# quarter 3 refers to August and September
sst.stat.q <- data.frame(Year = rep(1965:2015, 2), Quarter = rep(c(1, 3), each = length(1965:2015)))

sst.stat.q$MeanSST = c(
    with(subset(sst.stat, subset = Month%in%c(1,2)), aggregate(MeanSST, by = list(Year=Year), FUN=mean)$x),
    with(subset(sst.stat, subset = Month%in%c(8,9)), aggregate(MeanSST, by = list(Year=Year), FUN=mean)$x)
)

sst.stat.q$sdSST = c(
    with(subset(sst.stat, subset = Month%in%c(1,2)), aggregate(sdSST, by = list(Year=Year), FUN=mean)$x),
    with(subset(sst.stat, subset = Month%in%c(8,9)), aggregate(sdSST, by = list(Year=Year), FUN=mean)$x)
)

sst.stat.q$CVofSST = c(
    with(subset(sst.stat, subset = Month%in%c(1,2)), aggregate(CVofSST, by = list(Year=Year), FUN=mean)$x),
    with(subset(sst.stat, subset = Month%in%c(8,9)), aggregate(CVofSST, by = list(Year=Year), FUN=mean)$x)
)

sst.stat.q = sst.stat.q[with(sst.stat.q, order(Year, Quarter)), ]
sst.stat.q = subset(sst.stat.q, subset = Year%in%c(1991:2015)) # restrict data between 1991 and 2015
sst.stat.q = sst.stat.q[, c("Year", "Quarter", "MeanSST", "CVofSST")]

write.csv(sst.stat.q, file=paste0(wd, "output\\SST.csv"), row.names=FALSE)

