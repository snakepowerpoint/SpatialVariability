library(akima)
library(ggplot2)
library(ggmap)
library(maps)

wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data"))


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


### CPUE data
cpue = read.csv("compiled_cpue_data.csv", header = TRUE)

# align latitude and longitude to each subarea
num = rep(28:52, each = 19)
char = rep(c(paste0('E',5:9), paste0('F',0:9), paste0('G',0:3)), length(28:52))

coordinate = data.frame(Subarea = paste0(num, char), Lat = 0, Lon = 0)
coordinate$Lat = rep(seq(49.75, 61.75, 0.5), each = 19)
coordinate$Lon = rep(seq(-4.5, 13.5, 1), length(28:52))

pos = match(x = cpue$Subarea, table = coordinate$Subarea)
cpue$lon = coordinate$Lon[pos]
cpue$lat = coordinate$Lat[pos]

xmin <- round(min(cpue$lon)) - 1
xmax <- round(max(cpue$lon)) + 1
ymin <- round(min(cpue$lat)) - 1
ymax <- round(max(cpue$lat)) + 1

cpue = split(cpue, f = cpue$Species)


### Select year & quarter that have the highest and lowest age diversity
species_list = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus",
                 "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
                 "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")

sp = species_list[1]
age_sp = age[[sp]]
cpue_sp = cpue[[sp]]

time_max_diversity = age_sp[which.max(age_sp$Shannon.age), c("Year", "Quarter")]
time_min_diversity = age_sp[which.min(age_sp$Shannon.age), c("Year", "Quarter")]

cpue_sp_max = subset(cpue_sp, subset=cpue_sp$Year == time_max_diversity$Year & 
                         cpue_sp$Quarter == time_max_diversity$Quarter)
cpue_sp_min = subset(cpue_sp, subset=cpue_sp$Year == time_min_diversity$Year & 
                         cpue_sp$Quarter == time_min_diversity$Quarter)

# normalize to [0, 1]
cpue_sp_max$CPUE = cpue_sp_max$CPUE / max(cpue_sp_max$CPUE)
cpue_sp_min$CPUE = cpue_sp_min$CPUE / max(cpue_sp_min$CPUE)

# interpolation
inter_max = interp2xyz(interp(cpue_sp_max$lon, cpue_sp_max$lat, cpue_sp_max$CPUE), data.frame=TRUE)
inter_min = interp2xyz(interp(cpue_sp_min$lon, cpue_sp_min$lat, cpue_sp_min$CPUE), data.frame=TRUE)


### plot
# map
mp = fortify(map(fill = TRUE, plot = FALSE))

# figure parameters
max_diver = round(max(age_sp$Shannon.age), 4)
min_diver = round(min(age_sp$Shannon.age), 4)
tit_max = bquote(italic(.(paste0(sp, ','))) ~ "age diversity =" ~ .(max_diver))
tit_min = bquote(italic(.(paste0(sp, ','))) ~ "age diversity =" ~ .(min_diver))

base = ggplot() +
    coord_fixed(xlim=c(xmin, xmax), ylim=c(ymin, ymax), ratio=1) + 
    labs(x = "Longitude", y = "Latitude", fill = "CPUE") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(size = 16, colour = "black"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          panel.grid.major = element_line(colour = 'transparent'),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

# spatial distribution at the highest age diversity
base + labs(title = tit_max) +
    stat_contour(aes(x=inter_max$x, y=inter_max$y, z=inter_max$z, fill = ..level..), 
                 geom="polygon", binwidth=0.005, na.rm=TRUE) +
    scale_fill_gradient(low="white", high="red", breaks=seq(0, 1, 0.5), limits=c(0, 1)) +
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black")
ggsave(filename=paste0(sp, "_high.eps"), path=paste0(wd, "output\\figures\\suppl\\diversity"), 
       device="eps", scale=1.5)
ggsave(filename=paste0(sp, "_high.png"), path=paste0(wd, "output\\figures\\suppl\\diversity"), 
       device="png", scale=1.5)

# spatial distribution at the lowest age diversity
base + labs(title = tit_min) + 
    stat_contour(aes(x=inter_min$x, y=inter_min$y, z=inter_min$z, fill=..level..), 
                 geom="polygon", binwidth=0.005, na.rm=TRUE) +
    scale_fill_gradient(low="white", high="blue", breaks=seq(0, 1, 0.5), limits=c(0, 1)) +
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black")
ggsave(filename=paste0(sp, "_low.eps"), path=paste0(wd, "output\\figures\\suppl\\diversity"), 
       device="eps", scale=1.5)
ggsave(filename=paste0(sp, "_low.png"), path=paste0(wd, "output\\figures\\suppl\\diversity"), 
       device="png", scale=1.5)



### Appendix
## map
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-20, 59), ylim = c(35, 71), asp = 1)

## filled.contour
filled.contour(x,y,z, plot.axes={axis(1); axis(2); map(add=TRUE, interior=FALSE)} )

## ggplot 
# use original data instead of interpolated data
base + 
    stat_contour(aes(x=cpue_sp_max$lon, y=cpue_sp_max$lat, z=cpue_sp_max$CPUE, fill = ..level..), 
                 geom="polygon", binwidth=0.01, na.rm=TRUE) +
    scale_fill_gradientn(colors = c("white", "red")) +
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black")

base + 
    stat_contour(aes(x=cpue_sp_min$lon, y=cpue_sp_min$lat, z=cpue_sp_min$CPUE, fill = ..level..), 
                 geom="polygon", binwidth=0.01, na.rm=TRUE) +
    scale_fill_gradientn(colors = c("white", "blue")) +
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black")

# scale_fill_gradientn
base + labs(title = tit_min) + 
    stat_contour(aes(x=inter_min$x, y=inter_min$y, z=inter_min$z, fill=..level..), 
                 geom="polygon", binwidth=0.005, na.rm=TRUE) +
    scale_fill_gradientn(colors=c("white", "blue"), breaks=seq(0, 1, 0.5), limits=c(0, 1)) +
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black")
