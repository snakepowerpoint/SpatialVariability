library(akima)
library(ggplot2)
library(ggmap)
library(maps)
library(extrafont)
library(sysfonts)
library(showtext)

wd = "C:\\Users\\b9930\\Google ���ݵw��\\publication\\SpatialVariability\\"

loadfonts(device = "win")
# add the Arial font
font_add("Arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


### load compiled data
setwd(paste0(wd, "output"))
species_list = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus",
                 "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
                 "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")

compiled_data = list()
for (i in 1:length(species_list)){
    compiled_data[[species_list[i]]] = read.csv(paste0(species_list[i], ".csv"))
}


### CPUE data
setwd(paste0(wd, "data"))
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


### Select year & quarter that have the highest and lowest indicator
sp = species_list[2]
data_sp = compiled_data[[sp]]
cpue_sp = cpue[[sp]]

names(data_sp)
variable = "Shannon.age"  # Shannon.age, Total.CPUE, AMO, MeanBT, CVofBT, MeanSST, CVofSST
var_name = "age diversity"  # age diversity, abundance, AMO, temperature, cv of temperature

time_max = data_sp[which.max(data_sp[[variable]]), c("Year", "Quarter")]
time_min = data_sp[which.min(data_sp[[variable]]), c("Year", "Quarter")]

cpue_sp_max = subset(cpue_sp, subset=cpue_sp$Year == time_max$Year & 
                         cpue_sp$Quarter == time_max$Quarter)
cpue_sp_min = subset(cpue_sp, subset=cpue_sp$Year == time_min$Year & 
                         cpue_sp$Quarter == time_min$Quarter)

# normalize to [0,1]
cpue_sp_max$CPUE = cpue_sp_max$CPUE / max(cpue_sp_max$CPUE)
cpue_sp_min$CPUE = cpue_sp_min$CPUE / max(cpue_sp_min$CPUE)

# interpolation
inter_max = interp2xyz(interp(cpue_sp_max$lon, cpue_sp_max$lat, cpue_sp_max$CPUE), data.frame=TRUE)
inter_min = interp2xyz(interp(cpue_sp_min$lon, cpue_sp_min$lat, cpue_sp_min$CPUE), data.frame=TRUE)


### plot
# map
mp = fortify(map(fill = TRUE, plot = FALSE))

# figure parameters
max_val = round(max(data_sp[[variable]]), 2)
min_val = round(min(data_sp[[variable]]), 2)
tit_max = bquote(paste(italic(.(sp)), ', ', .(var_name), " = ", .(max_val)))
tit_min = bquote(paste(italic(.(sp)), ', ', .(var_name), " = ", .(min_val)))

# breaks on x and y axis
scale_x_longitude <- function(xmin=-180, xmax=180, step=1, ...) {
    xbreaks <- seq(xmin, xmax, step)
    xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*W")), ifelse(x > 0, parse(text=paste0(x,"^o", "*E")),x))))
    return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), ...))
}
scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, ...) {
    ybreaks <- seq(ymin, ymax, step)
    ylabels <- unlist(lapply(ybreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*S")), ifelse(x > 0, parse(text=paste0(x,"^o", "*N")),x))))
    return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), ...))
} 

# plot
base = ggplot() +
    coord_fixed(xlim=c(xmin, xmax), ylim=c(ymin, ymax), ratio=1) + 
    scale_x_longitude(xmin=xmin, xmax=xmax, step=5) +
    scale_y_latitude(ymin=ymin, ymax=ymax, step=5) +
    theme(plot.title = element_text(hjust = 0.5, size = 22),
          axis.title = element_blank(),
          #axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 20, colour = "black"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          text = element_text(family='Arial'))

# spatial distribution at the highest age diversity
setEPS(width=8, height=8)
save_path = paste0(wd, "output\\figures\\suppl\\distribution_index\\", sp, "_", variable, "_high.eps")
postscript(save_path)
showtext_begin() ## call this function after opening a device

base + labs(x = "", y = "", fill = "Abundance", title = tit_max) +
    stat_contour(aes(x=inter_max$x, y=inter_max$y, z=inter_max$z, fill = ..level..), 
                 geom="polygon", binwidth=0.005, na.rm=TRUE) +
    scale_fill_gradient(low="white", high="red", breaks=seq(0, 1, 0.5), limits=c(0, 1)) +
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black")

dev.off()

ggsave(filename=paste0(sp, "_", variable, "_high.png"), 
       path=paste0(wd, "output\\figures\\suppl\\distribution_index"), 
       device="png", scale=1.5)

### alternative saving
'''
ggsave(filename=paste0(sp, "_", variable, "_high.eps"), 
       path=paste0(wd, "output\\figures\\suppl\\distribution_index"), 
       device="eps", scale=1.5)
'''

# spatial distribution at the lowest age diversity
setEPS(width=8, height=8)
save_path = paste0(wd, "output\\figures\\suppl\\distribution_index\\", sp, "_", variable, "_low.eps")
postscript(save_path)
showtext_begin() ## call this function after opening a device

base + labs(x = "", y = "", fill = "Abundance", title = tit_min) +
    stat_contour(aes(x=inter_min$x, y=inter_min$y, z=inter_min$z, fill=..level..), 
                 geom="polygon", binwidth=0.005, na.rm=TRUE) +
    scale_fill_gradient(low="white", high="blue", breaks=seq(0, 1, 0.5), limits=c(0, 1)) +
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black")

dev.off()

ggsave(filename=paste0(sp, "_", variable, "_low.png"), 
       path=paste0(wd, "output\\figures\\suppl\\distribution_index"), 
       device="png", scale=1.5)

### alternative saving
'''
ggsave(filename=paste0(sp, "_", variable, "_low.eps"), 
       path=paste0(wd, "output\\figures\\suppl\\distribution_index"), 
       device="eps", scale=1.5)
'''
