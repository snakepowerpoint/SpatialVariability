wd = "C:\\Users\\b9930\\Google ���ݵw��\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data"))

age = read.csv("compiled_age_data.csv", header = TRUE)



### data processing
# align latitude and longitude to each subarea
num = rep(28:52, each = 19)
char = rep(c(paste0('E',5:9), paste0('F',0:9), paste0('G',0:3)), length(28:52))

coordinate = data.frame(Subarea = paste0(num, char), Lat = 0, Lon = 0)
coordinate$Lat = rep(seq(49.75, 61.75, 0.5), each = 19)
coordinate$Lon = rep(seq(-4.5, 13.5, 1), length(28:52))

pos = match(x = age$Subarea, table = coordinate$Subarea)
age$lon = coordinate$Lon[pos]
age$lat = coordinate$Lat[pos]


## prepare map
library(ggplot2)
library(maps)

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

# map
mp = fortify(map(fill = TRUE, plot = FALSE))

xmin <- round(min(age$lon)) - 1
xmax <- round(max(age$lon)) + 1
ymin <- round(min(age$lat)) - 1
ymax <- round(max(age$lat)) + 1

Amap <- ggplot() + 
    geom_polygon(aes(x=long, y=lat, group=group), data=mp, fill="grey", colour="black") + 
    coord_fixed(xlim=c(xmin, xmax), ylim=c(ymin, ymax), ratio=1) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())
#print(Amap)



### plot proportional symbols
# a function to plot spatial distribution for each age class
plot_map_per_age = function(data, species, age, map_plot){
    tit <- bquote(italic(.(paste0(species, ','))) ~ .(age))
    
    age_plot = map_plot + geom_point(data=data, aes_string(x="lon", y="lat", size=age)) + 
        scale_size(range=c(0,6)) + 
        labs(size='Abundance') + xlab("Longitude") + ylab("Latitude") + ggtitle(tit) + 
        scale_x_longitude(xmin=xmin, xmax=xmax, step=2) +
        scale_y_latitude(ymin=ymin, ymax=ymax, step=2) +
        theme(plot.title = element_text(hjust = 0.5, size = 20),
              axis.title = element_text(face = "bold"),
              axis.text = element_text(size = 16, colour = "black"),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              panel.grid.major = element_line(colour = 'transparent'),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
    
    return(age_plot)
}

# split data according to species
age_per_species = split(age, f=age$Species, drop=TRUE)

species = c("Gadus morhua")
age_data = age_per_species[[species]]
age_data_mean = with(age_data, aggregate(
    cbind(Age_0, Age_1, Age_2, Age_3, Age_4, Age_5, Age_6, Age_7, Age_8, Age_9, Age_10),
    by=list(Subarea=Subarea), FUN=mean, na.rm=TRUE))

pos = match(x = age_data_mean$Subarea, table = coordinate$Subarea)
age_data_mean$lon = coordinate$Lon[pos]
age_data_mean$lat = coordinate$Lat[pos]

age_class_name = names(age_data_mean)[grep("Age", names(age_data_mean))]

for (age_class in age_class_name){
    age_plot = plot_map_per_age(age_data_mean, 
                                species=species, 
                                age=age_class, 
                                map_plot=Amap)
    #print(age_plot)
    ggsave(filename=paste0(species, age_class, ".eps"),
           plot=age_plot,
           path=paste0(wd, "output\\figures\\suppl"), 
           device="eps",
           scale=1.5)
}


## Alternatives
# devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(extrafont)
library(showtext)


age_species = split(age, age$Species)  # split data by species

# sum up CPUE
age_species = lapply(age_species, FUN=function(x){
    col_age = grep("Age", colnames(x), value = TRUE)
    aggregate(
        cbind(
            Age_0, Age_1, Age_2, Age_3, Age_4, Age_5, 
            Age_6, Age_7, Age_8, Age_9, Age_10
        ) ~ Subarea + lon + lat,
        x, mean, na.action = NULL
    )
})

# standardize data
age_species = lapply(age_species, FUN=function(data){
    col_age = grep("Age", colnames(data), value = TRUE)
    age_scaled = apply(data[, col_age], 2, FUN=function(x){x/max(x)})
    
    scaled_data = data.frame(cbind(data[, c("Subarea", "lon", "lat")], age_scaled))
    return(scaled_data)
})

# write a function for plotting
# font_import()
loadfonts(device = "win")
world <- ne_countries(scale = "medium", returnclass = "sf")

plot.map = function(data, species, map_plot){
    tit <- bquote(italic(.(species)))
    # colors = c(1:11)
    colors = c('black', 'red', 'green3', 'blue3', 'cyan', 'chocolate4', 
               'orange', 'grey', 'blueviolet', 'yellow', 'magenta')
    col_to_plot = grep("Age", colnames(data), value = TRUE)
    
    age_plot = map_plot + 
        geom_scatterpie(aes(x=lon, y=lat, group=Subarea, r=0.25), data=data, cols=col_to_plot) +
        xlab("Longitude") + ylab("Latitude") + ggtitle(tit) +
        theme(panel.grid.major = element_line(colour = 'transparent'),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
              plot.title = element_text(hjust = 0.5, size = 24),
              axis.text.x = element_text(colour = "black", size = 18),
              axis.text.y = element_text(colour = "black", size = 18),
              axis.title = element_text(size = 18, face = "bold"), 
              legend.position = "none",  # add this when legends are needed
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 16)) + 
        scale_fill_manual(values=colors, labels=col_to_plot) +
        scale_x_longitude(xmin=xmin, xmax=xmax, step=2) +
        scale_y_latitude(ymin=ymin, ymax=ymax, step=2)
    
    return(age_plot)
}

# save example
data = age_species$`Pleuronectes platessa`
data[is.na(data)] = 0

age_plot = plot.map(data, 'Pleuronectes platessa', Amap)
#print(age_plot)
ggsave(filename=paste0('Pleuronectes platessa', ".png"),
       plot=age_plot,
       path=paste0(wd, "output\\figures\\suppl"),
       scale=1.5)

for (species in names(age_species)){
    data = age_species[[species]]
    data[is.na(data)] = 0
    
    age_plot = plot.map(data, species, Amap)
    ggsave(filename=paste0(species, ".eps"),
           plot=age_plot,
           path=paste0(wd, "output\\figures\\suppl"),
           scale=1)
}

# alternative method to plot
for (species in names(age_species)){
    data = age_species[[species]]
    data[is.na(data)] = 0
    plot.map(data, species)
    ggsave(filename = paste0(species, ".png"),
           path = paste0(wd, "output\\figures\\suppl"), 
           device = "eps")
}



### Alternatives
# pie plot on map
library(ggplot2)
library(sf)
library(rgeos)

library(maps)
library(mapplots)

##
image(x = xmin:xmax, y = ymin:ymax, z = outer(1:15, 1:15, "+"), 
      xlab = "lon", ylab = "lat")
map(database = "world", xlim = c(xmin, xmax), ylim = c(ymin, ymax), 
    add = TRUE, bg = 'lightblue')
add.pie(z = age[1, 5:9], x = age$lon[1], y = age$lat[1])




### legacy
# demo code
test = subset(age, subset = age$Species == "Clupea harengus")
mp1 <- fortify(map(fill = TRUE, plot = FALSE))
mp2 <- map_data('world2')
mp3 <- map_data('world')

xmin <- min(test$lon) - 2
xmax <- max(test$lon) + 2
ymin <- min(test$lat) - 2
ymax <- max(test$lat) + 2
zmin <- min(test[, grep("Age", names(test))])
zmax <- max(test[, grep("Age", names(test))])

#plot map
Amap <- ggplot() + 
    geom_polygon(aes(x = long, y = lat, group = group), data = mp1, fill = "grey", colour = "grey") + 
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) + 
    theme_bw()

#ratio = 0.8
Amap <- ggplot() + 
    geom_polygon(aes(x = long, y = lat, group = group), data = mp1, fill = "grey", colour = "grey") + 
    coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), ratio = 0.8) + 
    theme_bw()

#add point
tit <- bquote(italic(.(paste0(unique(test$Species), ','))) ~ 'age 0')
Amap + geom_point(data = test, aes(x = lon, y = lat, size = Age_0)) + 
    scale_size(range = c(0,6)) + 
    labs(size = 'Total CPUE') + ggtitle(tit) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))

