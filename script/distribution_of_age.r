wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data"))

age = read.csv("compiled_age_data.csv", header = TRUE)


### data processing
# latitude and longitude of each subarea
num = rep(28:52, each = 19)
char = rep(c(paste0('E',5:9), paste0('F',0:9), paste0('G',0:3)), length(28:52))
coordinate = data.frame(Subarea = paste0(num, char), Lat = 0, Lon = 0)
coordinate$Lat = rep(seq(49.75, 61.75, 0.5), each = 19)
coordinate$Lon = rep(seq(-4.5, 13.5, 1), length(28:52))

pos = match(x = age$Subarea, table = coordinate$Subarea)
age$lon = coordinate$Lon[pos]
age$lat = coordinate$Lat[pos]


##### plot proportional symbols
library(ggplot2)
library(maps)

scale_x_longitude <- function(xmin=-180, xmax=180, step=1, ...) {
    xbreaks <- seq(xmin,xmax,step)
    xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*W")), ifelse(x > 0, parse(text=paste0(x,"^o", "*E")),x))))
    return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), ...))
}
scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, ...) {
    ybreaks <- seq(ymin,ymax,step)
    ylabels <- unlist(lapply(ybreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*S")), ifelse(x > 0, parse(text=paste0(x,"^o", "*N")),x))))
    return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), ...))
}  
# plot function
plot.map = function(data, species){
    # map
    mp = fortify(map(fill = TRUE, plot = FALSE))
    
    xmin <- round(min(data$lon)) - 1
    xmax <- round(max(data$lon)) + 1
    ymin <- round(min(data$lat)) - 1
    ymax <- round(max(data$lat)) + 1
    zmin <- min(data[, grep("Age", names(data))])
    zmax <- max(data[, grep("Age", names(data))])
    
    age_class = names(data)[grep("Age", names(data))]
    for (age in age_class){
        tit <- bquote(italic(.(paste0(unique(data$Species), ','))) ~ .(age))
        Amap <- ggplot() + 
            geom_polygon(aes(x = long, y = lat, group = group), data = mp, fill = "grey", colour = "grey") + 
            coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), ratio = 0.8) + 
            theme_bw()
        print(Amap + geom_point(data = data, aes_string(x = "lon", y = "lat", size = age)) + 
                  scale_size(range = c(0,6)) + 
                  labs(size = 'Total CPUE') + xlab("Longitude") + ylab("Latitude") + ggtitle(tit) + 
                  scale_x_longitude(xmin=xmin, xmax=xmax, step=5) +
                  scale_y_latitude(ymin=ymin, ymax=ymax, step=2) +
                  theme(plot.title = element_text(hjust = 0.5, size = 20),
                        axis.title = element_text(face = "bold"),
                        axis.text = element_text(size = 16, colour = "black"),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 18),
                        legend.text = element_text(size = 16),
                        legend.title = element_text(size = 16),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.5)))
        ggsave(filename = paste0(species, age, ".eps"),
               path = paste0(wd, "output\\figures\\suppl"), 
               device = "eps")
    }
}

species = as.character(unique(age$Species))
for (sp in species){
    plot.map(data = subset(age, subset = age$Species == sp), species = sp)
}



### Alternatives
# devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(extrafont)
library(showtext)

xmin = min(age$lon)
xmax = max(age$lon)
ymin = min(age$lat)
ymax = max(age$lat)

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

plot.map = function(data, species){
    tit <- bquote(italic(.(species)))
    # colors = c(1:11)
    colors = c('black', 'red', 'green3', 'blue3', 'cyan', 'chocolate4', 
               'orange', 'grey', 'blueviolet', 'yellow', 'magenta')
    col_to_plot = grep("Age", colnames(data), value = TRUE)
    
    ggplot(data = world) +
        xlab("Longitude") + ylab("Latitude") + ggtitle(tit) +
        geom_sf(color = "black", fill = "gray") +
        coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
        theme(panel.grid.major = element_line(colour = 'transparent'),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, size = 2),
              plot.title = element_text(hjust = 0.5, size = 20),
              axis.text.x = element_text(colour = "black", size = 16),
              axis.text.y = element_text(colour = "black", size = 16),
              axis.title = element_text(size = 18, face = "bold"), 
              legend.position = "none",  # add this when legends are needed
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 16),
              text=element_text(family = 'Arial')) + 
        geom_scatterpie(aes(x=lon, y=lat, group=Subarea),
                        data=data, cols=col_to_plot) +
        scale_fill_manual(values=colors, labels=col_to_plot)
}

data = age_species$`Pleuronectes platessa`
data[is.na(data)] = 0

setEPS()
postscript(paste0(wd, "output\\figures\\suppl\\label.eps"))
showtext_begin() ## call this function after opening a device
print(plot.map(data, 'Pleuronectes platessa'))
dev.off()


# add the Arial font
font_add("Arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

for (species in names(age_species)){
    data = age_species[[species]]
    data[is.na(data)] = 0
    
    setEPS()
    postscript(paste0(wd, "output\\figures\\suppl\\", species, ".eps"))
    showtext_begin() ## call this function after opening a device
    print(plot.map(data, species))
    dev.off()
}

# alternative method to plot
for (species in names(age_species)){
    data = age_species[[species]]
    data[is.na(data)] = 0
    plot.map(data, species)
    ggsave(filename = paste0(species, ".eps"),
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

