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

# plot function
plot.map = function(data, species){
    # map
    mp = fortify(map(fill = TRUE, plot = FALSE))
    
    xmin <- min(data$lon) - 2
    xmax <- max(data$lon) + 2
    ymin <- min(data$lat) - 2
    ymax <- max(data$lat) + 2
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
                  labs(size = 'Total CPUE') + ggtitle(tit) + 
                  theme(plot.title = element_text(hjust = 0.5, size = 18),
                        axis.title = element_text(face = "bold"),
                        axis.text = element_text(size = 16, colour = "black"),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 18),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
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

