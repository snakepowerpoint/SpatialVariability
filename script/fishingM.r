library(XML)
library(plyr)
library(dplyr)
library(reshape2)

setwd("C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\data\\stock_assessment")

fishingM = data.frame(matrix(0, nrow=0, ncol=3))
colnames(fishingM) = c("Year", "F", "Species")

loadData = function(filename, species=NULL){
    
    temp = ldply(xmlToList(filename), data.frame)
    data = filter(temp, .id=="FishData")
    idx = sapply(data, is.factor)
    data[,idx] = lapply(data[,idx], function(x) as.numeric(as.character(x)))
    data = subset(data, select=c(Year, F))
    data$Species = species
    
    return(rbind(fishingM, data))
}

fishingM = loadData("herringAutumn7689.xml", "Clupea harengus")
fishingM = loadData("cod8052.xml", "Gadus morhua")
fishingM = loadData("haddock8068.xml", "Melanogrammus aeglefinus")
fishingM = loadData("whiting7483.xml", "Merlangius merlangus")
fishingM = loadData("plaice7445.xml", "Pleuronectes platessa")
fishingM = loadData("saithe8066.xml", "Pollachius virens")
fishingM = loadData("mackerel8120.xml", "Scomber scombrus")
fishingM = loadData("sprat7181.xml", "Sprattus sprattus")
fishingM = loadData("norwaypout7998.xml", "Trisopterus esmarkii")

fishingM = subset(fishingM, subset = Year %in% c(1991:2015))
write.csv(fishingM, file = "..\\..\\output\\fishingM.csv", row.names = FALSE)



#Clupea harengus, herring
#Gadus morhua, cod
#Melanogrammus aeglefinus, haddock
#Merlangius merlangus, whiting
#Pleuronectes platessa, plaice
#Pollachius virens, saithe
#Scomber scombrus, mackerel 
#Sprattus sprattus, sprat
#Trisopterus esmarkii, pout

