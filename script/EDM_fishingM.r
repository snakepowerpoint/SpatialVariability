library(Rcpp)
library(devtools)
library(rEDM)



### Data preparation
wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "output"))

species = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus",
            "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
            "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")

# load data
EDM_list = list()
for (i in 1:length(species)){
    EDM_list[[i]] = list()
    EDM_list[[i]]$species = species[i]
    EDM_list[[i]]$data = read.csv(paste0(species[i], ".csv"), header = TRUE)
}
names(EDM_list) = species

# standarization
EDM_list = lapply(EDM_list, FUN = function(item){
    item$data_std = cbind(item$data[, c(1,2)], scale(item$data[, -c(1,2)]))
    
    return(item)
})

# change column name
EDM_list = lapply(EDM_list, function(item){
    var.name = c("Year", "Quarter", "CV.CPUE", "AgeDiversity", "Abundance", 
                 "AMO", "SBT", "CVofSBT", "SST", "CVofSST")
    names(item$data) = var.name
    names(item$data_std) = var.name
    
    return(item)
})

# drop discontinuous data
EDM_list$`Pleuronectes platessa`$data_std = EDM_list$`Pleuronectes platessa`$data_std[19:48, ]

# a function to compute mean value for each year
compute_mean_per_year = function(data){
    data_wo_year = subset(data, select=-c(Year, Quarter))
    mean_per_year = aggregate(data_wo_year, by=list(Year=data$Year), FUN=mean)
    return(mean_per_year)
}

# apply above function to each item (species) in list
EDM_list = lapply(EDM_list, FUN=function(item){
    data_std = item$data_std
    data_per_year = compute_mean_per_year(data_std)
    item$data_per_year = data_per_year
    return(item)
})

# load fishing mortality data
setwd(paste0(wd, "output"))
fishingM = read.csv("fishingM.csv", header=TRUE)
fishingM_wo_species = subset(fishingM, select=-c(Species))
fishingM_wo_species = split(fishingM_wo_species, f=fishingM$Species)

# standardization
fishingM_wo_species = lapply(fishingM_wo_species, FUN=function(item){
    item$F = as.vector(scale(item$F))
    return(item)
})

# combined data with fishing mortality data
for (i in 1:length(EDM_list)){
    result = merge(EDM_list[[i]]$data_per_year, fishingM_wo_species[[i]], by=c("Year"))
    EDM_list[[i]]$data_per_year = result
}



### a function returning all local peaks
# (will be applied to get all local optimal embedding dimensions)
detectPeak = function(data){
    n = length(data)
    peaks = c()
    rho = c()
    for (idx in 2:(n-1)){
        if ((data[idx-1] < data[idx]) & (data[idx+1] < data[idx])){
            peaks = c(peaks, idx)
            rho = c(rho, data[idx])}
    }
    output = data.frame(peaks, rho)
    
    return(output[order(output$rho, decreasing = TRUE), ])
}

### a function using simplex projection to determine embedding dimension for each variable
computeEmbedding = function(data, show_result=TRUE){
    # create a list to record embedding dimensions
    dim.list = list()
    
    # determine embedding dimension of each variables
    for (column in colnames(data)){
        output <- simplex(data[[column]], E = 1:10, tau = 1, tp = 1) # embedding dimension 
        if (show_result){
            plot(output$E, output$rho, type = "l", main = column, lwd = 2,
                 xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")
        }
        dim.list[[column]] = detectPeak(output$rho)
        #dim.list[[column]] = output[, c("E", "rho")][order(output$rho, decreasing = TRUE), ]
    }
    
    return(dim.list)
}

EDM_list = lapply(EDM_list, FUN = function(item){
    data = subset(item$data_per_year, select=-c(Year))
    item$E = computeEmbedding(data, show_result = FALSE)
    
    return(item)
})



### CCM Analysis
# Fishing mortality -> Shannon age, Total CPUE
library(data.table)

# a function applying CCM to examine (lagged) causal effect of fishing mortality
# on age diversity and abundance 
determineCausality = function(
    data, 
    dim.list, 
    species,
    lags = 3, 
    num_samples = 100, 
    show_result = FALSE
){
    cols = c('length', 'lib.size', 'E', 'library', 'target', 'tar.lag', 
             'rho', 'sd.rho', 'kendall.tau', 'significance')
    cv.var = data.frame(matrix(0, length(lags:0), length(cols)))  # ccm results (cv xmap var)
    colnames(cv.var) = cols
    
    output = data.frame(matrix(0, nrow = 0, ncol = length(cols)))  # ccm results for all variables
    colnames(output) = colnames(cv.var)
    
    E = dim.list$`F`$peak[1]
    max_time = dim(data)[1]
    for (idx in 2:length(data)){
        lib_variable = colnames(data)[idx]
        for (lag in -lags:0){
            cv_x <- ccm(data, E = E, tp = lag, lib_sizes = c(seq(E, max_time, 3), max_time),
                        lib_column = idx, target_column = 1, num_samples = num_samples)
            cv_x_m <- data.frame(ccm_means(cv_x), sd.rho = with(cv_x, tapply(rho, lib_size, sd)))
            
            n = dim(cv_x_m)[1]  # last row: ccm results with maximal library size
            
            current_row = lag + lags + 1
            cv.var[current_row, "length"] = max_time
            cv.var[current_row, "lib.size"] = cv_x_m$lib_size[n]
            cv.var[current_row, "E"] = E
            cv.var[current_row, "library"] = lib_variable
            cv.var[current_row, "target"] = "fishingM"
            cv.var[current_row, "tar.lag"] = lag
            cv.var[current_row, "rho"] = cv_x_m$rho[n]
            cv.var[current_row, "sd.rho"] = cv_x_m$sd.rho[n]
            cv.var[current_row, "kendall.tau"] = cor.test(cv_x_m$lib_size, cv_x_m$rho, method = 'kendall')$p.value
            cv.var[current_row, "significance"] = t.test(subset(cv_x, subset = lib_size==cv_x_m$lib_size[n])$rho,
                                                         alternative = "greater")$p.value
        } 
        output = rbind(output, cv.var) 
    }
    
    return(output)  
}

EDM_list = lapply(EDM_list, function(item){
    data = subset(item$data_per_year, select=c(F, AgeDiversity, Abundance))
    
    ccm_result = determineCausality(data = data, 
                                    dim.list = item$E, 
                                    species = item$species,
                                    show_result = FALSE)
    item$ccm = ccm_result
    
    return(item)
})


# ggplot function to plot CCM results
library(ggplot2)
setwd(paste0(wd, 'script'))
source("utils/plot.r")

# plot CCM results for all species
setwd(paste0(wd, 'output/ccm'))
for (species_list in EDM_list){
    species = species_list$species
    gplot_ccm_result(species_list)
    ggsave(paste0(species, '_fishingM.png'), scale = 2)
}

# (doubled check if needed) plot CCM results for all species
for (item in EDM_list){
    plot_ccm_result(item$ccm, 'AgeDiversity', item$species)
    plot_ccm_result(item$ccm, 'Abundance', item$species)
}



### Appendix
# species name
''' 
Clupea harengus
Gadus morhua
Melanogrammus aeglefinus
Merlangius merlangus
Pleuronectes platessa
Pollachius virens
Scomber scombrus
Sprattus sprattus
Trisopterus esmarkii
'''
