### Data preparation
wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "script"))
source("EDM.r")


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
    num_samples = 100
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
                                    species = item$species)
    item$ccm = ccm_result
    
    return(item)
})


# ggplot function to plot CCM results
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



### CCM analysis to find causal variables of age diversity and abundance
# sort data frame: first are affected var



### Apply S-map on age diversity and abundance


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
