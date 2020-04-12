library(Rcpp)
library(devtools)
library(rEDM)
library(viridis)



### Data preparation
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

# standarization (data)
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

# drop discontinuous data (Pleuronectes platessa has missing year)
EDM_list$`Pleuronectes platessa`$data_std = EDM_list$`Pleuronectes platessa`$data_std[19:48, ]

# load fishing mortality data
setwd(paste0(wd, "output"))
fishingM = read.csv("fishingM.csv", header=TRUE)
fishingM_wo_species = subset(fishingM, select=-c(Species))
fishingM_wo_species = split(fishingM_wo_species, f=fishingM$Species)


## convert quarterly data to yearly data to meet with fishing mortality data
# a function to compute mean value for each year
compute_mean_per_year = function(data){
    data_wo_year = subset(data, select=-c(Year, Quarter))
    mean_per_year = aggregate(data_wo_year, by=list(Year=data$Year), FUN=mean)
    return(mean_per_year)
}

# apply above function to each item (species) in list
EDM_list = lapply(EDM_list, FUN=function(item){
    data = item$data
    data_per_year = compute_mean_per_year(data)
    item$data_per_year = data_per_year
    return(item)
})

# combine data with fishing mortality data
for (i in 1:length(EDM_list)){
    result = merge(EDM_list[[i]]$data_per_year, fishingM_wo_species[[i]], by=c("Year"))
    EDM_list[[i]]$data_per_year = result
}

# standarization (data_per_year)
EDM_list = lapply(EDM_list, FUN = function(item){
    item$data_per_year = cbind(item$data_per_year[c(1)], 
                               scale(item$data_per_year[, -c(1)]))
    return(item)
})

rm(result)
rm(fishingM)
rm(fishingM_wo_species)


## detrend 
# lienar regression to remove significant trend
detrend_sig = function(series){
    data = data.frame(x=1:length(series), y=series)
    fit = lm(y~x, data=data)
    fit_summary = summary(fit)
    p_value = fit_summary$coefficient[, 4]['x']
    
    if (p_value < 0.05){
        return(as.numeric(series - fit$fitted.values))
    } else {
        return(series)
    }
}



### Simplex projection 
# a function returning all local peaks 
# (will be applied to get all local optimal embedding dimensions)
detectPeak = function(data){
    n = length(data)
    peaks = c()
    rho = c()
    for (idx in 2:(n-1)){ # from E=2
        if ((data[idx-1] < data[idx]) & (data[idx+1] < data[idx])){
            peaks = c(peaks, idx)
            rho = c(rho, data[idx])}
    }
    output = data.frame(peaks, rho)
    
    if (nrow(output) < 1){
        rho = max(data[-1]) # ignore E=1
        peaks = which(data == rho)
        return(data.frame(peaks, rho))
    } else {
        return(output[order(output$rho, decreasing = TRUE), ])
    }
}

# a function using simplex projection to determine embedding dimension for each variable
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



### CCM
library(data.table)

# a function using CCM to determine causal variables (including lagged terms) of spatial CV
determineCausality = function(data, dim.list, species, lags=8, num_samples=100, seed=1234){
    cols = c('length', 'lib.size', 'E', 'library', 'target', 'tar.lag', 
             'rho', 'sd.rho', 'kendall.tau', 'significance')
    ccm_var = data.frame(matrix(0, length(lags:0), length(cols)))  # ccm results
    colnames(ccm_var) = cols
    
    output = data.frame(matrix(0, nrow = 0, ncol = length(cols)))  # ccm results for all variables
    colnames(output) = cols
    
    library_var = colnames(data)[1]
    max_time = dim(data)[1]  # number of data points
    for (idx in 2:length(data)){
        target_var = colnames(data)[idx]  # target variable
        E = dim.list[[target_var]]$peak[1]  # E of target variable
        
        for (lag in -lags:0){
            cv_x <- ccm(data, E = E, tp = lag, lib_sizes = c(seq(E, max_time, 3), max_time),
                        lib_column = 1, target_column = idx, num_samples = num_samples, RNGseed=seed)
            cv_x_m <- data.frame(ccm_means(cv_x), sd.rho = with(cv_x, tapply(rho, lib_size, sd)))
            
            n = dim(cv_x_m)[1]  # last row: ccm results with maximal library size
            
            current_row = lag + lags + 1
            ccm_var[current_row, "length"] = max_time
            ccm_var[current_row, "lib.size"] = cv_x_m$lib_size[n]
            ccm_var[current_row, "E"] = E
            ccm_var[current_row, "library"] = library_var
            ccm_var[current_row, "target"] = names(data)[idx]
            ccm_var[current_row, "tar.lag"] = lag
            ccm_var[current_row, "rho"] = cv_x_m$rho[n]
            ccm_var[current_row, "sd.rho"] = cv_x_m$sd.rho[n]
            ccm_var[current_row, "kendall.tau"] = cor.test(cv_x_m$lib_size, cv_x_m$rho, 
                                                           alternative = 'greater', method = 'kendall')$p.value
            ccm_var[current_row, "significance"] = t.test(subset(cv_x, subset = lib_size==cv_x_m$lib_size[n])$rho,
                                                          alternative = 'greater')$p.value
        } 
        output = rbind(output, ccm_var)
    }
    
    return(output) 
}



### S-map
## generate a data frame corresponding to lagged variables for S-map analysis 
generateSmapData = function(item, data, lib_var, is_full=FALSE){
    num_Smap_model = length(item$E_feasible$peaks)
    if (num_Smap_model < 1){
        item$Smap1 = NA
        return(item)
    }
    
    name_Smap_model = paste0("Smap", 1:num_Smap_model)
    num_sample = dim(data)[1]
    num_variable = item$E_feasible$peaks
    
    ccm_result = item$ccm_most_significant
    ccm_result = ccm_result[order(ccm_result$rho, decreasing = TRUE), ]
    
    for (i in 1:num_Smap_model){
        item[[name_Smap_model[i]]] = list()
        item[[name_Smap_model[i]]]$data = data.frame(matrix(0, nrow = num_sample, ncol = 0))
        item[[name_Smap_model[i]]]$data[, lib_var] = data[lib_var]
        
        for (idx_var in 1:(num_variable[i]-1)){
            if (is_full){
                vars = as.character(ccm_result$target[idx_var])
                item[[name_Smap_model[i]]]$data[, vars] = data[, vars]
            } else {
                lags = abs(ccm_result$tar.lag[idx_var])
                vars = as.character(ccm_result$target[idx_var])
                item[[name_Smap_model[i]]]$data[, vars] = c(rep(NA,lags), data[, vars][1:(num_sample-lags)])
            }
        }
    }
    
    return(item)
}

# function to perform S-map analysis
performSmap = function(data_for_smap){
    data = data.frame(data_for_smap)
    
    # determine the optimal theta (non-linearity parameter)
    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 
              0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
    rho.theta = data.frame(matrix(0, length(theta), 2))
    colnames(rho.theta) = c('theta', 'rho')
    rho.theta[, 'theta'] = theta
    
    # compute rho for each S-map model with different theta
    n_theta = length(theta)
    for (i in 1:n_theta){
        block_lnlp_output <- block_lnlp(
            data, columns = c(1:dim(data)[2]), target_column = 1, tp = 1,
            method = "s-map", num_neighbors = 0, stats_only = F,
            theta = rho.theta[i, 'theta'], save_smap_coefficients = T
        )
        rho.theta[i, "rho"] = block_lnlp_output$rho
    }
    theta.opt = rho.theta$theta[which.max(rho.theta$rho)]
    
    block_lnlp_output <- block_lnlp(
        data, columns = c(1:dim(data)[2]), target_column = 1, tp = 1, 
        method = "s-map", num_neighbors = 0, stats_only = F, 
        theta = theta.opt, save_smap_coefficients = T
    )
    
    # The first E columns are for the E lags or causal variables,
    # while the (E+1)th column is the constant
    coeff = data.frame(block_lnlp_output$smap_coefficients[[1]])  # s-map coefficient
    colnames(coeff) = c(colnames(data), "Constant")
    
    rho = block_lnlp_output$rho
    
    # test on the significance of rho
    n_pred = block_lnlp_output$num_pred
    t = rho*sqrt(n_pred-2)/sqrt(1-rho^2)
    if (t >= 0){
        pvalue = 1 - pt(t, df = n_pred-2)    
    } else {
        pvalue = pt(t, df = n_pred-2)
    }
    
    theta = round(theta.opt, 2)
    
    return(list(coefficients = coeff, rho = rho, pvalue = pvalue, theta = theta))
}


## plot S-map coefficients
# give each variable a unique color and shape
# variable: c('F', 'spatial CV', 'age diversity','abundance',
#             'AMO', 'SBT'(SST), 'CV of BT'(CV of SST))
# color: black, yellow, red, green, blue, light blue, purple : c(1,7,2,3,4,5,6)
# cl = c(1,7,2,3,4,5,6,5,6)
# shape: crossbox, cross, circle(s), circle, diamond(s), triangle(s), triangle: c(7,4,21,1,23,24,2)
variables = c("F", "CV.CPUE", "AgeDiversity", "Abundance", "AMO", 
              "SBT", "CVofSBT", "SST", "CVofSST")
cl = plasma(7, direction = -1)
cl = c(cl, cl[6:7])
sh = c(7,4,21,1,23,24,2,24,2)
