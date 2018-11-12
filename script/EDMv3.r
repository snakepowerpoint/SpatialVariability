wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "output"))

library(Rcpp)
library(devtools)
library(rEDM)

species = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus",
            "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
            "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")

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
    var.name = c("Year", "Quarter", "CV.CPUE", "AgeDiversity", "Abundance", "AMO", 
                 "SBT", "CVofSBT",
                 "SST", "CVofSST")
    names(item$data) = var.name
    names(item$data_std) = var.name
    
    return(item)
})


# drop temperature
col_to_drop = c("SST", "CVofSST")
#col_to_drop = c("SBT", "CVofSBT")
EDM_list = lapply(EDM_list, FUN = function(item, col_to_drop){
    df = item$data_std
    item$data_std = df[, -which(names(df) %in% col_to_drop)]
    
    return(item)
}, col_to_drop = col_to_drop)

# drop discontinuous data
EDM_list$`Pleuronectes platessa`$data_std = EDM_list$`Pleuronectes platessa`$data_std[19:48, ]



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
    item$E = computeEmbedding(item$data_std[, -c(1,2)], show_result = FALSE)
    
    return(item)
})


### a function using CCM to determine causal variables (including lagged terms) of spatial CV
# Spatial CV ~ Shannon age, Total CPUE, AMO, Mean BT, CV of BT
library(data.table)

determineCausality = function(data, dim.list, species,
                              lags = 8, num_samples = 100, show_result = FALSE){
    cols = c('length', 'lib.size', 'E', 'library', 'target', 'tar.lag', 
             'rho', 'sd.rho', 'kendall.tau', 'significance')
    cv.var = data.frame(matrix(0, length(lags:0), length(cols)))  # ccm results (cv xmap var)
    colnames(cv.var) = cols
    
    output = data.frame(matrix(0, ncol = length(cols), nrow = 0))  # ccm results for all variables
    colnames(output) = colnames(cv.var)
    
    for (idx in 2:length(data)){
        for (lag in -lags:0){
            cv_x <- ccm(data, E = dim.list[[idx]]$peaks[1], tau = 1, tp = lag, 
                        lib_sizes = seq(5, dim(data)[1], 3),
                        lib_column = 1, target_column = idx, num_samples = num_samples)
            cv_x_m <- data.frame(ccm_means(cv_x), sd.rho = with(cv_x, tapply(rho, lib_size, sd)))
            
            n = dim(cv_x_m)[1]  # last row: ccm results with maximal library size
            
            cv.var[(lag+9), "length"] = dim(data)[1]
            cv.var[(lag+9), "lib.size"] = cv_x_m$lib_size[n]
            cv.var[(lag+9), "E"] = dim.list[[idx]]$peaks[1]
            cv.var[(lag+9), "library"] = "cv"
            cv.var[(lag+9), "target"] = names(dim.list)[idx]
            cv.var[(lag+9), "tar.lag"] = lag
            cv.var[(lag+9), "rho"] = cv_x_m$rho[n]
            cv.var[(lag+9), "sd.rho"] = cv_x_m$sd.rho[n]
            cv.var[(lag+9), "kendall.tau"] = cor.test(cv_x_m$lib_size, cv_x_m$rho, method = 'kendall')$p.value
            cv.var[(lag+9), "significance"] = t.test(subset(cv_x, subset = lib_size==cv_x_m$lib_size[n])$rho,
                                                     alternative = "greater")$p.value
        } 
        
        if (show_result == TRUE){
            # plot ccm results
            # cv xmap variables
            plot(x = -lags:0, y = cv.var$rho, type = 'l', col = 'blue', ylim = c(-0.5,1), xaxt = 'n',
                 xlab = expression(paste('Cross map lag (', italic('l'), ')')), main = species, 
                 ylab = expression(paste('Correlation coefficient ( ', rho, ' )')))
            axis(1, at = seq(-8, 0, 1))
            segments(-lags:0, cv.var[, 'rho'] - cv.var[, 'sd.rho'],
                     -lags:0, cv.var[, 'rho'] + cv.var[, 'sd.rho'], col = 'blue')
            segments(-lags:0 - 0.2, cv.var[, 'rho'] - cv.var[, 'sd.rho'],
                     -lags:0 + 0.2, cv.var[, 'rho'] - cv.var[, 'sd.rho'], col = 'blue')
            segments(-lags:0 - 0.2, cv.var[, 'rho'] + cv.var[, 'sd.rho'],
                     -lags:0 + 0.2, cv.var[, 'rho'] + cv.var[, 'sd.rho'], col = 'blue')
            abline(h = 0)
            legend(x = -3, y = 0.98, legend = paste0('cv xmap ', names(dim.list)[idx]), 
                   text.col = c('blue'))
        }

        output = rbind(output, cv.var)
    }
    
    #output$target = factor(output$target, levels=unique(output$target))
    #output = split(output, output$target)
    
    return(output)  #combine embedding dimension of spatial CV with results
}

# run CCM 100 times to find the best lag for each variable
EDM_list = lapply(EDM_list, function(item){
    item$ccm = data.frame(matrix(0, nrow = 0, ncol = 10))
    for (i in 1:100){
        ccm_result = determineCausality(data = item$data_std[, -c(1,2)], 
                                        dim.list = item$E, 
                                        species = item$species,
                                        show_result = FALSE)
        item$ccm = rbind(item$ccm, ccm_result)
    }
    
    return(item)
})

# save ccm results
lapply(EDM_list, function(item, which_temp="BT"){
    write.csv(item$ccm, 
              file = paste0(wd, "output\\ccm\\ccm_", item$species, "_", which_temp, ".csv"),
              row.names = FALSE)
})

# load ccm results
EDM_list = lapply(EDM_list, function(item, which_temp="BT"){
    item$ccm = read.csv(paste0(wd, "output\\ccm\\ccm_", item$species, "_", which_temp, ".csv"),
                        header = TRUE)
        
    return(item)
})

# keep significant lagged terms
EDM_list = lapply(EDM_list, function(item, min_count = 100*0.95){
    # keep lags which pass the kendall's tau test and t-test, and sort lags by CCM rho
    ccm_sig_lag = subset(item$ccm, subset = item$ccm$kendall.tau < 0.05 & item$ccm$significance < 0.05)
    ccm_sig_lag_count = with(ccm_sig_lag, 
                             aggregate(tar.lag, 
                                       by = list(target = target, tar.lag = tar.lag), 
                                       length))
    ccm_sig_lag_meah_rho = with(ccm_sig_lag, 
                                aggregate(rho, 
                                          by = list(target = target, tar.lag = tar.lag), 
                                          mean))
    ccm_sig_lag_count = cbind(ccm_sig_lag_count, ccm_sig_lag_meah_rho$x)
    names(ccm_sig_lag_count) = c("target", "tar.lag", "count", "rho")
    
    ccm_sig_lag_count = ccm_sig_lag_count[ccm_sig_lag_count$count >= min_count, ]
    ccm_sig_lag_count = ccm_sig_lag_count[order(ccm_sig_lag_count$target,
                                                ccm_sig_lag_count$rho,
                                                decreasing = TRUE), ]
    
    # find the one having maximal rho
    ccm_sig_lag_count = split(ccm_sig_lag_count, f = ccm_sig_lag_count$target, drop = TRUE)
    ccm_sig_lag_count = lapply(ccm_sig_lag_count, function(x) x[1,])
    ccm_sig_lag_count = do.call(rbind, ccm_sig_lag_count)
    item$ccm_most_significant = ccm_sig_lag_count[order(ccm_sig_lag_count$rho, 
                                                        decreasing = TRUE), ]
    
    return(item)
})

# show all significant lagged variables
lapply(EDM_list, function(item){
    
    return(item$ccm_most_significant)
})


### Check if the local-optimal embedding dimension of spatial CV is less than the number of
### causal variables. If so, we can use S-map to reconstruct the state space
feasibleCVEmbedding = function(item){
    embeddings = item$E$CV.CPUE
    embeddings$feasible = (embeddings$peaks - 1) <= length(item$ccm_most_significant$target)
    
    item$E_feasible = embeddings[embeddings$feasible==TRUE, ]
    
    return(item)
}

EDM_list = lapply(EDM_list, feasibleCVEmbedding)
# check if each species has feasible embeddings
lapply(EDM_list, function(item){return(item$E_feasible)})  #Sprattus sprattus

### generate a data frame corresponding to lagged variables for S-map analysis 
generateSmapData = function(item){
    num_Smap_model = length(item$E_feasible$peaks)
    if (num_Smap_model < 1){
        item$Smap1 = NA
        
        return(item)
    }
    name_Smap_model = paste0("Smap", 1:num_Smap_model)
    
    num_sample = length(item$data_std$CV.CPUE)
    num_variable = item$E_feasible$peaks
    
    ccm_result = item$ccm_most_significant
    ccm_result = ccm_result[order(ccm_result$rho, decreasing = TRUE), ]
        
    for (i in 1:num_Smap_model){
        item[[name_Smap_model[i]]] = list()
        item[[name_Smap_model[i]]]$data = data.frame(matrix(0, nrow = num_sample, ncol = 0))
        item[[name_Smap_model[i]]]$data[, "CV.CPUE"] = item$data_std[, "CV.CPUE"]
        
        for (idx_var in 1:(num_variable[i]-1)){
            lags = abs(ccm_result$tar.lag[idx_var])
            vars = as.character(ccm_result$target[idx_var])
            item[[name_Smap_model[i]]]$data[, vars] = c(rep(NA,lags), item$data_std[, vars][1:(num_sample-lags)])
        }
    }
    
    return(item)
}

EDM_list = lapply(EDM_list, generateSmapData)


### perform S-map
performSmap = function(data_for_smap){
    data = data.frame(data_for_smap)
    
    # determine the optimal theta
    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 
              0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
    rho.theta = data.frame(matrix(0, length(theta), 2))
    colnames(rho.theta) = c('theta', 'rho')
    rho.theta[, 'theta'] = theta
    
    # compute rho under each theta
    for (i in 1:length(theta)){
        block_lnlp_output <- block_lnlp(data, columns = c(1:dim(data)[2]), target_column = 1, tp = 1,
                                        method = "s-map", num_neighbors = 0, stats_only = F, theta = rho.theta[i, 'theta'],
                                        save_smap_coefficients = T)
        rho.theta[i, "rho"] = block_lnlp_output$rho
    }
    theta.opt = rho.theta$theta[which.max(rho.theta$rho)]
    
    block_lnlp_output <- block_lnlp(data, columns = c(1:dim(data)[2]), target_column = 1, tp = 1, 
                                    method = "s-map", num_neighbors = 0, stats_only = F, theta = theta.opt,
                                    save_smap_coefficients = T)
    
    # The first E columns are for the E lags or causal variables,
    # while the (E+1)th column is the constant
    # we plot only the causal variables
    coeff = data.frame(block_lnlp_output$smap_coefficients[[1]])  # s-map coefficient
    colnames(coeff) = c(colnames(data), "Constant")
    
    return(list(coefficients = coeff, rho = round(block_lnlp_output$rho, 2),
                theta = round(theta.opt, 2)))
}

EDM_list = lapply(EDM_list, function(item){
    if (is.na(item$Smap1)){
        return(item)
    }
    
    smap_model_index = grep("Smap", names(item))
    for (i in 1:length(smap_model_index)){
        data = item[[smap_model_index[i]]]$data
        smap_result_list = performSmap(data)
        
        item[[smap_model_index[i]]]$coefficients = smap_result_list$coefficients
        item[[smap_model_index[i]]]$rho = smap_result_list$rho
        item[[smap_model_index[i]]]$theta = smap_result_list$theta
    }
    
    return(item)
})


### plot S-map coefficients
library(ggplot2)

plotSmapCoeff = function(data_of_coeff, lag_of_var, species=NULL, 
                         mode = "series", rho){
    # data_of_coeff contains coefficients of all variables excluding CV.CPUE and constant
    #sort data
    variables = c("AgeDiversity", "Abundance", "AMO", 
                  "SBT", "CVofSBT", "SST", "CVofSST")
    order.var = variables[sort(match(names(data_of_coeff), variables))]
    lag_of_var = lag_of_var[order(match(names(data_of_coeff), order.var))]
    data_of_coeff = data_of_coeff[, order.var, drop = FALSE]
    
    ntime = dim(data_of_coeff)[1]
    nvar = dim(data_of_coeff)[2]
    coeff.melt = cbind(date = rep(1:ntime, nvar), melt(data_of_coeff))
    
    #give each variable an unique color and shape
    #variable: c('age diversity','abundance','AMO','SBT'(SST),'CV of BT'(CV of SST))
    #color: red, green, blue, light blue, purple : c(2,3,4,5,6)
    #shape: circle(s),circle,diamond(s),triangle(s),triangle: c(21,1,23,24,2)
    whichvar = match(names(data_of_coeff), variables)
    cl = rep(c(2,3,4,5,6,5,6)[whichvar], each = ntime) 
    sh = rep(c(21,1,23,24,2,24,2)[whichvar], each = ntime)
    
    if (!mode %in% c("series", "boxplot")){
        stop("mode must be either 'series' or 'boxplot'")
    }
    
    scaleFUN = function(x){sprintf("%.2f", x)}
    if (mode == "series"){ 
        #plot time series
        smaptime = ggplot(data = coeff.melt, aes(x = date, y = value, group = variable))
        print(smaptime + geom_point(size = 2, col = cl, shape = sh, fill = cl) +
                  geom_line(col = cl) +
                  geom_hline(yintercept = 0, lty = 2) +
                  theme(plot.title = element_text(hjust = 0.5, size = 20),
                        axis.title = element_text(size = 18, face = "bold"),
                        axis.text = element_text(size = 18, colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
                        panel.background = element_blank()) + 
                  xlab('Time') +
                  ylab('S-map coefficients') +
                  scale_y_continuous(labels=scaleFUN) +
                  ggtitle(bquote(paste(italic(.(species)), ", ", rho, " = ", .(rho)))))
    }
    
    if (mode == "boxplot"){
        #plot box plot
        smapbox = ggplot(data = coeff.melt, aes(x = variable, y = value, group = variable))
        print(smapbox + geom_boxplot(na.rm = T, col = c(2,3,4,5,6,5,6)[whichvar], lwd = 1, width = 0.5*nvar/5) + 
                  geom_hline(yintercept = 0, lty = 2) +
                  theme(axis.title.x = element_blank(), 
                        plot.title = element_text(hjust = 0.5, size = 20),
                        axis.title = element_text(size = 18, face = "bold"),
                        axis.text = element_text(size = 18, colour = "black"),
                        panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
                        panel.background = element_blank()) + 
                  ylab('S-map coefficients') +
                  scale_y_continuous(labels=scaleFUN) +
                  scale_x_discrete(labels=paste0(colnames(data_of_coeff), "-", lag_of_var)) + 
                  ggtitle(bquote(paste(italic(.(species)), ", ", rho, " = ", .(rho)))))
    }
    
    return(NULL)
}

lapply(EDM_list, function(item, which_temp="BT", mode="series"){
    if (is.na(item$Smap1)){
        return(NULL)
    }
    
    smap_model_index = grep("Smap", names(item))
    for (i in 1:length(smap_model_index)){
        data_of_coeff = item[[smap_model_index[i]]]$coefficients
        data_of_coeff = data_of_coeff[, -which(names(data_of_coeff) %in% c("CV.CPUE", "Constant")), drop = FALSE]
        
        data_for_smap = item[[smap_model_index[i]]]$data
        data_for_smap = data_for_smap[, !names(data_for_smap) %in%  c("CV.CPUE"), drop = FALSE]
        
        lag_of_var = as.numeric(colSums(is.na(data_for_smap)))
        species = item$species
        rho = item[[smap_model_index[i]]]$rho
        
        plotSmapCoeff(data_of_coeff = data_of_coeff, lag_of_var = lag_of_var,
                      species = species, mode = mode, rho = rho)
        ggsave(paste0(wd, "output\\smap\\smap_", mode, species, i, "_", which_temp, ".png"),
               width = 9, height = 6, units = "in")
        ggsave(paste0(wd, "output\\smap\\smap_", mode, species, i, "_", which_temp, ".eps"))
    }
    
    return(NULL)
})

# mean of S-map coefficients over time
smap_results_list = lapply(EDM_list, function(item){
    smap_model_index = grep("Smap", names(item))
    
    meanSmap = list()
    for (idx in 1:length(smap_model_index)){
        if (length(item[[smap_model_index[idx]]]) > 1){
            meanSmap[[idx]] = colMeans(item[[smap_model_index[idx]]]$coefficients, na.rm = TRUE)
            meanSmap[[idx]] = data.frame(t(meanSmap[[idx]]))
            meanSmap[[idx]]$rho = item[[smap_model_index[idx]]]$rho
            meanSmap[[idx]]$theta = item[[smap_model_index[idx]]]$theta
        }
    }
    
    return(meanSmap)
})

# dimension of the system
lapply(EDM_list, function(item){
    
    return(item$E_feasible)
})

# store S-map coefficients as data frame
library(gtools)

smap_results_df = data.frame()
for (species in smap_results_list){
    for (list in species){
        data = data.frame(as.list(list))
        smap_results_df = smartbind(smap_results_df, data)
    }
}

variables = c("AgeDiversity", "Abundance", "AMO", 
              "SBT", "CVofSBT", "SST", "CVofSST", "theta", "rho")
order.var = variables[sort(match(names(smap_results_df), variables))]
smap_results_df = smap_results_df[, order.var, drop = FALSE]
smap_results_df = as.data.frame(apply(smap_results_df, 2, round, digits = 4))

model_names = c()
for (idx in 1:length(smap_results_list)){
    if (length(smap_results_list[[idx]]) > 0){
        num_model = length(smap_results_list[[idx]])
        model_name = paste0(names(smap_results_list)[idx], 1:num_model)
        model_names = c(model_names, model_name)
    }
}

smap_results_df$model_name = model_names
row.names(smap_results_df) = NULL

#write.csv(x = smap_results_df, 
#          file = paste0(wd, "output\\smap\\mean_coefficients_", "BT", ".csv"))



### robustness test on S-map
source(paste0(wd, "script\\robustness_test_function.r"))

robustness_list = lapply(EDM_list, FUN = function(species){
    data = species$data_std[, -c(1,2)] # exclude colunm "Year" and "Quarter"
    
    if (nrow(species$E_feasible) > 0){
        dim.cv = species$E_feasible$peaks[1]
        vars = species$ccm_most_significant$target
        vars = as.character(vars)
        lags = species$ccm_most_significant$tar.lag
        lags = abs(lags)
        
        data = subset(data, select = c("CV.CPUE", vars))
        
        robust(data, dim.cv = dim.cv, lags = lags)
    }
})

robustness_list = lapply(robustness_list, FUN=function(data){
    if (is.data.frame(data)){
        variables = c("AgeDiversity", "Abundance", "AMO", 
                      "SBT", "CVofSBT", "SST", "CVofSST", "theta", "rho")
        order.var = variables[sort(match(names(data), variables))]
        data = data[, order.var, drop = FALSE]
        
        return(round(data, digits = 4))
    }
    
    else {return(NULL)}
})

for (idx in 1:length(robustness_list)){
    write.csv(x = robustness_list[[idx]], 
              file = paste0(wd, "output\\robustness_test\\", names(robustness_list)[idx],
                            "_BT", ".csv"),
              row.names = FALSE)
}




### Appendix
john = EDM_list$`Melanogrammus aeglefinus`$data_std
robust(data = john[, -c(1,2)], dim.cv = 4, lags = c(7,6,7,4,4))

#species name
"
Clupea harengus
Gadus morhua
Melanogrammus aeglefinus
Merlangius merlangus
Pleuronectes platessa
Pollachius virens
Scomber scombrus
Sprattus sprattus
Trisopterus esmarkii
"
