### Data preparation
wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "script"))
source("EDM.r")


### Simplex projection to determine embedding dimension
EDM_age = lapply(EDM_list, FUN = function(item){
    data = subset(item$data_std, select=-c(Year, Quarter))
    item$E = computeEmbedding(data, show_result = FALSE)
    
    return(item)
})


### CCM analysis to determine causality
# Shannon age ~ Spatial CV, Total CPUE, AMO, Mean BT, CV of BT
# run CCM 100 times to find the best lag for each variable
EDM_age = lapply(EDM_age, function(item, name_var='AgeDiversity'){
    data = item$data_std[, -c(1, 2)]
    `%ni%` <- Negate(`%in%`)
    data = cbind(data[name_var], subset(data, select=names(data) %ni% name_var))
    
    item$ccm = data.frame(matrix(0, nrow = 0, ncol = 10))
    for (i in 1:2){
        ccm_result = determineCausality(data = data, 
                                        dim.list = item$E, 
                                        species = item$species)
        item$ccm = rbind(item$ccm, ccm_result)
    }
    
    return(item)
})

# save ccm results
subpath = "output\\ccm\\age_diversity\\"
dir.create(file.path(wd, subpath), showWarnings = FALSE)

lapply(EDM_age, function(item){
    write.csv(item$ccm, 
              file = paste0(wd, subpath, item$species, ".csv"),
              row.names = FALSE)
})

# load ccm results
EDM_age = lapply(EDM_age, function(item){
    item$ccm = read.csv(paste0(wd, subpath, item$species, ".csv"), header = TRUE)
    return(item)
})


## keep significant lagged terms
EDM_age = lapply(EDM_age, function(item, min_count = 2*0.95){
    # keep lags which pass the kendall's tau test and t-test, and sort lags by its CCM rho
    ccm_sig_lag = subset(item$ccm, subset = item$ccm$kendall.tau < 0.05 & item$ccm$significance < 0.05)
    ccm_sig_lag_count = with(
        ccm_sig_lag, 
        aggregate(tar.lag, by = list(target = target, tar.lag = tar.lag), length))
    ccm_sig_lag_mean_rho = with(
        ccm_sig_lag, 
        aggregate(rho, by = list(target = target, tar.lag = tar.lag), mean))
    ccm_sig_lag_count = cbind(ccm_sig_lag_count, ccm_sig_lag_mean_rho$x)
    names(ccm_sig_lag_count) = c("target", "tar.lag", "count", "rho")
    
    ccm_sig_lag_count = ccm_sig_lag_count[ccm_sig_lag_count$count >= min_count, ]
    ccm_sig_lag_count = ccm_sig_lag_count[
        order(ccm_sig_lag_count$target, ccm_sig_lag_count$rho, decreasing = TRUE), ]
    
    # find the one having maximal rho
    ccm_sig_lag_count = split(ccm_sig_lag_count, f = ccm_sig_lag_count$target, drop = TRUE)
    ccm_sig_lag_count = lapply(ccm_sig_lag_count, function(x) x[1,])
    ccm_sig_lag_count = do.call(rbind, ccm_sig_lag_count)
    item$ccm_most_significant = ccm_sig_lag_count[order(ccm_sig_lag_count$rho, decreasing = TRUE), ]
    
    return(item)
})

# show all significant lagged variables
lapply(EDM_age, function(item){
    return(item$ccm_most_significant)
})


## Check if the local-optimal embedding dimension of spatial CV is less than the number of
## causal variables. If so, we can use S-map to reconstruct the state space
feasibleCVEmbedding = function(item, name_var="AgeDiversity"){
    embeddings = item$E[[name_var]]
    embeddings$feasible = (embeddings$peaks - 1) <= length(item$ccm_most_significant$target)
    
    item$E_feasible = embeddings[embeddings$feasible==TRUE, ]
    
    return(item)
}

EDM_age = lapply(EDM_age, feasibleCVEmbedding)

# check if each species has feasible embeddings
lapply(EDM_age, function(item){return(item$E_feasible)})  


## generate a data frame corresponding to lagged variables for S-map analysis 
generateSmapData = function(item, data, name_var){
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
        item[[name_Smap_model[i]]]$data[, name_var] = data[name_var]
        
        for (idx_var in 1:(num_variable[i]-1)){
            lags = abs(ccm_result$tar.lag[idx_var])
            vars = as.character(ccm_result$target[idx_var])
            item[[name_Smap_model[i]]]$data[, vars] = c(rep(NA,lags), data[, vars][1:(num_sample-lags)])
        }
    }
    
    return(item)
}

EDM_age = lapply(EDM_age, FUN=function(item){
    name_var = "AgeDiversity"
    data = item$data_std
    generateSmapData(item, data=data, name_var=name_var)
})



### S-map
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
    
    rho = round(block_lnlp_output$rho, 2)
    
    # test on the significance of rho
    n_pred = block_lnlp_output$num_pred
    t = rho*sqrt(n_pred-2)/sqrt(1-rho^2)  
    pvalue = 2*pt(-abs(t), df = n_pred-2)
    
    theta = round(theta.opt, 2)
    
    return(list(coefficients = coeff, rho = rho, pvalue = pvalue, theta = theta))
}

EDM_age = lapply(EDM_age, function(item){
    if (is.na(item$Smap1)){
        return(item)
    }
    
    smap_model_index = grep("Smap", names(item))
    num_smap_model = length(smap_model_index)
    for (i in 1:num_smap_model){
        data = item[[smap_model_index[i]]]$data
        smap_result_list = performSmap(data)
        
        item[[smap_model_index[i]]]$coefficients = smap_result_list$coefficients
        item[[smap_model_index[i]]]$rho = smap_result_list$rho
        item[[smap_model_index[i]]]$pvalue = smap_result_list$pvalue
        item[[smap_model_index[i]]]$theta = smap_result_list$theta
    }
    
    return(item)
})


### plot S-map coefficients
setwd(paste0(wd, 'script'))
source("utils/plot.r")






lapply(EDM_age, function(item, which_temp="BT", mode="series"){
    if (is.na(item$Smap1)){
        return(NULL)
    }
    
    species = item$species
    smap_model_index = grep("Smap", names(item))
    for (i in 1:length(smap_model_index)){
        # extract lags for each variables
        data_for_smap = item[[smap_model_index[i]]]$data
        target_vars = names(data_for_smap)[-1]
        num_na = colSums(is.na(data_for_smap))
        lag_of_var = as.numeric(num_na[target_vars])
        
        # coefficients of S-map model without library variables and constant
        data_of_coeff = item[[smap_model_index[i]]]$coefficients
        data_of_coeff = data_of_coeff[target_vars] 
        
        rho = item[[smap_model_index[i]]]$rho
        
        plotSmapCoeff(data_of_coeff = data_of_coeff, lag_of_var = lag_of_var,
                      species = species, mode = mode, rho = rho)
        file_name_eps = paste0(wd, "output\\smap\\smap_", mode, species, i, "_", which_temp, ".png")
        ggsave(filename = file_name_eps, width = 9, height = 6, units = "in")
        file_name_png = paste0(wd, "output\\smap\\smap_", mode, species, i, "_", which_temp, ".eps")
        ggsave(filename = file_name_png)
    }
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
            meanSmap[[idx]]$pvalue = item[[smap_model_index[idx]]]$pvalue
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
              "SBT", "CVofSBT", "SST", "CVofSST", "theta", "rho", "pvalue")
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
                      "SBT", "CVofSBT", "SST", "CVofSST", "theta", "rho", "pvalue")
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
