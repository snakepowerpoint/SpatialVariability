### Data preparation
wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "script"))
source("edm.r")


# detrend
detrend_fun = detrend_sig
EDM_list = lapply(EDM_list, FUN=function(item){
    data_std = item$data_std[-c(1,2)]
    data_per_year = item$data_per_year[-c(1)]
    
    # apply linear regression
    data_std = apply(data_std, 2, FUN=detrend_fun)
    data_per_year = apply(data_per_year, 2, FUN=detrend_fun)
    
    item$data_std = data.frame(cbind(item$data_std[c(1,2)], data_std))
    item$data_per_year = data.frame(cbind(item$data_per_year[c(1)], data_per_year))
    
    return(item)
})



### Simplex projection to determine embedding dimension
`%ni%` <- Negate(`%in%`)

dataset = "data_std"
EDM_lib_var = lapply(EDM_list, FUN = function(item){
    data = item[[dataset]]
    data = subset(data, select=names(data) %ni% c("Year", "Quarter"))
    item$E = computeEmbedding(data, show_result = FALSE)
    
    return(item)
})


### CCM analysis to determine causality
# Spatial CV ~ Age diversity, Abundance, AMO, SBT/SST, CV of SBT/SST
lags = 8
time_ccm = 100
library_var = "CV.CPUE"

########## Warning!
# This step takes much time. 
# One can skip this step if the pre-run CCM results are provided.
# Please see the loading code below.
########## Warning!
# run CCM 100 times to find the best lag for each variable
EDM_lib_var = lapply(EDM_lib_var, function(item, lib_var=library_var, lag=lags, t_ccm=time_ccm){
    data = item[[dataset]]
    data = subset(data, select=names(data) %ni% c("Year", "Quarter"))
    data = cbind(data[lib_var], subset(data, select=names(data) %ni% lib_var))
    
    item$ccm = data.frame(matrix(0, nrow = 0, ncol = 10))
    for (i in 1:t_ccm){
        ccm_result = determineCausality(data = data, 
                                        dim.list = item$E, 
                                        species = item$species,
                                        lags = lag)
        item$ccm = rbind(item$ccm, ccm_result)
    }
    
    return(item)
})

# save ccm results
subpath = paste0("output\\ccm\\", library_var, "\\")
dir.create(file.path(wd, subpath), showWarnings = FALSE)
lapply(EDM_lib_var, function(item){
    write.csv(item$ccm, 
              file = paste0(wd, subpath, item$species, ".csv"),
              row.names = FALSE)
})

# load ccm results
subpath = paste0("output\\ccm\\", library_var, "\\")
dir.create(file.path(wd, subpath), showWarnings = FALSE)
EDM_lib_var = lapply(EDM_lib_var, function(item){
    item$ccm = read.csv(paste0(wd, subpath, item$species, ".csv"), header = TRUE)
    return(item)
})


## keep significant lagged terms
EDM_lib_var = lapply(EDM_lib_var, function(item, min_count = time_ccm*0.95){
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
    item$ccm_sig = do.call(rbind, ccm_sig_lag_count)
    
    ccm_sig_lag_count = lapply(ccm_sig_lag_count, function(x) x[1,])
    ccm_sig_lag_count = do.call(rbind, ccm_sig_lag_count)
    item$ccm_most_significant = ccm_sig_lag_count[order(ccm_sig_lag_count$rho, decreasing = TRUE), ]
    
    return(item)
})

# keep SBT for demersal species, and SST for pelagic species
drop_temperature = function(item){
    demersal_species = c("Gadus morhua", "Melanogrammus aeglefinus", 
                         "Merlangius merlangus", "Pleuronectes platessa",
                         "Pollachius virens", "Trisopterus esmarkii")
    if (item$species %in% demersal_species){
        data = item$ccm_most_significant
        row_to_drop = grep("SST", data$target)
        item$ccm_most_significant = data[-c(row_to_drop), ]
        item$ccm_sig = item$ccm_sig[-c(grep("SST", item$ccm_sig$target)), ]
    }
    else {
        data = item$ccm_most_significant
        row_to_drop = grep("SBT", data$target)
        item$ccm_most_significant = data[-c(row_to_drop), ]
        item$ccm_sig = item$ccm_sig[-c(grep("SBT", item$ccm_sig$target)), ]
    }
    return(item)
}

EDM_lib_var = lapply(EDM_lib_var, FUN=function(item){drop_temperature(item)})

# show all significant lagged variables
lapply(EDM_lib_var, function(item){
    return(item$ccm_most_significant)
})


## Check if the local-optimal embedding dimension of spatial CV is less than the number of
## causal variables. If so, we can use S-map to reconstruct the state space.
feasibleCVEmbedding = function(item, lib_var=library_var){
    embeddings = item$E[[lib_var]]
    embeddings$feasible = (embeddings$peaks - 1) <= length(item$ccm_most_significant$target)
    
    item$E_feasible = embeddings[embeddings$feasible==TRUE, ]
    
    return(item)
}

EDM_lib_var = lapply(EDM_lib_var, feasibleCVEmbedding)

# check if each species has feasible embeddings
lapply(EDM_lib_var, function(item){return(item$E_feasible)})  


## save summaried ccm results
ccm_data = data.frame(0)
for (i in 1:length(EDM_lib_var)){
    data = EDM_lib_var[[i]]$ccm_most_significant
    if (nrow(data) > 0){
        data$species = names(EDM_lib_var)[i]
        ccm_data = merge(ccm_data, data, all=TRUE, sort=FALSE)
    }
}

ccm_table = xtabs(rho ~ species + target, data=ccm_data, sparse=TRUE)
ccm_table = as.data.frame.matrix(ccm_table)

for (i in 1:length(EDM_lib_var)){
    print(EDM_lib_var[[i]]$E[[library_var]]$peaks)
}

write.csv(x=ccm_table, file=paste0(wd, "\\output\\ccm\\ccm.csv"))



### S-map
# generate a data frame corresponding to lagged variables for S-map analysis 
is_full = FALSE
EDM_lib_var = lapply(EDM_lib_var, FUN=function(item, lib_var=library_var){
    data = item[[dataset]]
    generateSmapData(item, data=data, lib_var=lib_var, is_full=is_full)
})

# perform S-map analysis
EDM_lib_var = lapply(EDM_lib_var, function(item){
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
# give each variable a unique color and shape
# variable: c('F', 'spatial CV', 'age diversity','abundance',
#             'AMO', 'SBT'(SST), 'CV of BT'(CV of SST))
# color: black, yellow, red, green, blue, light blue, purple : c(1,7,2,3,4,5,6)
# shape: crossbox, cross, circle(s), circle, diamond(s), triangle(s), triangle: c(7,4,21,1,23,24,2)
variables = c("F", "CV.CPUE", "AgeDiversity", "Abundance", "AMO", 
              "SBT", "CVofSBT", "SST", "CVofSST")
cl = c(1,7,2,3,4,5,6,5,6)
sh = c(7,4,21,1,23,24,2,24,2)
plot_mode = "series"  # 'series' or 'box'

# directory for saving S-map results
save_dir = paste0(wd, "output\\smap\\", library_var, "\\")
dir.create(save_dir, showWarnings = FALSE)

setwd(paste0(wd, 'script'))
source("utils/plot.r")

lapply(EDM_lib_var, function(item, colors=cl, shapes=sh, mode=plot_mode){
    if (is.na(item$Smap1)[1]){
        return(NULL)
    }
    
    species = item$species
    smap_model = grep("Smap", names(item))
    num_smap_model = length(smap_model)
    for (i in 1:num_smap_model){
        smapplot = plotSmapCoeff(smap_result_list=item[[smap_model[i]]],
                                 species=species,
                                 colors=colors,
                                 shapes=shapes,
                                 mode=mode)
        save_path = paste0(save_dir, mode, "_", species, i)
        
        file_name_eps = paste0(save_path, ".eps")
        ggsave(filename=file_name_eps, plot=smapplot, width=9, height=6, units="in")
        file_name_png = paste0(save_path, ".png")
        ggsave(filename=file_name_png, plot=smapplot, width=9, height=6, units="in")
    }
})

if (plot_mode == "series"){
    smap_timeseries_legend(lib_var=library_var, colors=cl, shapes=sh)
} else if (plot_mode == "box"){
    smap_boxplot_legend(lib_var=library_var, colors=cl)
}
file_name_eps = paste0(save_dir, plot_mode, "_legend", ".eps")
ggsave(filename = file_name_eps)



### Save S-map results
# mean of S-map coefficients over time
smap_results_list = lapply(EDM_lib_var, function(item){
    smap_model = grep("Smap", names(item))
    num_smap_model = length(smap_model)
    
    meanSmap = list()
    for (idx in 1:num_smap_model){
        model = item[[smap_model[idx]]]
        if (length(model) > 1){
            meanSmap[[idx]] = colMeans(model$coefficients, na.rm = TRUE)
            meanSmap[[idx]] = data.frame(t(meanSmap[[idx]]))
            meanSmap[[idx]]$rho = model$rho
            meanSmap[[idx]]$pvalue = model$pvalue
            meanSmap[[idx]]$theta = model$theta
        }
    }
    
    return(meanSmap)
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

# re-arrange columns
smap_var = c(variables, "theta", "rho", "pvalue")
order.var = smap_var[sort(match(names(smap_results_df), smap_var))]
smap_results_df = smap_results_df[, order.var, drop = FALSE]
smap_results_df = as.data.frame(apply(smap_results_df, 2, round, digits = 4))

# label each model
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

# embedding dimension of target variable
embed_dim = c()
for (species in EDM_lib_var){
    sub_model_names = gsub("[[:digit:]]+", "", x=model_names)
    if (species$species %in% sub_model_names){
        dims = species$E_feasible$peaks
        embed_dim = c(embed_dim, dims)
    }
}
smap_results_df$E = embed_dim

filename = paste0(wd, "output\\smap\\mean_coefficients_", library_var, ".csv")
write.csv(x=smap_results_df, file=filename)



### Robustness test on S-map
setwd(paste0(wd, 'script'))
source("robustness_test.r")

# remove robustness test on sub-optimal E
robustness_list = lapply(EDM_lib_var, 
                         FUN=function(species, lib_var=library_var, use_all_data=is_full){
    # exclude colunm of "Year" and "Quarter"
    data = species[[dataset]]
    data = subset(data, select=names(data) %ni% c("Year", "Quarter"))
    
    if (nrow(species$E_feasible) > 0){
        tar_vars = as.character(species$ccm_most_significant$target)
        tar_vars = as.character(tar_vars)
        if (use_all_data){
            lags = rep(0, length(tar_vars))
        } else {
            lags = species$ccm_most_significant$tar.lag
            lags = abs(lags)
        }
        data = subset(data, select = c(lib_var, tar_vars))
        
        robust_output = data.frame(0)
        for (emb_dim in species$E_feasible$peaks){
            output = robust(data, dim_lib_var=emb_dim, lags=lags)
            robust_output = merge(robust_output, output, all=TRUE)
        }
        return(robust_output)
    } else {
        return(NULL)
    }
    
})

# merge results of robustness test for all species
robustness_table = data.frame(0)
for (i in 1:length(robustness_list)){
    data = robustness_list[[i]]
    if (!is.null(data)){
        data = round(data, digits=4)
        data$species = names(robustness_list)[i]
        robustness_table = merge(robustness_table, data, all=TRUE, sort=FALSE)
    }
}

robust_var = c(variables, "theta", "rho", "pvalue", "species")
order.var = robust_var[sort(match(names(robustness_table), robust_var))]
robustness_table = robustness_table[, order.var, drop = FALSE]

# directory for saving robustness test results
subpath = paste0("output\\robustness_test\\", library_var, "\\")
dir.create(file.path(wd, subpath), showWarnings = FALSE)
filename = paste0(wd, subpath, "all_species.csv")
write.csv(x=robustness_table, file=filename, row.names = FALSE)



### Robustness test on different lags
robust_lag_list = lapply(EDM_lib_var, FUN=function(species){
    if (nrow(species$E_feasible) > 0){
        raw_data = species[[dataset]]
        embed_dim = species$E_feasible$peaks[1]
        lib_var = library_var
        ccm_most_sig = species$ccm_most_significant
        ccm_sig = species$ccm_sig
        
        print(species$species)
        return(smap_each_lag(raw_data=raw_data, 
                             embed_dim=embed_dim, 
                             lib_var=lib_var, 
                             ccm_most_sig=ccm_most_sig, 
                             ccm_sig=ccm_sig))
    } else {
        return(NULL)
    }
})

robust_lag_list = lapply(robust_lag_list, FUN=function(species){
    robust_var = c("lag", "lag_var", variables, "theta", "rho", "pvalue")
    order.var = robust_var[sort(match(names(species), robust_var))]
    sorted_table = species[, order.var, drop = FALSE]
    
    return(sorted_table)
})

# directory for saving robustness test results
subpath = paste0("output\\robustness_test\\", library_var, "_sig_detrend_each_lag\\")
dir.create(file.path(wd, subpath), showWarnings = FALSE)
for (i in 1:length(robust_lag_list)){
    species = names(robust_lag_list)[i]
    filename = paste0(wd, subpath, species, ".csv")
    write.csv(x=robust_lag_list[[i]], file=filename, row.names=FALSE)
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
