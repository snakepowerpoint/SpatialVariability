### CCM analysis to determine causality
# Spatial CV ~ Age diversity, Abundance, AMO, SBT/SST, CV of SBT/SST

########## Warning!
# This step takes much time. 
# One can skip this step if the pre-run CCM results are provided.
# Please see the loading code below.
########## Warning!
# run CCM several times to find the best lag for each variable
EDM_lib_var = lapply(EDM_lib_var, function(item, lib_var=library_var, lag=lags, t_ccm=time_ccm){
    data = item[[dataset]]
    data = subset(data, select=names(data) %ni% c("Year", "Quarter"))
    data = cbind(data[lib_var], subset(data, select=names(data) %ni% lib_var))
    
    item$ccm = data.frame(matrix(0, nrow = 0, ncol = 10))
    for (i in 1:t_ccm){
        seed = 1234 + i * 10
        ccm_result = determineCausality(data = data, 
                                        dim.list = item$E, 
                                        species = item$species,
                                        lags = lag,
                                        seed = seed)
        item$ccm = rbind(item$ccm, ccm_result)
    }
    
    return(item)
})

# save ccm results
lapply(EDM_lib_var, function(item){
    write.csv(item$ccm, 
              file = paste0(wd, ccm_path, item$species, ".csv"),
              row.names = FALSE)
})
