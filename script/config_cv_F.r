## config 
dataset = "data_per_year" 
lags = 4
is_robust_each_lag = FALSE 

time_ccm = 200
detrend_fun = detrend_sig
is_full = FALSE

library_var = "CV.CPUE"


## path for the results of ccm 
ccm_path = paste0("output\\ccm\\", library_var, "_sig_detrend_F\\")
dir.create(file.path(wd, ccm_path), showWarnings = FALSE)


## path for the results of smap
smap_path = paste0("output\\smap\\", library_var, "_sig_detrend_F\\")
dir.create(file.path(wd, smap_path), showWarnings = FALSE)


## path for the results of robustness tes
robust_path = paste0("output\\robustness_test\\", library_var, "_sig_detrend_F\\")
dir.create(file.path(wd, robust_path), showWarnings = FALSE)

