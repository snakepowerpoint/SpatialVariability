## config 
dataset = "data_std"  # swtich to "data_per_year" if yearly data is used
lags = 8  # swtich to "4" if time series is short (when fishing mortality is included)
is_robust_each_lag = FALSE  # swtich to "TRUE" if time series is short

time_ccm = 200
detrend_fun = detrend_sig
is_full = FALSE

library_var = "CV.CPUE"


## path for the results of ccm 
ccm_path = paste0("output\\ccm\\", library_var, "_sig_detrend\\")
dir.create(file.path(wd, ccm_path), showWarnings = FALSE)


## path for the results of smap
smap_path = paste0("output\\smap\\", library_var, "_sig_detrend\\")
dir.create(file.path(wd, smap_path), showWarnings = FALSE)


## path for the results of robustness tes
robust_path = paste0("output\\robustness_test\\", library_var, "_sig_detrend\\")
dir.create(file.path(wd, robust_path), showWarnings = FALSE)
