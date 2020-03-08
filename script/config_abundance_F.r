## config 
dataset = "data_per_year"
lags = 4
is_robust_each_lag = TRUE

time_ccm = 200
detrend_fun = detrend_sig
is_full = FALSE

library_var = "Aundance"


## path for the results of ccm 
ccm_path = paste0("output\\ccm\\", library_var, "_sig_detrend_F\\")
dir.create(file.path(wd, ccm_path), showWarnings = FALSE)


## path for the results of smap
smap_path = paste0("output\\smap\\", library_var, "_sig_detrend_F\\")
dir.create(file.path(wd, smap_path), showWarnings = FALSE)


## path for the results of robustness tes
robust_path = paste0("output\\robustness_test\\", library_var, "_sig_detrend_F\\")
dir.create(file.path(wd, robust_path), showWarnings = FALSE)


## plot S-map coefficients
# give each variable a unique color and shape
# variable: c('F', 'spatial CV', 'age diversity','abundance',
#             'AMO', 'SBT'(SST), 'CV of BT'(CV of SST))
# color: black, yellow, red, green, blue, light blue, purple : c(1,7,2,3,4,5,6)
# shape: crossbox, cross, circle(s), circle, diamond(s), triangle(s), triangle: c(7,4,21,1,23,24,2)
variables = c("F", "CV.CPUE", "AgeDiversity", "Abundance", "AMO", 
              "SBT", "CVofSBT", "SST", "CVofSST")
cl = c(1,7,2,3,4,5,6,5,6)
sh = c(7,4,21,1,23,24,2,24,2)