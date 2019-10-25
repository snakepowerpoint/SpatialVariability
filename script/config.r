# config 
dataset = "data_per_year"
lags = 4
time_ccm = 100
library_var = "AgeDiversity"
library_var = "Abundance"
library_var = "CV.CPUE"

# load results of ccm
subpath = paste0("output\\ccm\\", library_var, "_sig_detrend\\")
subpath = paste0("output\\ccm\\", library_var, "_sig_detrend_F\\")

# path for saving results of smap 
filename = paste0(wd, "output\\smap\\mean_coefficients_", library_var, "_sig_detrend.csv")
filename = paste0(wd, "output\\smap\\mean_coefficients_", library_var, "_w_fishingM_sig_detrend.csv")

# path for saving results of robustness test
subpath = paste0("output\\robustness_test\\", library_var, "_sig_detrend\\")
subpath = paste0("output\\robustness_test\\", library_var, "_sig_detrend_each_lag\\")
subpath = paste0("output\\robustness_test\\", library_var, "_w_fishingM_sig_detrend\\")
subpath = paste0("output\\robustness_test\\", library_var, "_w_fishingM_sig_detrend_each_lag\\")

