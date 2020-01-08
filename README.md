# Causal effects of age structure, abundance and environment on spatial variability of marine fishes

This is an implementaion of **Causal effects of age structure, abundance and environment on spatial variability of marine fishes.** For readers who want to quickly implement the empirical dynamic modelling (EDM), please refer to **Quickly Getting Started**. For readers who want to implement the EDM from scratch, please refer to **Step by Step Analysis**.

The repository includes:
* Source code

Raw data including *fishery survey data*, *Atlantic Multidecadal Oscillation (AMO)*, *sea surface temperature (SST)*, *sea bottom tempeature (SBT)*, and *fishing mortality* can be downloaded on Zenodo (https://doi.org/10.5281/zenodo.3518702) or at the links provided in the paper.

# Requirements 
R 3.5.0. Versions for all packages are listed below:

* data.table 1.11.2
* devtools 1.13.6
* extrafont 0.17
* ggplot2 3.1.1
* gtools 3.8.1
* Rcpp 1.0.1
* rEDM 0.6.9
* showtext 0.6
* showtextdb 2.0
* sysfonts 0.8

# Quickly Getting Started
1. Open ```run_edm.r```, and set the working directory ```wd``` to the path you save the repository. 
2. Open ```config.r```, and change the variables based on different experimental scenarios.
3. Run all the code in ```run_edm.r```. 

Results of CCM, S-map, and robustness test will be saved in the path you specified in ```config.r```.

# Step by Step Analysis
1. Download raw data from the repository on Zenodo, and put them in ```data``` directory accordingly.
2. Run ```compile_cpue_subarea_data.r```, ```compile_age_data.r```, ```amo.r```, ```sst.r```, ```bottomT.r```, ```fishingM.r```. 

These scripts process the raw data, and save the processed data in ```output``` directory.

3. Run ```combine_all_data_for_edm.r```. 

The script combines the biological data and environmental data together, and saves the combined data (by species name) in ```output``` directory.

4. Open ```config.r```, and change the variables based on different experimental scenarios. 

Because fishing mortality is yearly data, any analyses including fishing mortality should be based on the compiled yearly datasets. Morover, one should use a smaller lag (4 in this study) to do EDM since the short time-series length of yearly data may not be enough to do state-space reconstruction in EDM, especially when the lagged time series are considered.

5. Open ```run_edm.r```. Run all the code before CCM analysis.
6. Run ```run_ccm.r``` to do CCM analysis.

We run CCM 100 times to minimize the effect of bias resulting from the randomness in CCM analysis. Results of CCM will be saved in the path you specified in ```config.r```.

7. Go back to ```run_edm.r```, and continue the code of CCM to process the results of CCM analysis.
8. Run the code of S-map in ```run_edm.r```. 
9. Run the code of robustness test in ```run_edm.r```.
10. Change the variables in ```config.r``` (e.g. dataset, library variable or lags), and repeat step 5-9.


To be continued...

## Notes
If you apply EDM step by step, results might be slightly different in the numeric number because of the randomness in the sampling algorithm used in EDM. However, the results should be qualitatively similar.
