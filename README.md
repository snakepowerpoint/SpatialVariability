# Causal effects of population dynamics and environmental changes on spatial variability of marine fishes

This repository is an implementaion of the research work **Causal effects of population dynamics and environmental changes on spatial variability of marine fishes.** 

- For readers who want to quickly implement the empirical dynamic modelling (EDM), please refer to **Quickly Getting Started**. 
- For readers who want to implement the EDM from scratch, please refer to **Step by Step Analysis**.

The repository contains:
* Source code

Raw data including *fishery survey data*, *Atlantic Multidecadal Oscillation (AMO)*, *sea surface temperature (SST)*, *sea bottom tempeature (SBT)*, and *fishing mortality* can be downloaded on Zenodo (https://doi.org/10.5281/zenodo.3695703) or at the links provided in the paper.

# Requirements 
R program with version 3.5.0. R is a free software and publicly available at https://www.r-project.org/. To install R, please follow the procedures below:
1. Choose an CRAN Mirrors 
2. Choose your operating system
3. Download binaries for base distribution
4. Install R accordingly

Previous versions are available at https://cran.r-project.org/bin/windows/base/old/. All R-packages and the corresponding versions used in this study are listed below:

* data.table 1.11.2
* devtools 1.13.6
* gtools 3.8.1
* Rcpp 1.0.1
* rEDM 0.6.9
* showtextdb 2.0
* ggplot2 3.1.1
* ggmap 3.0.0
* maps 3.3.0
* extrafont 0.17
* sysfonts 0.8
* showtext 0.6
* ncdf4 1.16
* dplyr 0.8.0.1
* akima 0.6-2
* rnaturalearth 0.1.0
* rnaturalearthdata 0.1.0
* scatterpie 0.1.2

Please install these R-packages via instruction `install.packages(xxx)` in R, where `xxx` is package name.

# Quickly Getting Started
1. Open `run_edm.r` in R, and set the working directory `wd` to the path where you save the repository.
2. Change config in `run_emd.r` based on different experimental scenarios.

    In this study, we have four experimental scenarios:

    ```
    # First
    source("config_cv.r")
    ```
    ```
    # Second
    source("config_age_F.r")
    ```
    ```
    # Third
    source("config_abundance_F.r")
    ```
    ```
    # Fourth
    source("config_cv_F.r")
    ```
3. Run all the code in `run_edm.r` for each experimental scenario.

    Results of CCM, S-map, and robustness test will be saved in the path specified in config file (i.e., `ccm_path`, `smap_path`, and `robust_path`). Before conducting a new experimental scenario, remember to remove all defined variables.

# Step by Step Analysis
1. Download raw data from the repository on Zenodo, and put them in `data` directory accordingly.

2. Run `compile_cpue_subarea_data.r`, `compile_age_data.r`, `amo.r`, `sst.r`, `bottomT.r` and `fishingM.r` to process raw data.

    These scripts process the raw data, and save the compiled data in `output` directory.

3. Run `combine_all_data_for_edm.r`. 

    The script combines the biological data and environmental data together, and saves the compiled data (named by species name) in `output` directory.

4. Open `rum_edm.r`, and change the config based on different experimental scenarios. 

    As mentioned in **Quickly Getting Started** above, we have four experimental scenarios in this study (i.e. `config_cv.r`, `config_age_F.r`, `config_abundance_F.r` and `config_cv_F.r`).
    
    Because fishing mortality is yearly data, any analyses including fishing mortality should be based on the compiled yearly dataset. Morover, we use a smaller lag (4 in this study) to perform EDM on yearly dataset because the short time series in yearly dataset may not be enough for state-space reconstruction, especially when the lagged time series are considered.

5. Run all the code before CCM analysis part in `run_edm.r`.

6. Run `run_ccm.r` to perform CCM analysis.

    We run CCM 200 times to mitigate the effect of bias resulting from the randomness in CCM analysis. Results of CCM will be saved in the path specified in config file.

7. Go back to `run_edm.r`, and continue the code of CCM to process the results of CCM analysis.

8. Run the code of S-map in `run_edm.r` to perform S-map analysis. 

9. Run the code of robustness test in `run_edm.r` to perform robustness test.

# Contact
If you find any bugs or have any questions about the implementation, pelease contact us via r03241220@ntu.edu.tw
