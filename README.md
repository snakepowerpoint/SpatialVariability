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
2. Open ```config.r```, and change the parameters based on different experimental scenarios.
3. Run all the code in ```run_edm.r```. 

Results of CCM, S-map, and robustness test will be stored in the path you specified in ```config.r```.

# Step by Step Analysis
1. Download raw data from the repository on Zenodo, and put them in ```data``` directory accordingly.

To be continued...

## Notes
If you apply empirical dynamic modelling (EDM) step by step, results might be slightly different in the numeric number because of the randomness in the sampling algorithm used in EDM. However, the results should be qualitatively similar.



