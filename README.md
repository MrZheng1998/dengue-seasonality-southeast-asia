# dengue-seasonality-southeast-asia
Processed data and analysis code for dengue seasonality phenotypes in Southeast Asia.
# Dengue seasonality phenotypes in Southeast Asia


This repository contains processed data and analysis code for the manuscript:
"Latitude-associated climatic gradients characterize distinct dengue seasonal phenotypes in Southeast Asia".


## Repository structure


- data/processed/: processed analytical datasets
- code/: R scripts for DTW clustering, peak timing, XGBoost-SHAP, GAM, sensitivity analyses, and figure generation
- docs/: variable descriptions and data dictionary


## Data sources


Raw data were obtained from publicly available sources, including OpenDengue, World Bank Open Data, ASEAN Statistics Division, NCBI Virus, and meteorological data sources.


## How to reproduce the analysis


Run the scripts in the following order:


1. 01_dtw_clustering.R
2. 02_xgboost_shap_analysis.R
3. 03__gam_threshold_analysis.R



## Software environment


R version: R 4.2.3
Key packages: tidyverse, dtwclust, xgboost, SHAPforxgboost, mgcv, brms


## Citation


Please cite the manuscript and this repository if using these data or code.
