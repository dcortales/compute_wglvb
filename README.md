# Scrips and function for paper Cortés-Morales et al. 2026
This repository provides the codes for obtaining the data and figures presented in Cortes-Morales et al. 2026. It is subdivided in five folders:
1. computation_oliv3
    - `annual_means_ARMOR3D.m`: Annual means of ARMOR3D temperature, salinity, geostrophic meridional velocity and geostrophic zonal velocity
    - `computation_OLIV3_GLOB_annual.m`: Computation beta-plane vertical velocities from ARMOR3D and ERA5 data (OLIV3) and isopycnal interpolation
    - `smooth2a.m`: smoothing function
2. perfect_model_test
    - `ekman_pumping_comp.py`: functions to compute Ekman pumping and beta-plane geostrophic vertical velocity (w_glvb)
    - `OCCITENS_annual_means.py`: Annual averages of Ekman pumping, w_glvb and w_tot
    - `OCCITENS_isopycnal_interpolation.m`: Isopycnal interpolation and spatial smoothing of w_glvb and w_tot
    - `perfect_model_test_metrics.m`: Metrics for the perfect model test (Correlation coefficient, relative error and vertical gradient)
    - `smooth2a.m`: smoothing function
3. validation_1_annual_means
4. validation_2_low_resolution
5. figures

## Data
All the data used in this repository can be accessed via the links provided in Cortés-Morales et al. 2026:
- OCCITENS
- ARMOR3D
- ERA5
- OMEGA3D
- ECCOv4r4
- GLORYS12v1

## Installation and Usage
The neutral density fields for OCCITENS, ARMOR3D and ECCOv4r4, and the geostrophic meridional velocities for OCCITENS has been computed using the repository https://github.com/meom-group/CDFTOOLS.

## Requirements
- Matlab
- Python
