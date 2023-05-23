# Antarctic buoys: a repository of diagnostic tools for analysing sea-ice drifters in the Antarctic marginal ice zone
_Authors: Ashleigh Womack, Marcello Vichi_

This is a repository for the computation of diagnostics related to sea-ice drifters in the Antarctic marginal ice zone (MIZ). The diagnostics have been used for the data analysis in Womack et al. (submitted to The Cryosphere). Preprint DOI:  

## Brief description

This collection of scripts processes drifter data from the Antarctic MIZ, in the Atlantic and Indian Ocean sectors. It has been tested on drifter data collected during the SA Agulhas II expeditions, and public data from the AWI catalogue (see reference in the README.md). These scripts compute the buoys’ drift diagnostics, the drift response to ERA5 atmospheric forcing, the spectral analysis and wavelet analysis of the drift velocity, the cluster absolute and relative dispersion statistics, and the spectral analysis of the deformation proxy.
 

Atmospheric data will need to be downloaded separately from ERA5 reanalysis (hourly on single levels from 1940 to present). Variables should include mean sea level pressure, wind velocity components (u, v), and 2-m air temperature for the period and region of the buoys’ drift. 

This software uses the wavelet module developed by Torrence and Compo (1998), available at https://github.com/chris-torrence/wavelets. Torrence, C., and Compo, G.P., 1998. A practical guide to wavelet analysis. Bull. Amer. Meteor. Soc., 79, 61–78. 

## Workflow 
The workflow is divided into the preprocessing of the data and the analyses. The raw data and preprocessed data are both provided in the [data/](https://github.com/mvichi/antarctic-buoys/tree/main/data) directory of this repository. 

The scripts and functions are provided in the [scripts/](https://github.com/mvichi/antarctic-buoys/tree/main/scripts) directory. 

## Preprocessing of data
The preprocessing file is [Drifter_preprocessing.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Drifter_preprocess.py). This file has been used to generate the files found in the [data/](https://github.com/mvichi/antarctic-buoys/tree/main/data) directory. 

This is the sequence of operations:
1.	Reading in and preprocessing of all buoys
      *	Interpolation of all data to a regular time interval – data with missing, irregular or duplicates
      *	NOTE: Read in [Drifter.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Drifter.py) before the preprocessing
      *	NOTE: ISVP 1 is available at both the 30 minutes and hourly time interval
2.	Computation of geodesic distance
3.	Computation of drift velocity and speed

## Analysis
1.	Drifter diagnostics ([Drifter_diagnostics.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Drifter_diagnostics.py)) 
      *	Drift velocity
      *	Meander coefficient (MC)
      *	Single-buoy absolute dispersion
      *	Trajectory shape
2.	Drifter response to wind and oceanic forcing ([Drifter_diagnostics_metE5.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Drifter_diagnostics_metE5.py))
      *	Wind factor, turning angle, correlation (Pearson and vector), and residual current
      *	Power Spectral Analysis of the drifter and ERA5 wind velocity
3.	Wavelet analysis and Butterworth High-pass filter ([Drifter_wavelet_analysis.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Drifter_wavelet_analysis.py))
      *	NOTE: Read in [Functions_wavelet.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Functions_wavelet.py) file before
4.	Absolute(single-particle) dispersion of the buoy cluster ([Cluster_abs_dispersion.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Cluster_abs_dispersion.py))
5.	Relative (two-particle) dispersion of the buoy cluster ([Cluster_rel_dispersion.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Cluster_rel_dispersion.py))
      *	NOTE: Read in [Functions_rel_dispersion.py](https://github.com/mvichi/antarctic-buoys/blob/main/scripts/Functions_rel_dispersion.py) file before 
      *	Includes the computation of the deformation proxy and its spectral analysis

## Data source
The 2019 Winter and Spring Cruise buoy used in the analysis is published here:

* de Vos, M., Ramjukadh, C., de Villiers, M., Lyttle, C., Womack, A., and Vichi, M. 2023. Polar Iridium Surface Velocity Profilers (p-iSVP), and standard Iridium Surface Velocity Profilers (iSVP) during SCALE 2019 Winter and Spring Cruises. [data set]. Zenodo. https://doi.org/10.5281/zenodo.7954779.
* Womack, A., Verrinder, R., Vichi, M., Alberello, A., MacHutchon, K., Skatulla, S., and Toffoli, A. 2023. An ice-tethered, non-floating Trident Sensors Helix Beacon during SCALE 2019 Spring Cruise. [data set]. Zenodo. https://doi.org/10.5281/zenodo.7954841. 

In addition, a number of other data sources are used in this repository: 

* Eayrs, C., Holland, D., Mojica, J.F., Vichi, M., Alberello, A., Bekker, A., Bennetts, L.G., de Jong, E., Joubert, W., MacHutchon, K., Messori, G., Onorato, M., Saunders, C., Skatulla, S., and Toffoli, A. 2019. SA Agulhas II Winter 2017 Cruise: Waves In Ice Observation Systems (WIIOS). [data set]. Australian Antarctic Data Centre. https://doi.org/10.26179/5cc934992f065.
* MacHutchon, K., Vichi, M., Eayrs, C., Alberello, A., Bekker, A., Bennetts, L.G., Holland, D., de Jong, E., Joubert, W., Messori, G., Mojica, J.F., Onorato, M., Saunders, C., Skatulla, S., and Toffoli, A. 2019. SA Agulhas II Winter 2017 Cruise: ice edge drifters. [data set]. Australian Antarctic Data Centre. https://doi.org/10.26179/5cc937513bd6d.
* Nicolaus, M., Hoppmann, M., Arndt, S., Hendricks, S., Katlein, C., König-Langlo, G., Nicolaus, A., Rossmann, L., Schiller, M., Schwegmann, S., Langevin, D., and Bartsch, A. 2017. Snow height and air temperature on sea ice from Snow Buoy measurements.[data set]. Alfred Wegener Institute, Helmholtz Center for Polar and Marine Research, Bremerhaven. doi: 10.1594/PANGAEA.875638.
* Nicolaus, M., Arndt, S., Hoppmann, M., Krumpen, T., Nicolaus, A., and Bartsch, A. 2017. Sea ice drift, surface temperature, and barometric pressure on sea ice from Surface Velocity Profiler measurements. [data set]. Alfred Wegener Institute, Helmholtz Center for Polar and Marine Research, Bremerhaven. doi: 10.1594/PANGAEA.875652.
* Copernicus Climate Change Service (C3S). 2017. ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate. [data set]. Reading, UK: Copernicus Climate Change Service Climate Data Store (CDS). Retrieved from https://cds.climate.copernicus.eu/cdsapp#!/home.
