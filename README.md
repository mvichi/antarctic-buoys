# antarctic-buoys
# Introduction
This is a collection of scripts to process drifter data from the Antarctic marginal ice zone in the Atlantic and Indian sectors.
These scripts compute the drift diagnostics, drift response to atmospheric forcing, spectral analysis and wavelet analysis of drift velocity, absolute and relative dispersion statisctics, PSD of deformation proxy. 

Atmospheric data will need to be separtely downloaded from Era5 reanalysis (hourly on single levels from 1940 to present). 

# Used datasets
* Eayrs, C., Holland, D., Mojica, J.F., Vichi, M., Alberello, A., Bekker, A., Bennetts, L.G., de Jong, E., Joubert, W., MacHutchon, K., Messori, G., Onorato, M., Saunders, C., Skatulla, S., Toffoli, A., 2019. SA Agulhas II Winter 2017 Cruise: Waves In Ice Observation Systems (WIIOS). https://doi.org/10.26179/5cc934992f065
* MacHutchon, K., Vichi, M., Eayrs, C., Alberello, A., Bekker, A., Bennetts, L.G., Holland, D., de Jong, E., Joubert, W., Messori, G., Mojica, J.F., Onorato, M., Saunders, C., Skatulla, S., Toffoli, A., 2019. SA Agulhas II Winter 2017 Cruise: ice edge drifters. https://doi.org/10.26179/5cc937513bd6d
* Nicolaus, M.; Hoppmann, M.; Arndt, S.; Hendricks, S.; Katlein, C.; König-Langlo, G.; Nicolaus, A.; Rossmann, L.; Schiller, M.; Schwegmann, S.; Langevin, D.; Bartsch, A. (2017): Snow height and air temperature on sea ice from Snow Buoy measurements. Alfred Wegener Institute, Helmholtz Center for Polar and Marine Research, Bremerhaven, doi:10.1594/PANGAEA.875638.
* Nicolaus, M.; Arndt, S.; Hoppmann, M.; Krumpen, T.; Nicolaus, A.; Bartsch, A. (2017): Sea ice drift, surface temperature, and barometric pressure on sea ice from Surface Velocity Profiler measurements. Alfred Wegener Institute, Helmholtz Center for Polar and Marine Research, Bremerhaven, doi:10.1594/PANGAEA.875652.
* Copernicus Climate Change Service (C3S) (2017) ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate. Reading, UK: Copernicus Climate Change Service Climate Data Store (CDS). Retrieved from https://cds.climate.copernicus.eu/cdsapp#!/home.

# External software
Wavelet software was provided by C. Torrence and G. Compo, and is available at URL: http://atoc.colorado.edu/research/wavelets/
Torrence, C., and G. P. Compo, 1998: A practical guide to wavelet analysis. Bull. Amer. Meteor. Soc., 79, 61–78.
The python functions are available at:
https://github.com/chris-torrence/wavelets
