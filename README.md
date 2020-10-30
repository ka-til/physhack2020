# Schrodinger's Hack - UWaterloo
Modelling the Venusian atmosphere and generating the absorption spectra given a hypothetical microbe species that has made its home in the CO2 Venusian atmosphere

## Data
CO.txt, N2.txt, SO2.txt, CO2.txt, H2O.txt --- Files downloaded from HITRAN - https://hitran.org/

Venus Atmosphere.xlsx --- Contains parameters likes molecular mass, volume mixing ratio and coefficients to calculate total internal partition sum(TIPS) for different gases

pressure_profile_values.csv --- pressure and altitude data points from digitizer, plot used to collect data points: venus_pressure_profile.jpg

temperature_profile_values.csv --- temperature and altitude data points from digitizer, plot used to collect data points: venus_temperature_profile.jpg

Plots are taken from "Data analysis and simulations of VIRTIS Venus spectra" https://amslaurea.unibo.it/8746/1/magurno_davide_tesi.pdf

## Code

spectral_line.py ---  Python class, when called will generate the transmittance for a single wavenumber

generate_spectral_lines.py --- python code that generates transmittance for multiple spectrum lines and produces two csv files for wavenumbers and their associated transmittance

diffusion_model.ipynb --- IPython notebook for a simple calculation of emitted particles by microbes into the atmosphere.

Test.ipynb --- Where the code was initially tested
