import numpy as np
import matplotlib.pylab as plt
import scipy.constants as spc
import pandas as pd
from scipy import interpolate
import spectrum_line

'''
Defining Pressure and temperatures

pressure_profile_values.csv and temperature_profile_values.csv contain values obtained using digitizer
Interpolate functions are created for both pressures and temp, given an altitude the function
automatically calculates the corresponding pressures and temperatures.
The mean pressure and temperature values are calculated for every 5km.
The mean values are used to calculate the transmittance of a wavenumber.

'''
pressure_profile = np.loadtxt('pressure_profile_values.csv', delimiter=',')
pressure_altitude = pressure_profile[:, 1]
pressure = pressure_profile[:, 0]
pressure_fcn = interpolate.interp1d(pressure_altitude, pressure)

temperature_profile = np.loadtxt('temperature_profile_values.csv', delimiter=',')
temperature_altitude = temperature_profile[:, 1]
temperature = temperature_profile[:, 0]
temperature_fcn = interpolate.interp1d(temperature_altitude, temperature)

'''
Finding the mean pressure and temerature for every 5km
'''

altitude_range = np.arange(2, 99, 5)
mean_of_pressure, mean_of_temperature = ([]), ([])
for val in altitude_range:
    if val == 97:
        mean_range = np.linspace(val, 99, 10)
    else:
        mean_range = np.linspace(val, val+5, 10)
    pressure_mean = pressure_fcn(mean_range).mean()*0.000986923 #[changing from millibars to atm]
    temperature_mean = temperature_fcn(mean_range).mean()
    mean_of_pressure = np.append(mean_of_pressure, pressure_mean)
    mean_of_temperature = np.append(mean_of_temperature, temperature_mean)

'''
Calling the class here
'''
a = spectrum_line.spectrum_line()

'''
Creating a list of HITRAN files
'''
HITRAN_files = ['CO', 'N2', 'SO2', 'CO2', 'H2O', 'O2']


'''
There are four different loops below:
    The first loop, loops over the HITRAN_files containing data regarding different gases
    The second loop, loops over the wavenumbe centers given in the text file
    The third loop, loops over the wavenumbers 0.01cm-1 away from the center
    The fourth loop, loops over the mean pressures and temperatures
'''

for file in HITRAN_files:
    tau_aggregate = ([])
    wavenumbers = ([])
    molecule_id, v_center, S_0, delta_air, alphaLorentz_a, alphaDoppler_s, gamma, E_l = a.getParams(file+'.txt')
    i = 0
    for center in v_center:
        if center > 1000 and center < 6000:
            a.setParams(molecule_id[i], v_center[i], S_0[i], delta_air[i], alphaLorentz_a[i], alphaDoppler_s[i], gamma[i], E_l[i])
            v = np.linspace(center-0.01, center+0.01, 3)
            wavenumbers = np.append(wavenumbers, v)
            for wavenumber in v:
                tau = 1
                for j in range(0, len(mean_of_pressure)):
                    a.setLayersParams(mean_of_pressure[j], mean_of_temperature[j])
                    tau = a.transmittance(wavenumber)*tau
                tau_aggregate = np.append(tau_aggregate, tau)
        print(i)
        i = i+1

    np.savetxt('transmittance_spectrum_atm_'+file+'.csv', tau_aggregate, delimiter=',')
    np.savetxt('wavenumber_atm_'+file+'.csv', wavenumbers, delimiter=',')
