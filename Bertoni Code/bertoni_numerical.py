import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from bertoni_analytic import coupling_squared

pi = np.pi

GeV_to_Kelvin = 1.16e4 *1e9 # 1GeV in Kelvin
time_conversion = 6.58e-16 * 1e-9 #1GeV^-1 of time in seconds
length_conversion = 1.97e-7 *1e-9 # 1GeV^-1 of length in metres

Temp = 1e5/GeV_to_Kelvin # Neutron star temperature in GeV
neutron_mass = 0.939 # mass of neutron in GeV
therm_time = 1e10 *3.154e7 /time_conversion # in GeV^-1

def energy_averge(previous_energy):
    return (3/7) * previous_energy

def scattering_rate(DM_mass, energy):
    numerator = 2*coupling_squared(DM_mass) * neutron_mass**2 * DM_mass * energy**2
    denominator = 15* pi**3
    return numerator/denominator

def initial_energy(DM_mass):
    return 1.05*DM_mass



def thermalization_time(DM_mass):
    energy_array = np.empty(0)
    therm_time_array = np.empty(0)
    energy_array = np.append(energy_array, initial_energy(DM_mass))
    energy_array = np.append(energy_array, energy_averge(energy_array[0]))


    i=1
    while (energy_averge(energy_array[i-1] - energy_array[i]) > Temp):
        therm_time_array = np.append(therm_time_array, 1/scattering_rate(DM_mass, energy_array[i]))
        new_energy = energy_averge(energy_array[i])
        energy_array = np.append(energy_array, new_energy)
        i += 1 
    
    return np.sum(therm_time_array)

print(thermalization_time(1)*time_conversion/3.154e7)