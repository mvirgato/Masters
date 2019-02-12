import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from bertoni_analytic import coupling_squared

pi = np.pi
hbar = cn.hbar
speed_of_light = cn.c
electric_charge = cn.elementary_charge
k_b = cn.Boltzmann

GeV_to_Kelvin = electric_charge/k_b *1e9 # 1GeV in Kelvin
time_conversion = hbar/electric_charge * 1e-9 #1GeV^-1 of time in seconds
length_conversion = (hbar*speed_of_light)/electric_charge *1e-9 # 1GeV^-1 of length in metres

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

    return np.sum(therm_time_array)*time_conversion/3.154e7


def make_plot_time():
    mass_range = np.logspace(-6, 1, num = 100)
    therm_array = np.empty(0)

    for i in mass_range:
        dummy = thermalization_time(i)
        therm_array = np.append(therm_array, dummy)

    fig, ax = plt.subplots(figsize = (10, 7), dpi = 100)
    ax.plot(mass_range, therm_array)
    ax.set_xscale('log')
    #ax.axis([1e-6, 1e1, 1e-62, 1e-51])
    ax.set(xlabel = r'Mass of DM [GeV]', ylabel = r'Thermalization Time')
    plt.savefig('Thermalzation time.png')
    plt.show()

make_plot_time()
