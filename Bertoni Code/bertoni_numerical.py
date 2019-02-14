import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as cn
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

def scattering_rate_no_coupling(DM_mass, energy):
    numerator = 2* neutron_mass**2 * DM_mass * energy**2
    denominator = 15* pi**3
    return numerator/denominator

def initial_energy(DM_mass):
    return 1.05*DM_mass



def numerical_coupling(DM_mass): #returns therm in 1/Gev units
    energy_array = np.empty(0)
    num_coupling_array = np.empty(0)
    energy_array = np.append(energy_array, initial_energy(DM_mass))
    energy_array = np.append(energy_array, energy_averge(energy_array[0]))


    i=1
    while (energy_averge(energy_array[i-1] - energy_array[i]) > Temp):
        num_coupling_array = np.append(num_coupling_array, 1/scattering_rate_no_coupling(DM_mass, energy_array[i]))
        new_energy = energy_averge(energy_array[i])
        energy_array = np.append(energy_array, new_energy)
        i += 1

    return np.sum(num_coupling_array)/therm_time#*time_conversion/3.154e7

def numerical_cross_section(DM_mass):
    numerator = numerical_coupling(DM_mass)*neutron_mass**2 * DM_mass**2
    denominator = pi*(neutron_mass + DM_mass)**2
    cross_section_GeV2 = numerator/denominator

    return cross_section_GeV2 *(length_conversion**2) * (100**2)


def make_plot_numeric_coupling():
    mass_range = np.logspace(-6, 1, num = 100)
    numerical_cross_section_array = np.empty(0)

    for i in mass_range:
        dummy = numerical_cross_section(i)
        numerical_cross_section_array = np.append(numerical_cross_section_array, dummy)

    fig, ax2 = plt.subplots(figsize = (10, 7), dpi = 100)
    ax2.loglog(mass_range, numerical_cross_section_array)
    ax2.axis([1e-6, 1e1, 1e-62, 1e-51])
    ax2.set(xlabel = r'Mass of DM [GeV]', ylabel = r'Cross Section [cm$^2$]')
    plt.savefig('Cross Section plot - numerical.png')
    plt.show()


def make_thing():
    mass_range = np.logspace(-6, 1, num=1000)
    array = np.empty(0)

    for mass in mass_range:
        dummy = numerical_coupling(mass)
        array = np.append(array, dummy)

    plt.plot(mass_range, array)
    plt.xscale('log')
    plt.show()

make_thing()
