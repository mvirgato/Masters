import numpy as np
import matplotlib.pyplot as plt


pi = np.pi

Temp = (1e5 / (1.16e4) * 1e-9)  # Neutron star temperature in GeV
neutron_mass = (0.939)  # mass of neutron in GeV
therm_time = (1e10 * 3.154e7 / (6.58e-16) * 1e9)  # in GeV^-1

mass_array = np.logspace(-6, 1, 500)

def energy_averge(previous_energy):
    return (3 / float(7) * previous_energy) 


def scattering_rate_no_coupling(DM_mass, energy):
    numerator = 2 * neutron_mass**2 * DM_mass * energy**2
    denominator = 15 * pi**3
    return numerator/float(denominator)


def initial_energy(DM_mass):
    return 1.05*DM_mass


def Energy_List(DM_mass):
    energies = np.empty(0)
    energies = np.append(energies, initial_energy(DM_mass))
    energies = np.append(energies, energy_averge(energies[0]))

    i = 0
    while(energy_averge(energies[i] -energies[i+1]) > Temp):
        i += 1
        energies = np.append(energies, energy_averge(energies[i]))

    return energies


def numerical_coupling(DM_mass):
    coupling_array = np.empty(0)
    scatt_array = np.zeros(0)
    temp_energy = Energy_List(DM_mass)

    for energy in temp_energy:
        scatt_array = np.append(scatt_array, scattering_rate_no_coupling(DM_mass, energy))

    for scatter in scatt_array:
        dummy = 1/float(therm_time * scatter)
        coupling_array = np.append(coupling_array, dummy)

    return sum(coupling_array)


def cross_section_numeric(DM_mass):
    numerator = numerical_coupling(DM_mass) * neutron_mass**2 * DM_mass
    denominator = pi * (neutron_mass + DM_mass)**2

    return numerator / float(denominator)


def plot_coupling():
    temp_array = np.zeros(0)

    for mass in mass_array:
        temp_array = np.append(temp_array, numerical_coupling(mass))

    plt.loglog(mass_array, temp_array)
    plt.show()

def plot_energy():
    temp_array = np.empty(0)

    for mass in mass_array:
        temp_array = np.append(temp_array, 1/sum(Energy_List(mass)))
    
    plt.loglog(temp_array, mass_array)
    plt.show()


def plot_cross_section():
    temp_array = np.zeros(0)

    for mass in mass_array:
        temp_array = np.append(temp_array, cross_section_numeric(mass))

    plt.loglog(mass_array, temp_array)
    plt.show()


plot_cross_section()