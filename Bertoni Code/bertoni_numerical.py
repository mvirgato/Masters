import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as cn


pi = np.pi

Temp =  (1e5/(1.16e4) *1e-9) # Neutron star temperature in GeV
neutron_mass =  (0.939) # mass of neutron in GeV
therm_time =  (1e10 *3.154e7 /(6.58e-16) * 1e9)# in GeV^-1

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
        i +=1

    return np.sum(num_coupling_array)/therm_time#*time_conversion/3.154e7

def numerical_cross_section(DM_mass):
    numerator = numerical_coupling(DM_mass)*neutron_mass**2 * DM_mass**2
    denominator = pi*(neutron_mass + DM_mass)**2
    cross_section_GeV2 = numerator/denominator

    return cross_section_GeV2 *(1/(1.97e7) *1e-9)**2 * 1e4


def make_plot_numeric_cross_section():
    mass_range = np.logspace(-6, 1, num = 1000)
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

make_plot_numeric_cross_section()
