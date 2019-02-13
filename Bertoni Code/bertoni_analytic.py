import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as cn

pi = np.pi
hbar = cn.hbar
speed_of_light = cn.c
electric_charge = cn.elementary_charge

GeV_to_Kelvin = electric_charge/k_b *1e9 # 1GeV in Kelvin
time_conversion = hbar/electric_charge * 1e-9 #1GeV^-1 of time in seconds
length_conversion = (hbar*speed_of_light)/electric_charge *1e-9 # 1GeV^-1 of length in metres

Temp = 1e5/GeV_to_Kelvin # Neutron star temperature in GeV
neutron_mass = 0.93956 # mass of neutron in GeV
therm_time = 1e10 *3.154e7 /time_conversion # in GeV^-1

def coupling_squared(DM_mass):

    k_n = np.sqrt(4*DM_mass*Temp)
    k_0 = DM_mass/3


    bracket = 1/(k_n)**4 - 1/(k_0)**4
    numerical_factor = (105*(pi**3))/(4*(neutron_mass**2)*therm_time)
    return numerical_factor * DM_mass * bracket

def cross_section(DM_mass):
    numerator = coupling_squared(DM_mass)*neutron_mass**2 * DM_mass**2
    denominator = pi*(neutron_mass + DM_mass)**2
    cross_section_GeV2 = numerator/denominator

    return cross_section_GeV2 *(length_conversion**2) * (100**2)


def make_plot():
    mass_range = np.logspace(-6, 1, num = 1000)
    cross_section_array = np.empty(0)

    for i in mass_range:
        dummy = cross_section(i)
        cross_section_array = np.append(cross_section_array, dummy)

    fig, ax = plt.subplots(figsize = (10, 7), dpi = 100)
    ax.loglog(mass_range, cross_section_array)
    ax.axis([1e-6, 1e1, 1e-62, 1e-51])
    ax.set(xlabel = r'Mass of DM [GeV]', ylabel = r'Cross Section [cm$^2$]')
    plt.savefig('Cross Section plot - Analytic.png')
    plt.show()

make_plot()
