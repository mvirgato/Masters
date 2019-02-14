import numpy as np
import matplotlib
import matplotlib.pyplot as plt

pi = np.pi

Temp =  (1e5/(1.16e4) *1e-9) # Neutron star temperature in GeV
neutron_mass =  (0.939) # mass of neutron in GeV
therm_time =  (1e10 *3.154e7 /(6.58e-16) * 1e9)# in GeV^-1

def coupling_squared(DM_mass):
    k_n = np.sqrt(4*DM_mass*Temp)
    k_0 = DM_mass/3

    bracket = (1/(k_n)**4 - 1/(k_0)**4)
    numerator = 105 * pi**3 * DM_mass * bracket
    denominator = 4* neutron_mass**2 * therm_time

    return numerator/denominator


def cross_section(DM_mass):

    numerator = coupling_squared(DM_mass) * neutron_mass**2 * DM_mass**2
    denominator = pi *(neutron_mass + DM_mass)**2
    cross_section_GeV2 = numerator/denominator

    return cross_section_GeV2 * (1/(1.97e7) *1e-9)**2 * 1e4


def make_plot_cross_section():
    mass_range = np.logspace(-6, 1, num = 1000)
    cross_section_array = np.empty(0)

    for mass in mass_range:
        dummy = cross_section(mass)
        cross_section_array = np.append(cross_section_array, dummy)

    fig, ax1 = plt.subplots(figsize = (10, 7), dpi = 100)
    ax1.loglog(mass_range, cross_section_array, color='black')
    ax1.axis([1e-6, 1e1, 1e-62, 1e-51])
    ax1.set(xlabel = r'Mass of DM [GeV]', ylabel = r'Cross Section [cm$^2$]')
    ax1.fill_between(mass_range, cross_section_array, facecolor='blue', alpha=0.3)
    ax1.grid()
    ax1.text(0.5e-3, 4e-60, 'No thermalization', fontsize=20)
    plt.savefig('Cross Section plot - Analytic.png')


make_plot_cross_section()
