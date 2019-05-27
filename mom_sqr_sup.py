import numpy as np
from scipy import constants as cnts
from scipy import integrate
from scipy import special
import matplotlib.pyplot as plt


SMALL_SIZE = 30
MEDIUM_SIZE = 40
BIGGER_SIZE = 50

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

pi = np.pi




Temp = (1e3 / (1.16e4) * 1e-9)  # Neutron star temperature in GeV
neutron_mass = (0.939)  # mass of neutron in GeV
therm_time = (1e10 * 3.154e7 / (6.58e-16) * 1e9)  # in GeV^-1
escape_vel = np.sqrt( (2 * cnts.G * 1.4 * 2e30) / 10e3 ) / (3e8)
v_d = 270e3 / 3e8
ns_vel = 200e3 / 3e8


def coupling_squared(DM_mass):
    # k_n = np.sqrt(4*DM_mass*Temp)
    k_0 = DM_mass/ 3.0

    bracket = (1/(4.0 * DM_mass * Temp)**3 - 1/(k_0)**6)
    numerator = 63.0 * pi**3 * DM_mass**3
    denominator = neutron_mass**2 * therm_time

    return (numerator/denominator) * bracket


def cross_section(DM_mass):

    second_moment = (v_d**2 + ns_vel**2 + escape_vel**2 )

    numerator = coupling_squared(DM_mass) * DM_mass**2 * second_moment * neutron_mass**4
    denominator = 2.0* pi * (DM_mass  +  neutron_mass)**4
    cross_section_GeV2 = numerator / denominator

    return cross_section_GeV2 * (1/(1.97e7) *1e-9)**2 * 1e4

def make_plot_cross_section():
    mass_range = np.logspace(-6, 6, num = 1000)
    cross_section_array = np.empty(0)

    for mass in mass_range:
        dummy = cross_section(mass)
        cross_section_array = np.append(cross_section_array, dummy)

    fig, ax1 = plt.subplots(figsize = (20, 11), dpi = 500)
    ax1.loglog(mass_range, cross_section_array, color='black')
    ax1.axis([1e-6, 1e6, 1e-52, 1e-40])
    ax1.set(xlabel = r'$m_\chi$ [GeV]', ylabel = r'$\sigma$ [cm$^2$]')
    ax1.grid(linestyle='--')
    ax1.fill_between(mass_range, cross_section_array, facecolor='forestgreen', alpha = 0.9)
    ax1.text(1e-2, 4e-48, 'No thermalization', fontsize=40)
    ax1.text(5e2, 5e-44, r'$G[\bar{\chi}\gamma^5\chi][\bar{\psi}\psi]$', fontsize=40)
    plt.savefig('mom_2.png')


make_plot_cross_section()
