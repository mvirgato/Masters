import numpy as np
import matplotlib.pyplot as plt

file_names = ['mass_ns', 'muFnchempot', 'nbdensity', 'Ynabund', 'esc_vel', 'n_density']



def eos_plots( ):
    """
    Plots all quantities in EoS
    """

    for name in file_names:

        lines = []
        domain = np.empty(0)
        range = np.empty(0)

        current_file = open(name + '.dat', 'r')
        lines = current_file.readlines()


        for x in lines:
            x = x.split('\t')
            domain = np.append(domain, float(x[0]))
            range = np.append(range, float(x[1]) )


        current_file.close

        fig, ax1 = plt.subplots()
        ax1.plot(domain, range, color='blue')
        ax1.set(xlabel = r'$radius$ [km]', ylabel = name)
        # ax1.set_yscale('log')
        plt.savefig(name + '.png')

def single_plot( chosen_file ):
    """
    Plots chosen data
    """

    lines = []
    domain = np.empty(0)
    range = np.empty(0)

    current_file = open( chosen_file + '.dat', 'r')
    lines = current_file.readlines()

    for x in lines:
        x = x.split('\t')
        domain = np.append(domain, float(x[0]))
        range = np.append(range, float(x[1]))

    current_file.close

    fig, ax1 = plt.subplots()
    ax1.plot(domain, range, color='blue')
    ax1.set(xlabel = r'$r/R$ [km]', ylabel = chosen_file)
    ax1.set_yscale('log')
    plt.savefig(chosen_file + '.png')



def cap_rate_plots():
    """
    Plots capture rate
    """
    lines = []
    domain = np.empty(0)
    range = np.empty(0)

    current_file = open('cap_rate_MC.dat', 'r')
    lines = current_file.readlines()


    for x in lines:
        x = x.split('\t')
        domain = np.append(domain, float(x[0]))
        range = np.append(range, float(x[1]))

    current_file.close

    fig, ax1 = plt.subplots()
    ax1.loglog(domain, range, color='blue')
    ax1.set(xlabel = r'$m$ [GeV]', ylabel = r'$C$')
    ax1.axis([1e-6, 1e2, 1e14, 1e38])
    plt.savefig('muFnchempot.png')

single_plot('muFnchempot')
