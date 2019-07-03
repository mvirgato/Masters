import numpy as np
import matplotlib.pyplot as plt


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
    # ax1.set_yscale('log')
    plt.savefig(chosen_file + '.png')



def cap_rate_plots():
    """
    Plots capture rate
    """
    lines = []
    domain = np.empty(0)
    range = np.empty(0)

    current_file = open('cap_rate_rad.dat', 'r')
    lines = current_file.readlines()


    for x in lines:
        x = x.split('\t')
        domain = np.append(domain, float(x[0]))
        # range = np.append(range, float(x[1]))
        range = np.append(range, float(x[2]))

    current_file.close

    fig, ax1 = plt.subplots()
    ax1.plot(domain, range, color='blue')
    ax1.set(xlabel = r'$Radius$ [km]', ylabel = r'$C/V$')
    # ax1.axis([0, 11.5, 0.4, 1])
    # plt.savefig('cap_rate_plot.png')
    plt.savefig('zeta.png')

cap_rate_plots()
