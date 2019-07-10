import numpy as np
import matplotlib.pyplot as plt



files = ['misner_full_cap', 'misner_no_vel_cap', 'vegas_no_vel_cap', 'vegas_full']

def cap_rate_plots():
    """
    Plots capture rate
    """
    lines = []
    domain = np.empty(0)
    range = np.empty(0)
    range2 = np.empty(0)
    # range3 = np.empty(0)

    current_file = open('cap_rate_rad.dat', 'r')
    lines = current_file.readlines()


    for x in lines:
        x = x.split('\t')
        domain = np.append(domain, float(x[0]))
        range = np.append(range, float(x[1]))
        range2 = np.append(range2, float(x[2]))
        # range3 = np.append(range3, float(x[3]))

    current_file.close

    fig, ax1 = plt.subplots()
    ax1.plot(domain, range, color='blue')
    ax1.set(xlabel = r'$Radius$ [km]', ylabel = r'$C/V$')
    # ax1.set_yscale('log')
    # ax1.axis([0, 12, 0, 1.2e56])
    plt.savefig('C_per_V.png')

    fig, ax2 = plt.subplots()
    ax2.plot(domain, range2, color='blue')
    ax2.set(xlabel = r'$Radius$ [km]', ylabel = r'$zeta n$')
    # ax2.set_yscale('log')
    # ax1.axis([0, 11.5, 0.4, 1])
    plt.savefig('zeta_times_n.png')

    # fig, ax3 = plt.subplots()
    # ax3.plot(domain, range3, color='blue')
    # ax3.set(xlabel = r'$Radius$ [km]', ylabel = r'$C/V$')
    # # ax3.set_yscale('log')
    # # ax1.axis([0, 11.5, 0.4, 1])
    # plt.savefig('no_density_cap.png')

def cap_interp_plts():

    lines = []
    domain = np.empty(0)
    range = np.empty(0)

    current_file = open('cap_interp_test.dat', 'r')
    lines = current_file.readlines()

    for x in lines:
        x = x.split('\t')
        domain = np.append(domain, float(x[0]))
        range = np.append(range, float(x[1]))

    current_file.close

    fig, ax1 = plt.subplots()
    ax1.plot(domain, range, color='blue')
    ax1.set(xlabel = r'$Radius$ [km]', ylabel = r'$C/V$')
    # ax1.set_yscale('log')
    # ax1.axis([0, 12, 0, 1.2e56])
    plt.savefig('cap_interp_plot.png')


def full_plotter():

    for name in files:
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
        ax1.loglog(domain, range, color='blue')
        ax1.set(xlabel = r'$m$ [GeV]', ylabel = r'$C$', title = name)
        # ax1.set_yscale('log')
        # ax1.axis([0, 12, 0, 1.2e56])
        plt.savefig(name + '.png')

full_plotter()
cap_rate_plots()
