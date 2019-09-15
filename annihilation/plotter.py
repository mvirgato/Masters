import numpy as np
import matplotlib.pyplot as plt



def ann_plots():

    lines = []
    domain = np.empty(0)
    range1 = np.empty(0)
    range2 = np.empty(0)
    range3 = np.empty(0)
    range4 = np.empty(0)
    range5 = np.empty(0)
    range6 = np.empty(0)

    current_file = open('d8ACS.dat', 'r')
    lines = current_file.readlines()

    for x in lines:
        x = x.split('\t')
        domain = np.append(domain, float(x[0]))
        range1 = np.append(range1, float(x[1]))
        range2 = np.append(range2, float(x[2]))
        range3 = np.append(range3, float(x[3]))
        range4 = np.append(range4, float(x[4]))
        range5 = np.append(range5, float(x[5]))
        range6 = np.append(range6, float(x[6]))

    current_file.close

    fig, ax1 = plt.subplots()
    ax1.loglog(domain, range1)
    ax1.loglog(domain, range2)
    ax1.loglog(domain, range3)
    ax1.loglog(domain, range4)
    ax1.loglog(domain, range5)
    ax1.loglog(domain, range6)
    ax1.set(xlabel = r'$m_\chi$ [GeV]', ylabel = r'$G^2$ $[GeV^{-2}]$')
    # ax1.set_yscale('log')
    # ax1.axis([0, 12, 0, 1.2e56])
    plt.savefig('d8ACS.pdf')


ann_plots()
