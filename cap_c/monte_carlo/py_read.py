import numpy as np
import matplotlib.pyplot as plt

file_names = ['mass_ns', 'muFnchempot', 'nbdensity', 'Ynabund']



for name in file_names:

    lines = []
    domain = np.empty(0)
    range = np.empty(0)

    current_file = open(name + '.dat', 'r')
    lines = current_file.readlines()


    for x in lines:
        x = x.split('\t')
        domain = np.append(domain, float(x[0]))
        range = np.append(range, float(x[1]))


    current_file.close

    fig, ax1 = plt.subplots()
    ax1.plot(domain, range, color='blue')
    ax1.set(xlabel = r'$radius$ [km]', ylabel = name)
    # ax1.set_yscale('log')
    plt.savefig(name + '.png')
