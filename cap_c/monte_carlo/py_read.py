import numpy as np
import matplotlib.pyplot as plt


current_file = open('cap_rate_MC.dat', 'r')
lines = current_file.readlines()
mass = np.empty(0)
integraleval = np.empty(0)


for x in lines:
    x = x.split('\t')
    mass = np.append(mass, float(x[0]))
    integraleval = np.append(integraleval, float(x[1]))
current_file.close

def make_plot():
    fig, ax1 = plt.subplots(figsize = (20, 11   ), dpi = 500)
    ax1.loglog(mass, integraleval, color='black')
    plt.savefig('initial_integrals.png')


make_plot()
