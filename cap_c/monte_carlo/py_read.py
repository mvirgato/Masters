import numpy as np
import matplotlib.pyplot as plt


nb_dens_file = open('nbdensity.dat')
lines = nb_dens_file.readlines()
radius = np.empty(0)
nb = np.empty(0)


for x in lines:
    x = x.split('\t')
    radius = np.append(radius, float(x[0]))
    nb = np.append(nb, float(x[1]))
nb_dens_file.close

def make_plot():
    fig, ax1 = plt.subplots(figsize = (20, 11), dpi = 500)
    ax1.plot(radius, nb, color='black')
    plt.show()


make_plot()
