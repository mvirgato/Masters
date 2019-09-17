import numpy as np
import matplotlib.pyplot as plt

plt.style.use('science')

SMALL_SIZE = 10
MEDIUM_SIZE = 14
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

names = []

for i in range(11)[1:]:
    names.append('d'+ str(i))

def ann_plots(figsize = (20, 11), dpi = 500):

    for j in names:
        lines = []

        mydict = {}

        current_file = open(j + 'ACS.dat', 'r')
        lines = current_file.readlines()

        for i in range(len(lines[0].split('\t'))):
            mydict['range'+ str(i)] = np.empty(0)

        for x in lines:
            x = x.split('\t')

            for i in range(len(x)):
                mydict['range'+ str(i)] = np.append(mydict['range'+ str(i)], float(x[i]))

        fig, ax1 = plt.subplots()
        for i in range(len(lines[0].split('\t')))[1:]:
            ax1.loglog(mydict['range0'], mydict['range'+ str(i)])
        ax1.set(xlabel = r'$m_\chi$ [GeV]', ylabel = r'$G^2$ $[GeV^{-2}]$')
        plt.savefig(j + 'ACS.pdf')



ann_plots()
