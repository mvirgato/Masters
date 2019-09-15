import numpy as np
import matplotlib.pyplot as plt

names = []

for i in range(11)[1:]:
    names.append('d'+ str(i))

def ann_plots():

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
