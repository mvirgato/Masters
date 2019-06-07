import numpy as np
import matplotlib.pyplot as plt


current_file = open('eos_24_lowmass.dat', 'r')
lines = current_file.readlines()

del lines[0]

radius = np.empty(0)
mass = np.empty(0)

for x in lines:
    x = x.split('\t')
    radius = np.append(radius, float(x[0]))
    mass = np.append(mass, float(x[1]))


current_file.close
