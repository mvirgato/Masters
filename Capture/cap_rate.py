import numpy as np
from scipy import constants as cnts
from scipy import integrate
import matplotlib.pyplot as plt

'''
Global Constants
'''

w = escape_vel = np.sqrt( (2 * cnts.G * 1.4 * 2e30) / 1e3 )

temp = (1e3/ 1.16e4) * 1e-9

therm_time = (1e10 * 3.145e7 * 1e9) / 6.58e-16

nm = 0.939

fermi_energy = 0.085

vd = 270e3

vs = 200e3

'''
Set up Functions
'''

def mu(dm):
    return dm/nm

def muplus(dm):
    return (1 + mu(dm)) / 2

def init_vel_sqr(s, t, dm):
    return (2 * mu(dm) * muplus(dm) * t**2 + 2* muplus(dm) * s**2 - mu(dm) * w**2) / cnts.speed_of_light**2

def final_vel_sqr(s, t, v, dm):
    return (2 * mu(dm) * muplus(dm) * t**2 + 2* muplus(dm) * s**2 - mu(dm) * v**2) / cnts.speed_of_light**2

def init_energy(s, t, dm):
    return 0.5 * nm * init_vel_sqr(s, t, dm)

def final_energy(s,t,v,dm):
    return 0.5 * nm * final_vel_sqr(s,t,v,dm)

def init_FD(s, t, dm):
    try:
        a = 1 / (1 + np.exp((init_energy(s, t, dm) - fermi_energy) / temp))
        return a

    except OverflowError:
        if init_energy(s, t, dm) > fermi_energy:
            return 0
        elif init_energy(s, t, dm) < fermi_energy:
            return 1

def final_FD(s, t, v, dm):
    try:
        a = 1 / (1 + np.exp((final_energy(s, t, v, dm) - fermi_energy) / temp))
        return a

    except OverflowError:
        if final_energy(s, t, v,  dm) > fermi_energy:
            return 0
        elif final_energy(s, t, v, dm) < fermi_energy:
            return 1




def make_plot():
    u_range = np.linspace(0.0849985, 0.0850010, 10000)
    dist = np.empty(0)

    for mass in u_range:
        dummy = (mass)
        dist = np.append(dist, dummy)

    fig, ax1 = plt.subplots()
    ax1.plot(u_range, dist)
    #ax1.set_xscale('log')
    plt.show()

print(10**400)