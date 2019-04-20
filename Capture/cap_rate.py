import numpy as np
from scipy import constants as cnts
from scipy import integrate
import matplotlib.pyplot as plt

'''
Global Constants
'''

w = escape_vel = np.sqrt( (2 * cnts.G * 1.4 * 2e30) / 1e3 )

temp = (1e7/ 1.16e4) * 1e-9

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
    if abs(init_energy(s, t, dm)/fermi_energy) > 10:
        if (init_energy(s, t, dm) - fermi_energy) > 0:
            return 0
        elif  (init_energy(s, t, dm) - fermi_energy) < 0:
            return 1
    else:
        return 1 / (1 + np.exp((init_energy(s, t, dm) - fermi_energy) / temp))

def final_FD(s, t, v, dm):
    if abs(final_energy(s, t, v, dm)/fermi_energy) > 10:
        if (final_energy(s, t, v, dm) - fermi_energy) > 0:
            return 0
        elif  (final_energy(s, t, v, dm) - fermi_energy) < 0:
            return 1
    else:
        return 1 / (1 + np.exp((-final_energy(s, t, v, dm) + fermi_energy) / temp))

'''
trig functions
'''
def incoming_cos(s, t):
    return (w**2 - s**2 -t**2 ) / (2 * s * t)

def outgoing_cos(s, t, v):
    return (v**2 - s**2 -t**2 ) / (2 * s * t)

'''
Step functions
'''
def incoming_step(s, t):
    if abs(incoming_cos(s, t)) > 1:
        return 0
    elif abs(incoming_cos(s, t)) <= 1:
        return 1

def outgoing_step(s, t, v):
    if abs(outgoing_cos(s, t, v)) > 1:
        return 0
    elif abs(outgoing_cos(s, t, v)) <= 1:
        return 1

def heaviside_product(s, t, v):
    return incoming_step(s, t) * outgoing_step(s, t, v)


'''
Integrand
'''
def integrand(s, t, v, dm):
    return v * t   * heaviside_product(s, t, v) * init_FD(s, t, dm) * final_FD(s, t, v, dm)


'''
main
'''
def init_FD_test(E):
    if abs(E/fermi_energy) > 100:
        if (E - fermi_energy) > 0:
            return 0
        elif  (E - fermi_energy) < 0:
            return 1
    else:
        return 1 / (1 + np.exp((E - fermi_energy) / temp))

def final_FD_test(E):
    if abs(E/fermi_energy) > 10:
        if (E - fermi_energy) > 0:
            return 0
        elif  (E - fermi_energy) < 0:
            return 1
    else:
        return 1 / (1 + np.exp((-E + fermi_energy) / temp))


def make_plot():
    v_range = np.linspace(0.08499, 0.08501, 1000)
    dist1 = np.empty(0)
    dist2 = np.empty(0)

    for x in v_range:
        dummy1 = init_FD_test(x)
        dist1 = np.append(dist1, dummy1)

        dummy2 = final_FD_test(x)
        dist2 = np.append(dist2, dummy2)


    fig, ax1 = plt.subplots()
    ax1.plot(v_range, dist1)
    ax1.plot(v_range, dist2)
    #ax1.set_xscale('log')
    plt.show()

print(init_FD_test(-100*fermi_energy))