# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 20:51:40 2020

@author: Patrick
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
import pandas as pd

def input_validation(variable):
    try:
        value = float(variable)
    except ValueError:
        print('The input must be a number')
        return (0, False)
    if value < 0:
        print('The input must be a number')
        return (0, False)
    else:
        return (float(variable), True)


def euler(m, k, b, nsteps, xe, ve):
    for i in range(nsteps - 1):
        # calculate a at ste i
        # no array of these as these are not needed for later
        ae = -(k / m) * xe[i] - (b / m) * ve[i]
        # add next value for x and v
        xe[i + 1] = xe[i] + ve[i] * h
        ve[i + 1] = ve[i] + ae * h

    eulerdata = np.column_stack((xe, ve))

    if path.exists('EulerData.csv'):
        f = open('EulerData.csv', 'w')
    else:
        f = open('EulerData.csv', 'x')

    df = pd.DataFrame(data=eulerdata)
    df.to_csv('EulerData.csv')
    
    f.close()
    
    plt.figure()
    plt.title('euler')
    plt.plot(xe, 'b')
    plt.plot(ve, 'r')
    plt.show()

    plt.figure()
    plt.title('euler')
    plt.plot(xe, ve, 'b')
    plt.show()




def improv_euler(m, k, b, nsteps, xie, vie):
    for i in range(nsteps - 1):
        # calculate a at ste i
        # no array of these as these are not needed for later
        aie = -(k / m) * xie[i] - (b / m) * vie[i]
        # add next value for x and v
        xie[i + 1] = xie[i] + vie[i] * h + 0.5 * h ** 2 * aie
        vie[i + 1] = vie[i] + aie * h

    impeulerdata = np.column_stack((xie, vie))

    if path.exists('ImprovedEulerData.csv'):
        f = open('ImprovedEulerData.csv', 'w')
    else:
        f = open('ImprovedEulerData.csv', 'x')
            
    df = pd.DataFrame(data=impeulerdata)
    df.to_csv('ImprovedEulerData.csv')
    
    f.close()
    
    plt.figure()
    plt.title('Ieuler')
    plt.plot(xie, 'b')
    plt.plot(vie, 'r')
    plt.show()

    plt.figure()
    plt.title('Ieuler')
    plt.plot(xie, vie, 'b')
    plt.show()



def verlet(m, k, b, nsteps, xv, vv):
    D = 2 * m + b * h
    A = 2 * ((2 * m - k * h ** 2) / D)
    B = (b * h - 2 * m) / D

    # use euler for first step
    xv[1] = xv[0] + (h * vv[0])

    for i in range(nsteps - 2):
        xv[i + 2] = A * xv[i + 1] + B * xv[i]
        vv[i + 1] = (xv[i + 2] - xv[i]) / (2 * h)

        verletdata = np.column_stack((xv, vv))

        if path.exists('VerletData.csv'):
            f = open('VerletData.csv', 'w')
        else:
            f = open('VerletData.csv', 'x')

    df = pd.DataFrame(data=verletdata)
    df.to_csv('VerletData.csv')
    
    f.close()
    
    plt.figure()
    plt.title('Verlet')
    plt.plot(xv, 'b')
    plt.plot(vv, 'r')
    plt.show()

    plt.figure()
    plt.title('Verlet')
    plt.plot(xv, vv, 'b')
    plt.show()




def euler_cormer(m, k, b, nsteps, xec, vec):
    for i in range(nsteps - 1):
        aec = -(k / m) * xec[i] - (b / m) * vec[i]
        vec[i + 1] = vec[i] + aec * h
        xec[i + 1] = xec[i] + h * vec[i + 1]


    eulercromerdata = np.column_stack((xec, vec))

    if path.exists('EulerCromerData.csv'):
        f = open('EulerCromerData.csv', 'w')
    else:
        f = open('EulerCromerData.csv', 'x')

    df = pd.DataFrame(data=eulercromerdata)
    df.to_csv('EulerCromerData.csv')
    
    f.close()
    
    plt.figure()
    plt.title('EC')
    plt.plot(xec, 'b')
    plt.plot(vec, 'r')
    plt.show()

    plt.figure()
    plt.title('EC')
    plt.plot(xec, vec, 'b')
    plt.show()



def analytical(m, b, h, k, T):
    gamma = b / m
    determinant = gamma**2 - 4 * k / m
    time_array = np.arange(0, T, h)
    A = 0
    B = 0
    C = 0
    D = 0

    if determinant == 0: #Condition for critical damping
        C = 1/gamma
        D = - C

        xa = C * np.exp((-gamma/2) * time_array) + D * np.exp((-gamma/2) * time_array)
        va = C * (-gamma / 2) * np.exp((-gamma/2) * time_array) + D * (-gamma/2) * np.exp((-gamma/2) * time_array)
    else: #Non-critical damping version
        alpha1 = (-gamma + np.sqrt(complex(determinant))) / 2
        alpha2 = (-gamma - np.sqrt(complex(determinant))) / 2
        A = -1 / (alpha1 - alpha2)
        B = - A

        xa = A * np.exp(alpha1 * time_array) + B * np.exp(alpha2 * time_array)
        va = A * alpha1 * np.exp(alpha1 * time_array) + B * alpha2 * np.exp(alpha2 * time_array)



    analyticaldata = np.column_stack((xa, va))
    analyticaldata = analyticaldata.real
    

    if path.exists('AnalyticalData.csv'):
        f = open('AnalyticalData.csv', 'w')
    else:
        f = open('AnalyticalData.csv', 'x')

    df = pd.DataFrame(data=analyticaldata)
    df.to_csv('AnalyticalData.csv')
    
    f.close()
    
    plt.figure()
    plt.title('Analytical')
    plt.plot(xa, 'b')
    plt.plot(va, 'r')
    plt.show()

    plt.figure()
    plt.title('Analytical')
    plt.plot(xa, va, 'b')
    plt.show()


# user inputs of physical parameters


m = 3.82
k = 0.61
b = 0.9
# valid = False
#
# while valid == False:
#     m, valid = input_validation(input('What is the value of the mass'))
#
# valid = False
# while valid == False:
#     k, valid = input_validation(input('What is the value of the spring constant'))
#
# valid = False
# while valid == False:
#     b, valid = input_validation(input('What is the value of the damping'))

# length of integration
T = 221
# step size
h = 1
# num of steps (i) needs to be an int
nsteps = int(T / h)

x = np.zeros(nsteps)
v = np.zeros(nsteps)
x[0] = 0
v[0] = -1


euler(m, k, b, nsteps, x, v)

improv_euler(m, k, b, nsteps, x, v)

verlet(m, k, b, nsteps, x, v)

euler_cormer(m, k, b, nsteps, x, v)

analytical(m, b, h, k, T)

t = np.arange(0, T, h)
#t= t[:-1]

EulerData = np.zeros((0, 2))
ImpEulerData = np.zeros((0, 2))
VerletData = np.zeros((0, 2))
EulerCromerData = np.zeros((0, 2))
AnalyticalData = np.zeros((0, 2))


df = pd.read_csv('EulerData.csv', index_col=0)
EulerData = df.to_numpy()

df = pd.read_csv('ImprovedEulerData.csv', index_col = 0)
ImpEulerData = df.to_numpy()

df = pd.read_csv('VerletData.csv', index_col=0)
VerletData = df.to_numpy()

df = pd.read_csv('EulerCromerData.csv', index_col=0)
EulerCromerData = df.to_numpy()

df = pd.read_csv('AnalyticalData.csv', index_col=0)
AnalyticalData = df.to_numpy()





EEuler = 0.5 * k * EulerData[:,0] ** 2 + 0.5 * m * EulerData[:,1] ** 2
EIEuler = 0.5 * k * ImpEulerData[:,0] ** 2 + 0.5 * m * ImpEulerData[:,1] ** 2
EVerlet = 0.5 * k * VerletData[:,0] ** 2 + 0.5 * m * VerletData[:,1] ** 2
EEc = 0.5 * k * EulerCromerData[:,0] ** 2 + 0.5 * m * EulerCromerData[:,1] ** 2
EAnalytical = 0.5 * k * AnalyticalData[:,0] ** 2 + 0.5 * m * AnalyticalData[:,1] ** 2

#EAnalytical = EAnalytical[:-1]

plt.figure()
plt.plot(t, EEuler, label='Euler')
plt.plot(t, EIEuler, label='Improved Euler')
plt.plot(t, EVerlet, label='Verlet')
plt.plot(t, EEc, label='Euler Cromer')
plt.plot(t, EAnalytical, label='Analytical')
plt.xlabel('Time (s)')
plt.ylabel('Energy(J)')
plt.legend()
plt.show()
