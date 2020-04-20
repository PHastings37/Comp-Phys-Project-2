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
    """
    Validates the three inputs of the program makes sure that they are numbers 
    above 0
    """
    
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
    """
    Function for calculating Euler's method
    """
    for i in range(nsteps - 1):
        # calculate a at ste i
        # no array of these as these are not needed for later
        ae = -(k / m) * xe[i] - (b / m) * ve[i]
        # add next value for x and v
        xe[i + 1] = xe[i] + ve[i] * h
        ve[i + 1] = ve[i] + ae * h

    eulerdata = np.column_stack((xe, ve))
    #makes it easier to manage if theres only one array being handled

    if path.exists('EulerData.csv'):
        f = open('EulerData.csv', 'w')
    else:
        f = open('EulerData.csv', 'x')

    df = pd.DataFrame(data=eulerdata)
    df.to_csv('EulerData.csv')
    
    f.close()
    

    plt.figure()
    plt.title('euler')
    plt.plot(xe, 'b', label='position')
    plt.plot(ve, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title('euler')
    plt.plot(xe, ve, 'b')
    plt.show()





def improv_euler(m, k, b, nsteps, xie, vie):
    """
    Function for plotting the improved Euler's method
    """
    
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
    plt.plot(xie, 'b', label='position')
    plt.plot(vie, 'r', label='velocity')
    plt.legend()
    plt.show()
  
    plt.figure()
    plt.title('Ieuler')
    plt.plot(xie, vie, 'b')
    plt.show()



def verlet(m, k, b, nsteps, xv, vv):
    """
    Function for plotting Verlet's method
    """
    D = 2 * m + b * h
    A = 2 * ((2 * m - k * h ** 2) / D)
    B = (b * h - 2 * m) / D
    a_nought = -(k / m) *xv[0] - (b / m) * vv[0]
    # use euler-cromer for first step, gives the a higher accuracy than imp Euler or Euler
    xv[1] = xv[0] + (h * vv[0]) + 0.5 * h**2 * a_nought

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
    plt.plot(xv, 'b', label='position')
    plt.plot(vv, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title('Verlet')
    plt.plot(xv, vv, 'b')
    plt.show()
 




def euler_cromer(m, k, b, nsteps, xec, vec):
    """
    Function for plotting the Euler-Cromer method
    """
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
    plt.plot(xec, 'b', label='position')
    plt.plot(vec, 'r', label='velocity')
    plt.legend()
    plt.show()
  
    plt.figure()
    plt.title('EC')
    plt.plot(xec, vec, 'b')
    plt.show()
  
 

def analytical(m, b, h, k, T):
    """
    Function for the analytical method
    """
    
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
        va = A * alpha1 * np.exp(alpha1 * time_array)+ B * alpha2 * np.exp(alpha2 * time_array)



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
    plt.plot(xa, 'b', label='position')
    plt.plot(va, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title('Analytical')
    plt.plot(xa, va, 'b')
    plt.show()



# user inputs of physical parameters



F = 0.1 # can be anything, looks good using this value though

valid = False

while valid == False:
     m, valid = input_validation(input('What is the value of the mass'))

valid = False
while valid == False:
     k, valid = input_validation(input('What is the value of the spring constant'))

valid = False
while valid == False:
     b, valid = input_validation(input('What is the value of the damping'))

# length of integration
T = 100
# step size
h = 1
"""
The length and step size can be varied, but this combination shows the methods
come to an almost complete stop before the integration ends, as well as 
demonstrating the difference in accuracys of the different methods
"""

# num of steps (i) needs to be an int
nsteps = int(T / h)

x = np.zeros(nsteps)
v = np.zeros(nsteps)
x[0] = 0
v[0] = -1


euler(m, k, b, nsteps, x, v)

improv_euler(m, k, b, nsteps, x, v)

verlet(m, k, b, nsteps, x, v)

euler_cromer(m, k, b, nsteps, x, v)

analytical(m, b, h, k, T)

t = np.arange(0, T, h)

frequency = 1 # frequency was varied here to investigate the resonance
F_sinusoidal = 20 * np.sin(frequency * t)

#initalises all the arrays for the energy calculations
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


"""
Finding the best method
Done by finding the average distance between the analytical method and each of 
the numerical methods, Verlet is the best method, with an average distance 
nearly 100 times smaller than the next best one(Euler-Cromer)
"""

energy_diffs = []

energy_diffs = np.append(energy_diffs, np.mean(np.abs(EAnalytical - EEuler)))

energy_diffs = np.append(energy_diffs, np.mean(np.abs(EAnalytical - EIEuler)))

energy_diffs = np.append(energy_diffs, np.mean(np.abs(EAnalytical - EVerlet)))

energy_diffs = np.append(energy_diffs, np.mean(np.abs(EAnalytical - EEc)))


#Small algorithm to determine which is the most accurate
minimum = np.min(energy_diffs)
for i in range(len(energy_diffs)):
    if energy_diffs[i] == minimum:
        pos = i
        break
        
"""
The below investigates the use of forces and pushes with different methods, 
but only for the 'best method' it was found that for extereme cases for high 
step size and damping, the best method changed, so this will still investigate 
for the appropriate method
"""
if pos == 0:
    xe = np.zeros(nsteps)
    ve = np.zeros(nsteps)
    print('The most accurate method is Euler')
    euler(m, k, np.sqrt(k * m), nsteps, x, v)
    
    euler(m, k, 2 * np.sqrt(k * m), nsteps, x, v)
    
    euler(m, k, 4 * np.sqrt(k * m), nsteps, x, v)
    
    for i in range(nsteps - 1):
        # Only adds the force in the region specified, resumes as before after
        if (nsteps / 2) < i < ((5 * nsteps) / 8):
            #The +F/m is the external force being added
            ae = -(k / m) * xe[i] - (b / m) * ve[i] + F/m 
            # add next value for x and v
            xe[i + 1] = xe[i] + ve[i] * h
            ve[i + 1] = ve[i] + ae * h
        else:
             ae = -(k / m) * xe[i] - (b / m) * ve[i] 
             # add next value for x and v
             xe[i + 1] = xe[i] + ve[i] * h
             ve[i + 1] = ve[i] + ae * h
        
    plt.figure()     
    plt.title('Forced euler')
    plt.plot(xe, 'b', label='position')
    plt.plot(ve, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title(' Forced euler')
    plt.plot(xe, ve, 'b')
    plt.show()
    
    for i in range(nsteps -1):
          ae = -(k / m) * xe[i] - (b / m) * ve[i] + F_sinusoidal[i]/m
          # add next value for x and v
          xe[i + 1] = xe[i] + ve[i] * h
          ve[i + 1] = ve[i] + ae * h
    
    plt.figure()     
    plt.title('Sinusoidal Forced euler')
    plt.plot(xe, 'b', label='position')
    plt.plot(ve, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title('Sinusoidal Forced euler')
    plt.plot(xe, ve, 'b')
    plt.show()
elif pos == 1:
    xie = np.zeros(nsteps)
    vie = np.zeros(nsteps)
    print('The most accurate method is Improved Euler')
    improv_euler(m, k, np.sqrt(k * m), nsteps, x, v)
    
    improv_euler(m, k, 2 * np.sqrt(k * m), nsteps, x, v)
    
    improv_euler(m, k, 4 * np.sqrt(k * m), nsteps, x, v)
    
    for i in range(nsteps - 1):
        # calculate a at ste i
        # no array of these as these are not needed for later
        if (nsteps / 2) < i < ((5 * nsteps) / 8):
            aie = -(k / m) * xie[i] - (b / m) * vie[i]
            # add next value for x and v
            xie[i + 1] = xie[i] + vie[i] * h + 0.5 * h ** 2 * aie + F/m
            vie[i + 1] = vie[i] + aie * h
        else:
            aie = -(k / m) * xie[i] - (b / m) * vie[i] 
            # add next value for x and v
            xie[i + 1] = xie[i] + vie[i] * h + 0.5 * h ** 2 * aie
            vie[i + 1] = vie[i] + aie * h
    plt.figure()     
    plt.title('Forced Improved euler')
    plt.plot(xie, 'b', label='position')
    plt.plot(vie, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title(' Forced Improved euler')
    plt.plot(xie, vie, 'b')
    plt.show()
    
    for i in range(nsteps - 1):
         aie = -(k / m) * xie[i] - (b / m) * vie[i] + F_sinusoidal[i] / m
         # add next value for x and v
         xie[i + 1] = xie[i] + vie[i] * h + 0.5 * h ** 2 * aie
         vie[i + 1] = vie[i] + aie * h
         
    plt.figure()     
    plt.title('Sinusodial Forced Improved euler')
    plt.plot(xie, 'b', label='position')
    plt.plot(vie, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title(' Forced Improved euler')
    plt.plot(xie, vie, 'b')
    plt.show()
        
elif pos == 2:
    xv = np.zeros(nsteps)
    vv = np.zeros(nsteps)
    xv[0] = 0
    vv[0] = -1
    
    print('The most accurate method is Verlet')
    """
    the below code prints graphs for Verlet's method for the cases of half, 
    double and exactly crtical damping
    """
    verlet(m, k, np.sqrt(k * m), nsteps, x, v)
    
    verlet(m, k, 2 * np.sqrt(k * m), nsteps, x, v)
    
    verlet(m, k, 4 * np.sqrt(k * m), nsteps, x, v)
    
    
    D = 2 * m + b * h
    A = 2 * ((2 * m - k * h ** 2) / D)
    B = (b * h - 2 * m) / D
    a_nought = -(k / m) *xv[0] - (b / m) * vv[0]
    # use euler for first step
    xv[1] = xv[0] + (h * vv[0]) + 0.5 * h**2 * a_nought

    for i in range(nsteps - 2):
        if (nsteps / 2) < i < ((5 * nsteps) / 8):
            xv[i + 2] = A * xv[i + 1] + B * xv[i] + ((2 * h**2) / D) * F
            vv[i + 1] = (xv[i + 2] - xv[i]) / (2 * h)
        else:
            xv[i + 2] = A * xv[i + 1] + B * xv[i] 
            vv[i + 1] = (xv[i + 2] - xv[i]) / (2 * h)

    """
    The above code works  to add a force for a few oscillations in the middle 
    of the time. It is observed that the postion jumps drastically while the 
    velocity does change, but only slightly. after the 'push', the graph 
    returns to normal very quickly without the driving force. The same was done 
    with a negative force, and the exact same was seen but inverted
    """
    plt.figure()     
    plt.title('Forced Verlet')
    plt.plot(xv, 'b', label='position')
    plt.plot(vv, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title(' Forced Verlet')
    plt.plot(xv, vv, 'b')
    plt.show()
        
    
    """
    This last part plots a graph is there is a sinusoidal force acting 
    throughout the duration of the run time, it is seen that at the start there
    is a transient period where the amplitude of the position and velocity 
    varies, but after a few oscillations it becomes constant. Measuring this 
    amplitude after the transient period can be used to investigate resonance
    """
    for i in range(nsteps - 12):
        xv[i + 2] = A * xv[i + 1] + B * xv[i] + ((2 * h**2) / D) * F_sinusoidal[i]
        vv[i + 1] = (xv[i + 2] - xv[i]) / (2 * h)

    plt.figure()     
    plt.title('Sinusoidal Forced Verlet')
    plt.plot(xv, 'b', label='position')
    plt.plot(vv, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title('Sinusoidal Forced Verlet')
    plt.plot(xv, vv, 'b')
    plt.show()
    
else:
    xec = np.zeros(nsteps)
    vec = np.zeros(nsteps)
    print('The most accurate method is Euler Cromer')
    euler_cromer(m, k, np.sqrt(k * m), nsteps, x, v)
    
    euler_cromer(m, k, 2 * np.sqrt(k * m), nsteps, x, v)
    
    euler_cromer(m, k, 4 * np.sqrt(k * m), nsteps, x, v)
    
    for i in range(nsteps - 1):
        if (nsteps / 2) < i < ((5 * nsteps) / 8):
            aec = -(k / m) * xec[i] - (b / m) * vec[i] + F/m
            vec[i + 1] = vec[i] + aec * h
            xec[i + 1] = xec[i] + h * vec[i + 1]
        else:
            aec = -(k / m) * xec[i] - (b / m) * vec[i] 
            vec[i + 1] = vec[i] + aec * h
            xec[i + 1] = xec[i] + h * vec[i + 1]
            
    plt.figure()     
    plt.title('Forced euler cromer')
    plt.plot(xec, 'b', label='position')
    plt.plot(vec, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title(' Forced euler cromer')
    plt.plot(xec, vec, 'b')
    plt.show()

    for i in range(nsteps - 1):
        aec = -(k / m) * xec[i] - (b / m) * vec[i] + F_sinusoidal[i] / m
        vec[i + 1] = vec[i] + aec * h
        xec[i + 1] = xec[i] + h * vec[i + 1]

    plt.figure()     
    plt.title('Sinusoidal Forced euler cromer')
    plt.plot(xec, 'b', label='position')
    plt.plot(vec, 'r', label='velocity')
    plt.legend()
    plt.show()
 
    plt.figure()
    plt.title('Sinusoidal Forced euler cromer')
    plt.plot(xec, vec, 'b')
    plt.show()
