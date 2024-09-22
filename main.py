###########################################
# coded by Pablo N. Zitelli
# https://www.linkedin.com/in/pablozitelli
###########################################

import sys
import numpy as np
import matplotlib.pyplot as plt
from utils import *

# aerodynamic parameters
alphaDeg = 10.0
vInf = 1.0

# actual code

# airfoils initialization
airfoils = []

airfoils.append(Airfoil('Front', 'coord_parabolic_front_20.dat'))
airfoils.append(Airfoil('Rear', 'coord_parabolic_rear_30.dat'))
#airfoils.append(Airfoil('Flat plate', 'coord_flat_plate_10.dat'))

# induced velocity on panel i collocation point due to unit vortex on panel j (all airfoils)
for am in airfoils:
    for pi in am.panels:
        for an in airfoils:
            for pj in an.panels:
                pi.unitVelIndCalc(pj.qPoint)

# RHS (independent vector) definition
N = 0
for a in airfoils:
    N = N + a.panelNum

b = np.zeros(N)
vectorInf = np.array([vInf*np.cos(alphaDeg*np.pi/180.0), vInf*np.sin(alphaDeg*np.pi/180.0), 0.0])

i = 0
for a in airfoils:
    for p in a.panels:
        b[i] = -np.dot(vectorInf, p.n)
        i = i + 1

# influence coefficients matrix definition
A = np.zeros((N, N))

i = 0
for a in airfoils:
    for pi in a.panels:
        for j in range(N):
            A[i,j] = np.dot(pi.unitVelInd[j], pi.n)

        i = i + 1

# linear system solution
gamma = np.linalg.solve(A, b)

# results
i = 0
for a in airfoils:
    for p in a.panels:
        p.setGamma(gamma[i])
        p.dCpCalc(vInf)
        i = i + 1

# total quantity of panels in system
totalPanels = []
for a in airfoils:
    totalPanels = totalPanels + a.panels

# airfoil lift
for a in airfoils:
    for p in a.panels:
        p.totalIndVelCalc(totalPanels)
        p.liftCalc(vectorInf)

    a.liftCalc(vInf)

writeResults(airfoils)

# plot airfoils
airfoilsNameList = []
fig, ax = plt.subplots()
for a in airfoils:
    airfoilsNameList.append(a.name)
    x = [p.nodes[0].coord[0] for p in a.panels]
    x.append(a.panels[-1].nodes[1].coord[0])
    y = [p.nodes[0].coord[1] for p in a.panels]
    y.append(a.panels[-1].nodes[1].coord[1])
    ax.plot(x, y, label=a.name)
    
plt.axis('equal')
plt.title('airfoils setup')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.legend(airfoilsNameList)
plt.grid()

# delta Cp along airfoils
f = 2
for a in airfoils:
    # chord vector
    cVec = a.panels[-1].nodes[1].coord - a.panels[0].nodes[0].coord

    sC = [np.dot(p.mPoint-a.panels[0].nodes[0].coord, cVec)/a.c for p in a.panels]
    dCpPlot = [p.dCp for p in a.panels]
    
    plt.figure(f)
    plt.title(a.name)
    plt.plot(sC, dCpPlot, 'k') 
    plt.xlabel('x/c')
    plt.ylabel('dCp')
    f = f + 1

plt.show()
