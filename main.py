import sys
import numpy as np
import matplotlib.pyplot as plt
from utils import *

# aerodynamic parameters
alphaDeg = 10.0
vInf = 1.0

# actual code
#nodes = []
#panels = []
#
## nodes coordinates reading
#nodesFile = np.loadtxt('coord_parabolic_20.dat')
#
#count = 1
#for i in nodesFile:
#    nodes.append(Node(np.array([i[0], i[1], 0.0]), count))
#    count += 1
#
#for i in range(len(nodes)-1):
#    panels.append(Panel([nodes[i], nodes[i+1]], i+1))

airfoils = []
airfoils.append(Airfoil('coord_parabolic_front_20.dat'))
airfoils.append(Airfoil('coord_parabolic_rear_30.dat'))
#airfoils.append(Airfoil('coord_flat_plate_10.dat'))
#sys.exit()

# induced velocity on panel i collocation point due to unit vortex on panel j (all airfoils)
for am in airfoils:
    for pi in am.panels:
        for an in airfoils:
            for pj in an.panels:
                pi.unitVelIndCalc(pj.qPoint)

#print(airfoils[0].panels[0].unitVelInd)

#sys.exit()

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

#sys.exit()

# influence coefficients matrix definition
A = np.zeros((N, N))

i = 0
for a in airfoils:
    for pi in a.panels:
        for j in range(N):
            A[i,j] = np.dot(pi.unitVelInd[j], pi.n)

        i = i + 1

#sys.exit()

# linear system solution
gamma = np.linalg.solve(A, b)

# results
i = 0
for a in airfoils:
    for p in a.panels:
        p.setGamma(gamma[i])
        p.dCpCalc(vInf)
        i = i + 1

for aa in airfoils:
    for p in aa.panels:
        for ab in airfoils:
            p.totalIndVelCalc(ab.panels)

        p.liftCalc(vectorInf)

    aa.liftCalc(vInf)

writeResults(airfoils)

#sys.exit()


# plot
# airfoils
#plt.figure(1)
fig, ax = plt.subplots()
for a in airfoils:
    x = [p.nodes[0].coord[0] for p in a.panels]
    x.append(a.panels[-1].nodes[1].coord[0])
    y = [p.nodes[0].coord[1] for p in a.panels]
    y.append(a.panels[-1].nodes[1].coord[1])
    ax.plot(x, y)
    
plt.axis('equal')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.grid()
#plt.show()

# delta Cp along airfoil
f = 2
for a in airfoils:
    # chord vector
    cVec = a.panels[-1].nodes[1].coord - a.panels[0].nodes[0].coord

# airfoil chord
#c = np.linalg.norm(cVec) 
#c = 1.0
#print('airfoil chord: {0:.2f} m\n'.format(c))

    sC = [np.dot(p.mPoint-a.panels[0].nodes[0].coord, cVec)/a.c for p in a.panels]
    dCpPlot = [p.dCp for p in a.panels]
    
    #fig, ax = plt.subplots()
    plt.figure(f)
    plt.plot(sC, dCpPlot, 'k') 
    plt.xlabel('x/c')
    plt.ylabel('dCp')
    #plt.show()
    f = f + 1
plt.show()
# aerodynamic characteristics
# L and Cl
#Gamma = 0.0
#
#for g in gamma:
#    Gamma += g
#
#L = 1.2*vInf*Gamma
#cl = 2.0*L/(1.2*vInf**2*c)


