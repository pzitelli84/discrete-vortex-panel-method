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
airfoils.append(Airfoil('coord_parabolic_rear_20.dat'))
#airfoils.append(Airfoil('coord_flat_plate_10.dat'))
#sys.exit()

# induced velocity on panel i collocation point due to unit vortex on panel j (all airfoils)
for am in airfoils:
    for pi in am.panels:
        for an in airfoils:
            for pj in an.panels:
                pi.velIndCalc(pj.qPoint)

#print(airfoils[0].panels[0].velInd)

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
            A[i,j] = np.dot(pi.velInd[j], pi.n)

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

    a.liftCalc(vInf)

    print('Airfoil #{0:d} aerodynamic characteristics:\n'.format(list(airfoils).index(a)+1))
    print('Gamma = {0:.4f} m2/s\n'.format(a.Gamma))
    print('L = {0:.4f} N/m\n'.format(a.L))
    print('Cl = {0:.4f}\n'.format(a.Cl))

#sys.exit()


# plot
# chord vector
#cVec = panels[-1].nodes[1].coord - panels[0].nodes[0].coord

# airfoil chord
#c = np.linalg.norm(cVec) 
#c = 1.0
#print('airfoil chord: {0:.2f} m\n'.format(c))

#sC = [np.dot(p.mPoint-panels[0].nodes[0].coord, cVec)/c for p in panels]
#dCpPlot = [p.dCp for p in panels]
#
#fig, ax = plt.subplots()
#ax.plot(sC, dCpPlot, 'k') 
#plt.xlabel('x/c')
#plt.ylabel('dCp')
#plt.show()

# aerodynamic characteristics
# L and Cl
#Gamma = 0.0
#
#for g in gamma:
#    Gamma += g
#
#L = 1.2*vInf*Gamma
#cl = 2.0*L/(1.2*vInf**2*c)


