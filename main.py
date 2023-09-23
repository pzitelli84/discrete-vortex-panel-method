import sys
import numpy as np
import matplotlib.pyplot as plt
from utils import *

# aerodynamic parameters
alphaDeg = 10.0
vinf = 1.0

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
airfoils.append(Airfoil('coord_parabolic_20.dat'))
#sys.exit()

# induced velocity on panel i collocation point due to unit vortex on panel j (all airfoils)
for am in airfoils:
    for pi in am.panels:
        for an in airfoils:
            for pj in an.panels:
                pi.velIndCalc(pj.qPoint)

#sys.exit()

# RHS (independent vector) definition
N = 0
for a in airfoils:
    N += a.panelNum

b = np.zeros(N)
vectorInf = np.array([vinf*np.cos(alphaDeg*np.pi/180.0), vinf*np.sin(alphaDeg*np.pi/180.0), 0.0])

ind = 0
for a in airfoils:
    for p in a.panels:
        b[ind] = -np.dot(vectorInf, p.n)
        ind += 1

#sys.exit()

# influence coefficients matrix definition
a = np.zeros((N, N))

i = 0
for am in airfoils:
    for pi in am.panels:
        for j in range(N):
            a[i,j] = np.dot(pi.velInd[j], pi.n)

        i += 1

#sys.exit()

# linear system solution
gamma = np.linalg.solve(a, b)

i = 0
for a in airfoils:
    for p in a.panels:
        p.setGamma(gamma[i])
        p.dCpCalc(vinf)
        i =+ 1

# plot
# chord vector
#cVec = panels[-1].nodes[1].coord - panels[0].nodes[0].coord

# airfoil chord
#c = np.linalg.norm(cVec) 
c = 1.0
print('airfoil chord: {0:.2f} m\n'.format(c))

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
Gamma = 0.0

for g in gamma:
    Gamma += g

L = 1.2*vinf*Gamma
cl = 2.0*L/(1.2*vinf**2*c)

print('Gamma = {0:.4f} m2/s\n'.format(Gamma))
print('L = {0:.4f} N/m\n'.format(L))
print('Cl = {0:.4f}\n'.format(cl))

