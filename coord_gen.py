import numpy as np

# node qty
N = 51

# initial and final coordinates
iCoord = np.array([0.0, 0.0])
fCoord = np.array([1.0, 0.0])

# airfoil chord
c = np.linalg.norm(fCoord - iCoord)

# dx
dx = c/(N-1)

# camber line function
ep = 0.1

with open('coord_parabolic_50.dat', 'w') as outF:
    for i in range(N):
        # flat plate
        # outF.write('{0:.3f}\t{1:.3f}\n'.format(iCoord[0]+i*dx, 0.0))

        # parabolic camber line
        outF.write('{0:.3f}\t{1:.3f}\n'.format(iCoord[0]+i*dx, -(4*ep/c)*(iCoord[0]+i*dx)**2+4*ep*(iCoord[0]+i*dx)))
