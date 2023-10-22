import numpy as np

# node class definition
class Node:
    def __init__(self, coord, name):
        self.coord = coord
        self.id = name

# panel class definition
class Panel:
    def __init__(self, nodes, name):
        self.unitVelInd = []
        #self.totalVelInd = []
        self.nodes = nodes
        self.id = name
        self.lift = 0.0

        # tangent unit vector and panel length
        tVec = self.nodes[1].coord - self.nodes[0].coord
        self.len = np.linalg.norm(tVec)
        self.t = tVec/self.len

        # normal unit vector
        kVersor = np.array([0.0, 0.0, 1.0])
        nVec = np.cross(kVersor, self.t)
        nMag = np.linalg.norm(nVec)
        self.n = nVec/nMag

        # quarter chord point coordinates
        self.qPoint = self.nodes[0].coord + 0.25*tVec

        # collocation point coordinates
        self.cPoint = self.nodes[0].coord + 0.75*tVec

        # midpoint coordinates
        self.mPoint = self.nodes[0].coord + 0.5*tVec

    # calculation of induced velocity produced by unit vortex 
    def unitVelIndCalc(self, quarterPoint):
        vector =  self.cPoint - quarterPoint 
        r = np.linalg.norm(vector)
        u = (0.5/(np.pi*r**2))*vector[1]
        v = - (0.5/(np.pi*r**2))*vector[0]
        self.unitVelInd.append(np.array([u, v, 0.0]))

    def setGamma(self, gamma):
        self.gamma = gamma

    # delta Cp calculation across panel
    def dCpCalc(self, velInf):
        self.dCp = 2*self.gamma/(velInf*self.len)

    # total velocity induced by all vortices
    def totalIndVelCalc(self, panels):
        self.totalVelInd = np.array([0.0, 0.0, 0.0])
        for p in panels:
            distVec =  self.qPoint - p.qPoint 
            r = np.linalg.norm(distVec)

            if r != 0.0:
                u = p.gamma*(0.5/(np.pi*r**2))*distVec[1]
                v = - p.gamma*(0.5/(np.pi*r**2))*distVec[0]
                self.totalVelInd = self.totalVelInd + np.array([u, v, 0.0])

    # panel lift calculation (generalized Kutta-Joukowski)
    def liftCalc(self, velInf):
        vInf = np.linalg.norm(velInf)
        self.lift = 1.2*vInf*self.gamma*(1 + np.dot(velInf, self.totalVelInd)/vInf**2)

# airfoil class definition
class Airfoil:
    def __init__(self, coordFile):
        self.panels = []

        # nodes coordinates reading
        nodesFile = np.loadtxt(coordFile)

        # nodes creation
        nodes = []
        count = 1
        for i in nodesFile:
            nodes.append(Node(np.array([i[0], i[1], 0.0]), count))
            count = count + 1

        # panels creation
        for i in range(len(nodes)-1):
            self.panels.append(Panel([nodes[i], nodes[i+1]], i+1))

        self.panelNum = len(self.panels)

        # chord calculation
        cVec = self.panels[-1].nodes[1].coord - self.panels[0].nodes[0].coord
        self.c = np.linalg.norm(cVec) 

    # aerodynamic characteristics: Gamma, L and Cl
    def liftCalc(self, vInf):
        self.Gamma = 0.0
        self.L = 0.0

        for p in self.panels:
            self.Gamma = self.Gamma + p.gamma
            self.L = self.L + p.lift

        self.Cl = 2.0*self.L/(1.2*vInf**2*self.c)

# domain (grid) class definition for streamlines
class Domain:
    def __init__(self, nx, ny, xMin, xMax, yMin, yMax):

        self.grid = []

        x = np.linspace(xMin, xMax, nx)
        y = np.linspace(yMin, yMax, ny)

# grid point class definition
class Point:
    def __init__(self, x, y):
        self.coord = np.array([x, y, 0.0])
        
# write results to file and terminal
def writeResults(airfoils):
    with open('results.dat', 'w') as out:
        print(40*'=' + '\n')
        out.write(40*'=' + '\n')

        for a in airfoils:
            print('Airfoil #{0:d} aerodynamic characteristics:'.format(list(airfoils).index(a)+1))
            print('Gamma = {0:.4f} m2/s'.format(a.Gamma))
            print('L = {0:.4f} N/m'.format(a.L))
            print('Cl = {0:.4f}\n'.format(a.Cl))
            print(40*'=' + '\n')
            
            out.write('Airfoil #{0:d} aerodynamic characteristics:\n'.format(list(airfoils).index(a)+1))
            out.write('Gamma = {0:.4f} m2/s\n'.format(a.Gamma))
            out.write('L = {0:.4f} N/m\n'.format(a.L))
            out.write('Cl = {0:.4f}\n'.format(a.Cl))
            out.write(40*'=' + '\n')
