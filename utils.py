import numpy as np

# panel class definition
class Panel:

    def __init__(self, nodes, name):
        self.velInd = []
        self.nodes = nodes
        self.id = name

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

    def velIndCalc(self, quarterPoint):
        vector =  self.cPoint - quarterPoint 
        r = np.linalg.norm(vector)
        u = (0.5/(np.pi*r**2))*vector[1]
        v = - (0.5/(np.pi*r**2))*vector[0]
        self.velInd.append(np.array([u, v, 0.0]))

    def setGamma(self, gamma):
        self.gamma = gamma

    def dCpCalc(self, velInf):
        self.dCp = 2*self.gamma/(velInf*self.len)

# node class definition
class Node:

    def __init__(self, coord, name):
        self.coord = coord
        self.id = name

