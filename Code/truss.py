"""
class implementation of a 2D/3D truss class

creator: Peter Mackenzie-Helnwein
date: Oct 29, 2019
revised: April 5, 2020
"""

# import libraries
import numpy as np

# the element master class
from Element import *


# class definition
class TrussElement(Element):
    """
    !class: TrussElement

    variables:
        self.N0  ... unit normal vector from 1 to 2
        self.n  ... unit normal vector from 1 to 2

    inherited variables:
        self.nnode ......... number of nodes per element
        self.ndof .......... number of degrees of freedom per node
        self.X = (X1,X2) ... tuple of nodal position vectors (as np.array)
        self.U ............. array of nodal displacement vectors
        self.force ......... internal forces at nodes
        self.stiffness ..... the stiffness matrix

    overloaded methods:
        def init(self) ....... element specific initialization steps
        compute(self) ...... does the actual computation

    inherited methods:
        __init__(self, X, params)
        setDisp(self, U) ... U is an array of nodal displacement vectors
        getFe(self) ........ return the internal force vector as array of nodal vectors
        getKe(self) ........ return the stiffness matrix as array of nodal matrices
        getFeAsMatrix(self) .. return the internal force vector as nx1 matrix
        getKeAsMatrix(self) .. return the stiffness matrix as nxn matrix
    """

    def init(self):

        # make sure all needed parameters exist
        if not 'E' in self.params.keys():
            self.params['E'] = 1.0

        if not 'A' in self.params.keys():
            self.params['A'] = 1.0

        # compute reference base vectors
        L0vec = self.X[1] - self.X[0]
        L02 = L0vec @ L0vec
        self.L0 = np.sqrt(L02)
        self.N0 = L0vec / self.L0


    def compute(self):

        # compute deformed base vectors
        x1 = self.X[0] + self.U[0]
        x2 = self.X[1] + self.U[1]

        lvec = x2 - x1
        l2 = lvec @ lvec
        l  = np.sqrt(l2)
        n = lvec / l

        # compute strain
        stretch = l / self.L0
        strain = np.log(stretch)

        # compute internal force
        k = self.params['E'] * self.params['A']
        fe = k * strain

        # compute nodal force
        self.force = (-fe * n, fe * n)

        # compute tangent stiffness
        #ke = (k ) / l * np.outer(n, n)
        ke = (k - fe)/l * np.outer(n, n) + fe/l * np.identity(self.ndof)
        self.stiffness = np.array([[ke,-ke],[-ke,ke]])

# defining main execution procedure

def main():
    # create a demo element
    X11 = np.array([0.,0.])
    X12 = np.array([5.5,0.5])
    e1 = TrussElement((X11,X12))

    X21 = np.array([9.5,0.5])
    X22 = np.array([5.5,0.5])
    e2 = TrussElement((X21,X22))

    # now set some displacement and check out changes
    e1.setDisp((np.array([0.0,0.0]),np.array([0.0,-.3])))
    # # print('=======')
    print('U = ',e1.U)
    # # print('---')
    # print('Fe = ',e1.getFe())
    # # print('Kte = \n',e1.getKe())
    # # print('---')
    print('Fe = ',e1.getFeAsMatrix())
    # # print('Kte = \n',e1.getKeAsMatrix())

    e2.setDisp((np.array([0.0, 0.0]), np.array([0.0, -.3])))
    # print('=======')
    print('U = ',e2.U)
    # print('---')
    # print('Fe = ',e2.getFe())
    # print('Kte = \n',e2.getKe())
    # print('---')
    print('Fe = ',e2.getFeAsMatrix())
    # print('Kte = \n',e2.getKeAsMatrix())




# main execution ****************************************

if __name__ == "__main__":
    main()
    sys.exit(0)
