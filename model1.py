import dadi
from dadi import PhiManip
import numpy

def split((param),(n1,n2),pts):

"""
Model with two splits, with bottlenecks in pop 2 and 3, and migration between all three populations

define parameters here
"""

#create grid
pts=20
xx= dadi.Numerics.default_grid(pts)

#phi for equilibrium ancestral population
phi=PhiManip.phi_1D(xx)

#first split of population 2 from 1
phi=PhiManip.phi_1D_to_2D(xx,phi)

#split of population 3 from 1
phi=PhiManip.phi_2D_to_3D_split_1(xx,phi)



