import dadi
from dadi import PhiManip, Integration
import numpy

def parallel_domestication((nu1F,T1),(n1,n2),pts):

    """
    Model with two splits, with bottlenecks in pop 2 and 3, and migration between all three populations
    
    nu1F: ancestral population size after growth (nu=pop size 1=population 1 F=final)(initial size is defined to be 1)
    T1: the scaled time between ancestral population growth and split 1
    """

    #create grid
    xx = yy =dadi.Numerics.default_grid(pts)

    #phi for equilibrium ancestral population
    phi=PhiManip.phi_1D(xx)

    # Ancestral population growth
    phi = Integration.one_pop(phi, xx, T=Tp, nu=nu1F)

    #first split of population 2 from 1
    phi=PhiManip.phi_1D_to_2D(xx,phi)

    #split of population 3 from 1
    phi=PhiManip.phi_2D_to_3D_split_1(xx,phi)

    #calculate spectrum
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))

    return sfs
