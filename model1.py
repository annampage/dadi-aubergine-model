import dadi
from dadi import PhiManip, Integration
import numpy

def parallel_domestication((nu1F, nu2B, nu2F, nu3B, nu3F, T1, T2, T3, T2B, T3B, m12, m13, m21, m23, m31, m32),(n1,n2,n3),pts):

    """
    Model with two splits, with bottlenecks in pop 2 and 3, and migration between all three populations
    
    nu1F: ancestral population size after growth (initial size is defined to be 1)
    nu2B: Bottleneck size in population 2
    nu2F: Final population size for population 2
    nu3B: Bottleneck size in population 3
    nu3F: Final population size in population 3
    T1: the scaled time between ancestral population growth and split 1
    T2: time between split 1 and split 2
    T3: time between split 3 and present
    T2B: length of bottleneck in population 2
    T3B: length of bottleneck in population 3
    m12: Rate of migration into 1 from 2
    m13: Rate of migration into 1 from 3
    m21: Rate of migration into 2 from 1
    m23: Rate of migration into 2 from 3
    m31: Rate of migration into 3 from 1
    m32: Rate of migration into 3 from 2
    """

    #create grid
    xx = dadi.Numerics.default_grid(pts)

    #phi for equilibrium ancestral population
    phi=PhiManip.phi_1D(xx)

    # Ancestral population growth
    phi = Integration.one_pop(phi, xx, T=T1, nu=nu1F)

    #split of population 2 from 1
    phi=PhiManip.phi_1D_to_2D(xx, phi)

    #bottleneck in population 2
    phi=Integration.two_pops(phi=phi, xx=xx, T=T2B, nu1=nu1F, nu2=nu2B)

    #migration between population 2 and 1
    phi=Integration.two_pops(phi, xx, T=T2, nu1=nu1F, nu2=nu2F, m12=m12, m21=m21)

    #split of population 3 from 1
    phi=PhiManip.phi_2D_to_3D_split_1(xx, phi)

    #bottleneck in population 3
    phi=Integration.three_pops(phi, xx, T=T3B, nu1=nu1F, nu2=nu2F, nu3=nu3B)
    
    #migration between all three populations
    phi=Integration.three_pops(phi, xx, T=T3, nu1=nu1F, nu2=nu2F, nu3=nu3F, m12=m12, m13=m13, m21=m21, m23=m23, m31=m31, m32=m32)

    #calculate spectrum
    sfs = dadi.Spectrum.from_phi_3D_direct(nx=n1, ny=n2, nz=n3, xx=xx, yy=xx, zz=xx, phi=phi)

    return sfs
