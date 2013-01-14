#!/usr/bin/env python

from numpy import *

def compute_normal_coord(h,mass_l):
    """Compute normal coordinates from Cartesian Hessian."""

    from scipy import linalg
    from utils.conversion import amu2au

    #number of Cartesian coordinates
    N=size(h,0)
    if N==6:
        K=5 #diatomic molecule
    elif N>6:
        K=6 #polyatomic molecule
    else:
        print 'dimension of Hessian is incorrect, N=',N
        sys.exit(2)

    m=[]
    for i in mass_l:
        for j in range(3):
            m.append(i)

    G=1/sqrt(outer(m,m))
    #eigenvalues w are real and ordered (smallest first)
    #eigenvectors are columns in V
    w,V = linalg.eigh(h*G)

    #the first K modes are translation and rotation
    freq=sqrt(w[K:N]/amu2au)
    # array / vector is elementwise on rows in array
    Q=transpose(V.transpose()[K:N,:]/sqrt(m))

    return Q,freq
