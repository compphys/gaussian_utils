#!/usr/bin/env python

from numpy import *
from atomic_data import label,mass

def get_gaussian_results(filename):
    """Reads information from Gaussian formatted checkpoint files."""

    input = open(filename,'r')
    
    s=input.readline()
    while s:
        if s.find('Number of atoms')!=-1:
            natoms=int(s.split()[4])
            N=3*natoms
        if s.find('Atomic numbers')!=-1:
            atom_l=[]
            mass_l=[]
            nrows=int( (natoms+5)/6 )
            for irow in range(nrows):
                for Z in input.readline().split():
                    atom_l.append(label(int(Z)))
                    mass_l.append(mass(label(int(Z))))
        if s.find('Current cartesian coordinates')!=-1:
            cart_l=[]
            nrows=int( (N+4)/5 )
            for irow in range(nrows):
                cart_l+=map(float,input.readline().split())
            cart_l=array(cart_l)
        if s.find('Cartesian Force Constants')!=-1:
            h=zeros((N,N))
            r=0; c=0
            nrows=int( (N*(N+1)/2+4)/5 )
            for irow in range(nrows):
                hessian_elements=input.readline().split()
                for j in hessian_elements:
                    if c<r:
                        h[r,c]=float(j)
                        h[c,r]=h[r,c]
                        c+=1
                    elif c==r:
                        h[r,r]=float(j)
                        c=0;r+=1
        s=input.readline()

    return h, atom_l, mass_l, cart_l
